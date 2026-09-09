"""
module_design_skill.py
=======================

Lets the chat model (1) identify Nextflow modules that a pipeline is
*missing*, (2) design them as structured specs, and (3) refine those
specs using feedback from the user -- all *before* any code is written
to disk. Actual implementation (emitting .nf files into the pipeline) is
deliberately a separate, later step; nothing here modifies the pipeline.

Why "missing" has a concrete meaning
------------------------------------
Two signals are detected statically, so the model reasons from ground
truth instead of guessing:

  1. called-but-undefined -- a process is invoked inside a `workflow {}`
     block but is not defined in any .nf file and is not imported from a
     module file that actually exists.
  2. broken include       -- `include { FOO } from './modules/foo'` whose
     target file does not exist, so FOO must be created.

Anything beyond that (does this pipeline *conceptually* need a step it
doesn't have?) is left to the model's judgement; the context tools feed
it enough call-site information to make that judgement well.

Tools exposed
-------------
    identify_missing_modules(path)
        Static scan of a pipeline dir/file: reports called-but-undefined
        processes and broken includes -- the candidate modules to create.

    gather_module_context(path, process_name)
        For one candidate: every call site, inferred inputs (the channel
        expressions passed in), upstream producers, downstream consumers,
        and a couple of existing processes as house-style templates.

    save_module_design(design)
        Persist a structured design spec (created by the model) to disk
        as JSON so it survives across chat turns / feedback rounds.

    update_module_design(name, changes, feedback="")
        Apply the user's feedback: overwrite the given fields, log the
        feedback + prior version in the design's revision history.

    get_module_design(name) / list_module_designs()
        Reload stored designs (e.g. after the user replies with feedback).

    render_module_design(name)
        Human-readable review view, including a *draft* preview of the
        Nextflow process the spec would produce. For review only -- this
        does not write anything into the pipeline.

    record_design_feedback(name, feedback, approved=False)
        Log feedback (or an approval) against a design without otherwise
        changing it.

Security / storage
------------------
Pipeline reads go through fs_skill's ALLOWED_ROOTS check (_resolve_safe).
Designs are written to a separate, writable output directory
(MODULE_DESIGN_DIR, default ./module_designs); design names are
sanitised to a safe filename so a name can never escape that directory.
"""

import json
import os
import re
import time
from pathlib import Path
from typing import Optional

from langchain_core.tools import tool

from .fs_skill import _resolve_safe
from .nextflow_skill import _extract_processes, _extract_includes

# --- Configuration -----------------------------------------------------------

MODULE_DESIGN_DIR = Path(
    os.getenv("MODULE_DESIGN_DIR", "./module_designs")
).expanduser()

MAX_NF_FILES = int(os.getenv("MD_MAX_NF_FILES", "300"))

# Directory names we never descend into when scanning a pipeline tree.
_SKIP_DIRS = {
    ".git", ".nextflow", "work", "__pycache__", "node_modules",
    ".cache", "results", "output", "outputs", "logs",
    ".singularity", ".conda", "test-datasets",
}

# Names that look like calls in a workflow body but are Nextflow/Groovy
# built-ins, channel factories, or operators -- never user processes.
NEXTFLOW_BUILTINS = {
    "Channel", "channel", "file", "files", "tuple", "path", "val", "env",
    "stdin", "stdout", "params", "workflow", "process", "value", "of",
    "empty", "create", "watchPath", "fromPath", "fromFilePairs", "fromSRA",
    "fromList", "map", "filter", "flatMap", "flatten", "collect",
    "collectFile", "reduce", "groupTuple", "groupBy", "join", "combine",
    "cross", "mix", "concat", "first", "last", "take", "unique", "distinct",
    "count", "min", "max", "sum", "toList", "toSortedList", "toInteger",
    "set", "view", "ifEmpty", "branch", "multiMap", "splitCsv", "splitText",
    "splitFasta", "splitFastq", "transpose", "dump", "println", "print",
    "error", "exit", "sleep", "tap", "buffer", "collate", "until", "spread",
    "phase", "merge", "into", "separate", "subscribe", "close", "randomSample",
    "if", "else", "for", "while", "return", "def", "new", "assert", "throw",
    "try", "catch", "log", "emit", "main", "as", "in", "size", "name",
}


# --- Static analysis helpers --------------------------------------------------

def _collect_nf_files(base: Path) -> list[Path]:
    """All .nf files under `base`, skipping noisy/output directories."""
    out = []
    for p in sorted(base.rglob("*.nf")):
        if any(part in _SKIP_DIRS for part in p.parts):
            continue
        if p.is_file():
            out.append(p)
        if len(out) >= MAX_NF_FILES:
            break
    return out


def _strip_scripts_and_strings(text: str) -> str:
    """Remove triple-quoted script blocks and quoted strings so their
    contents (e.g. a bash `foo(` inside a script) aren't mistaken for
    Nextflow process calls."""
    text = re.sub(r'""".*?"""', " ", text, flags=re.DOTALL)
    text = re.sub(r"'''.*?'''", " ", text, flags=re.DOTALL)
    text = re.sub(r'"[^"\n]*"', " ", text)
    text = re.sub(r"'[^'\n]*'", " ", text)
    return text


def _workflow_bodies(content: str) -> list[tuple[str, str]]:
    """Yield (workflow_name, full_body) for each `workflow [NAME] { ... }`
    using brace matching so nested closures are handled."""
    content = re.sub(r'""".*?"""', " ", content, flags=re.DOTALL)
    content = re.sub(r"'''.*?'''", " ", content, flags=re.DOTALL)
    bodies = []
    for m in re.finditer(r"workflow\s*(\w*)\s*\{", content):
        name = m.group(1) or "(entry)"
        start = m.end()
        depth, pos = 1, start
        while pos < len(content) and depth > 0:
            c = content[pos]
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
            pos += 1
        bodies.append((name, content[start:pos - 1]))
    return bodies


def _extract_calls(body: str) -> list[str]:
    """Process-like call tokens in a workflow body: `NAME(` not preceded
    by a `.` (excludes `chan.map(`) and not a known built-in/operator."""
    body = _strip_scripts_and_strings(body)
    calls = re.findall(r"(?<![.\w])([A-Za-z_]\w*)\s*\(", body)
    seen, out = set(), []
    for c in calls:
        if c in NEXTFLOW_BUILTINS or c in seen:
            continue
        seen.add(c)
        out.append(c)
    return out


def _resolve_include_source(including_file: Path, source: str) -> Optional[Path]:
    """Resolve an include's `from '<source>'` relative to the including
    file. Nextflow allows the .nf extension to be omitted."""
    src = source.strip()
    cand = including_file.parent / src
    for p in (cand, cand.with_suffix(".nf"), Path(str(cand) + ".nf")):
        try:
            if p.exists() and p.is_file():
                return p
        except OSError:
            continue
    return None


def _split_top_level(s: str) -> list[str]:
    """Split on commas that are not inside (), [], or {}."""
    parts, depth, cur = [], 0, ""
    for ch in s:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        parts.append(cur.strip())
    return parts


def _balanced_args(body: str, open_pos: int) -> tuple[str, int]:
    """Given the index just after a `(`, return (arg_text, end_index)."""
    depth, pos = 1, open_pos
    while pos < len(body) and depth > 0:
        if body[pos] == "(":
            depth += 1
        elif body[pos] == ")":
            depth -= 1
        pos += 1
    return body[open_pos:pos - 1], pos


def _analyze_pipeline(base: Path) -> dict:
    """Core static scan. Returns defined/imported/called sets plus the
    computed missing modules and broken includes."""
    files = _collect_nf_files(base)
    defined: set[str] = set()
    imported_resolvable: set[str] = set()
    workflow_names: set[str] = set()
    broken_includes: list[dict] = []
    called: dict[str, list[dict]] = {}

    for f in files:
        try:
            content = f.read_text(encoding="utf-8", errors="replace")
        except Exception:
            continue

        for p in _extract_processes(content):
            defined.add(p["name"])

        for inc in _extract_includes(content):
            names = [(it.get("alias") or it["name"]) for it in inc["items"]]
            if _resolve_include_source(f, inc["source"]):
                imported_resolvable.update(names)
            else:
                broken_includes.append({
                    "file": str(f),
                    "source": inc["source"],
                    "items": names,
                })

        for wname, body in _workflow_bodies(content):
            if wname != "(entry)":
                workflow_names.add(wname)
            for c in _extract_calls(body):
                called.setdefault(c, []).append(
                    {"file": f.name, "workflow": wname}
                )

    accounted = defined | imported_resolvable | workflow_names
    missing_called = {n: s for n, s in called.items() if n not in accounted}

    broken_items: set[str] = set()
    for b in broken_includes:
        broken_items.update(b["items"])

    return {
        "files_scanned": len(files),
        "defined": sorted(defined),
        "imported_resolvable": sorted(imported_resolvable),
        "workflow_names": sorted(workflow_names),
        "missing_called": missing_called,
        "broken_includes": broken_includes,
        "broken_items": sorted(broken_items),
        "candidates": sorted(set(missing_called) | broken_items),
    }


# --- Detection / context tools -----------------------------------------------

@tool
def identify_missing_modules(path: str) -> str:
    """
    Statically scan a Nextflow pipeline directory (or a single .nf file's
    directory) and report which modules appear to be MISSING and therefore
    are candidates to be designed and created.

    Two concrete signals are reported:
      • called-but-undefined: a process invoked in a workflow block that
        is not defined anywhere and not imported from an existing module.
      • broken includes: `include { X } from '<file>'` where <file> does
        not exist -- X must be created.

    This is READ-ONLY and the correct FIRST step when the user asks to
    find, design, or add missing modules. `path` must be inside an
    allow-listed root. After this, call gather_module_context for each
    candidate before drafting a design.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"

    base = target if target.is_dir() else target.parent
    if not base.exists():
        return f"Error: '{base}' does not exist."

    info = _analyze_pipeline(base)

    lines = [f"=== Missing-module scan: {base.name} ===",
             f"(scanned {info['files_scanned']} .nf file(s))\n"]

    if not info["candidates"]:
        lines.append("✅ No missing modules detected.")
        lines.append(
            f"  Processes defined: {len(info['defined'])}, "
            f"imported: {len(info['imported_resolvable'])}."
        )
        lines.append(
            "\nNote: this only catches processes that are *called but "
            "absent*. If the user believes a conceptual step is missing "
            "(e.g. no QC stage at all), reason about that from the "
            "pipeline's processes instead."
        )
        return "\n".join(lines)

    lines.append(f"⚠️  {len(info['candidates'])} candidate module(s) to design:\n")

    if info["missing_called"]:
        lines.append("CALLED BUT NOT DEFINED / NOT IMPORTED:")
        for name, sites in info["missing_called"].items():
            where = ", ".join(
                f"{s['file']}:workflow {s['workflow']}" for s in sites[:3]
            )
            lines.append(f"  • {name}   (called in {where})")
        lines.append("")

    if info["broken_includes"]:
        lines.append("BROKEN INCLUDES (imported from a file that doesn't exist):")
        for b in info["broken_includes"]:
            items = ", ".join(b["items"])
            lines.append(
                f"  • {items}   ← include from '{b['source']}' "
                f"(in {Path(b['file']).name})"
            )
        lines.append("")

    lines.append("NEXT STEPS:")
    lines.append("  For each candidate, call:")
    lines.append(f"    gather_module_context('{path}', '<PROCESS_NAME>')")
    lines.append("  then draft a spec and save it with save_module_design(...).")
    lines.append(
        "  Present each design to the user and ask for feedback on the "
        "interface (inputs/outputs), resources, and container/conda "
        "BEFORE implementing anything."
    )
    return "\n".join(lines)


@tool
def gather_module_context(path: str, process_name: str) -> str:
    """
    Collect everything known about a single missing process so a good
    interface can be designed for it. Reports, from the pipeline's
    workflow blocks:
      • every call site and the exact call expression
      • the arguments passed in (these become the process INPUTS)
      • upstream processes whose `.out` feeds this one (depends_on)
      • downstream processes that consume this one's `.out` (feeds_into)
      • one or two existing processes as house-style templates

    Call this after identify_missing_modules and before drafting a design
    with save_module_design. `path` must be inside an allow-listed root.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"

    base = target if target.is_dir() else target.parent
    files = _collect_nf_files(base)

    call_sites: list[dict] = []
    upstream: set[str] = set()
    downstream: set[str] = set()

    for f in files:
        try:
            content = f.read_text(encoding="utf-8", errors="replace")
        except Exception:
            continue
        for wname, body in _workflow_bodies(content):
            # Call sites of the target itself + inferred inputs / upstream.
            for m in re.finditer(
                r"(?<![.\w])" + re.escape(process_name) + r"\s*\(", body
            ):
                args_raw, _ = _balanced_args(body, m.end())
                args = _split_top_level(args_raw)
                for a in args:
                    um = re.match(r"([A-Za-z_]\w*)\s*\.\s*out", a)
                    if um:
                        upstream.add(um.group(1))
                call_sites.append({
                    "file": f.name,
                    "workflow": wname,
                    "call": f"{process_name}({args_raw.strip()})",
                    "args": args,
                })
            # Downstream consumers: any other call whose args mention NAME.out
            for m in re.finditer(r"(?<![.\w])([A-Za-z_]\w*)\s*\(", body):
                callee = m.group(1)
                if callee == process_name or callee in NEXTFLOW_BUILTINS:
                    continue
                callargs, _ = _balanced_args(body, m.end())
                if re.search(
                    r"(?<![.\w])" + re.escape(process_name) + r"\s*\.\s*out",
                    callargs,
                ):
                    downstream.add(callee)

    if not call_sites:
        return (
            f"No call sites for '{process_name}' were found in workflow "
            f"blocks under '{base}'. Double-check the name (case-sensitive), "
            f"or run identify_missing_modules first."
        )

    # A couple of existing defined processes as style templates.
    templates = []
    for f in files:
        try:
            content = f.read_text(encoding="utf-8", errors="replace")
        except Exception:
            continue
        for m in re.finditer(r"process\s+(\w+)\s*\{", content):
            pname = m.group(1)
            body, _ = _balanced_args_brace(content, m.end())
            inp = _first_section(body, "input")
            out = _first_section(body, "output")
            templates.append({"name": pname, "input": inp, "output": out})
            if len(templates) >= 2:
                break
        if len(templates) >= 2:
            break

    lines = [f"=== Context for missing process: {process_name} ===\n"]

    n_inputs = max((len(cs["args"]) for cs in call_sites), default=0)
    lines.append(f"CALL SITES ({len(call_sites)}):")
    for cs in call_sites[:8]:
        lines.append(f"  • {cs['call']}   [{cs['file']}: workflow {cs['workflow']}]")
    lines.append("")
    lines.append(
        f"INFERRED INPUTS: {n_inputs} channel argument(s) at the widest "
        f"call site."
    )
    if call_sites:
        for i, a in enumerate(call_sites[0]["args"], 1):
            lines.append(f"    input {i}: {a}")
    lines.append("")
    lines.append(f"UPSTREAM (produces this process's inputs): "
                 f"{', '.join(sorted(upstream)) or 'none detected'}")
    lines.append(f"DOWNSTREAM (consumes this process's .out): "
                 f"{', '.join(sorted(downstream)) or 'none detected'}")

    if templates:
        lines.append("\nEXISTING PROCESSES (house-style templates):")
        for t in templates:
            lines.append(f"  ── {t['name']} ──")
            if t["input"]:
                lines.append(f"     input:  {t['input'][:120]}")
            if t["output"]:
                lines.append(f"     output: {t['output'][:120]}")

    lines.append(
        "\nUse this to draft a spec (name, purpose, inputs, outputs, "
        "directives, script) and save it with save_module_design(...). "
        "The downstream consumer tells you what the OUTPUT must look like."
    )
    return "\n".join(lines)


def _balanced_args_brace(content: str, open_pos: int) -> tuple[str, int]:
    """Like _balanced_args but for `{ }` (process bodies)."""
    depth, pos = 1, open_pos
    while pos < len(content) and depth > 0:
        if content[pos] == "{":
            depth += 1
        elif content[pos] == "}":
            depth -= 1
        pos += 1
    return content[open_pos:pos - 1], pos


def _first_section(body: str, section: str) -> str:
    """First line of an input:/output: section within a process body."""
    pat = re.compile(
        rf"{section}\s*:\s*\n?\s*([^\n]+)", re.IGNORECASE
    )
    m = pat.search(body)
    return m.group(1).strip() if m else ""


# --- Design store -------------------------------------------------------------

_REQUIRED_FIELDS = ("name", "purpose")


def _safe_design_name(name: str) -> str:
    """Sanitise a design name into a safe filename stem. Prevents any
    path traversal -- the name can only ever land inside MODULE_DESIGN_DIR."""
    stem = re.sub(r"[^A-Za-z0-9_]+", "_", str(name)).strip("_")
    if not stem:
        raise ValueError("design name is empty after sanitising.")
    return stem


def _design_path(name: str) -> Path:
    MODULE_DESIGN_DIR.mkdir(parents=True, exist_ok=True)
    return MODULE_DESIGN_DIR / f"{_safe_design_name(name)}.json"


def _validate_design(design: dict) -> tuple[bool, str]:
    if not isinstance(design, dict):
        return False, "design must be a JSON object."
    for field in _REQUIRED_FIELDS:
        if not design.get(field):
            return False, f"missing required field: '{field}'."
    for list_field in ("inputs", "outputs"):
        val = design.get(list_field)
        if val is not None and not isinstance(val, list):
            return False, f"'{list_field}' must be a list if provided."
    return True, ""


def _render_nextflow_process(d: dict) -> str:
    """Render a design spec into DRAFT DSL2 process text for review.
    This is a preview only -- it is NOT written into the pipeline."""
    name = d.get("name", "UNNAMED")
    lines = [f"process {name} {{"]

    directives = d.get("directives") or {}
    # Emit directives in a conventional order.
    order = ["tag", "label", "cpus", "memory", "time", "container",
             "conda", "publishDir", "errorStrategy", "maxRetries"]
    for key in order:
        if key in directives and directives[key] not in (None, ""):
            val = directives[key]
            if key in ("cpus", "maxRetries"):
                lines.append(f"    {key} {val}")
            else:
                lines.append(f"    {key} '{val}'")
    for key, val in directives.items():
        if key not in order and val not in (None, ""):
            lines.append(f"    {key} {val}")

    def _decl(item):
        if isinstance(item, str):
            return item
        if item.get("declaration"):
            return item["declaration"]
        typ = item.get("type", "path")
        nm = item.get("name") or item.get("pattern") or "value"
        emit = item.get("emit")
        base = f"{typ} {nm}" if typ in ("val", "env") else f'{typ} "{nm}"' \
            if typ == "path" and any(ch in nm for ch in "*.?") else f"{typ} {nm}"
        if emit:
            base += f", emit: {emit}"
        return base

    inputs = d.get("inputs") or []
    if inputs:
        lines.append("")
        lines.append("    input:")
        for item in inputs:
            lines.append(f"    {_decl(item)}")

    outputs = d.get("outputs") or []
    if outputs:
        lines.append("")
        lines.append("    output:")
        for item in outputs:
            lines.append(f"    {_decl(item)}")

    when = d.get("when")
    if when:
        lines.append("")
        lines.append(f"    when:\n    {when}")

    lines.append("")
    lines.append("    script:")
    script = d.get("script") or "// TODO: fill in command"
    lines.append('    """')
    for ln in str(script).splitlines() or [""]:
        lines.append(f"    {ln}")
    lines.append('    """')
    lines.append("}")
    return "\n".join(lines)


@tool
def save_module_design(design: dict) -> str:
    """
    Persist a module DESIGN SPEC (which you, the assistant, author) to
    disk so it survives across chat turns and feedback rounds. Call this
    once you've drafted a design from gather_module_context output.

    `design` is a JSON object with:
      name        (str, required) e.g. "SORT_BAM"
      purpose     (str, required) one-line description
      tool        (str) underlying CLI tool, e.g. "samtools"
      inputs      (list) each: {"type": "path|val|tuple", "name": "...",
                                "description": "..."}  OR
                          {"declaration": "tuple val(id), path(bam)"}
      outputs     (list) each: {"type": "path", "name": "*.bam",
                                "emit": "bam", "description": "..."} OR
                          {"declaration": "path '*.sorted.bam', emit: bam"}
      directives  (obj)  {"cpus": 4, "memory": "8 GB",
                          "container": "...", "publishDir": "..."}
      script      (str)  the command(s) for the script: block
      depends_on  (list) upstream process names
      feeds_into  (list) downstream process names
      notes       (str)  rationale / open questions for the user
      status      (str)  defaults to "draft"

    Saving does NOT write any .nf file into the pipeline -- it only
    records the design for review. After saving, show the user the design
    (render_module_design) and ask for feedback before implementing.
    """
    ok, err = _validate_design(design)
    if not ok:
        return f"Error: {err}"

    design = dict(design)
    design.setdefault("status", "draft")
    name = design["name"]

    path = _design_path(name)
    existing = None
    if path.exists():
        try:
            existing = json.loads(path.read_text())
        except Exception:
            existing = None

    if existing:
        history = existing.get("_history", [])
        snapshot = {k: v for k, v in existing.items() if k != "_history"}
        history.append({"at": time.strftime("%Y-%m-%d %H:%M:%S"),
                        "event": "overwritten", "snapshot": snapshot})
        design["_history"] = history
        design["revision"] = existing.get("revision", 1) + 1
    else:
        design["_history"] = []
        design["revision"] = 1

    path.write_text(json.dumps(design, indent=2))
    preview = _render_nextflow_process(design)
    return (
        f"Saved design '{name}' (revision {design['revision']}, "
        f"status: {design['status']}) → {path}\n\n"
        f"DRAFT preview (for review only, not written to the pipeline):\n"
        f"{'-'*50}\n{preview}\n{'-'*50}\n\n"
        f"Now present this to the user and ask whether the inputs/outputs, "
        f"resources, and container/conda look right before implementing."
    )


@tool
def update_module_design(name: str, changes: dict, feedback: str = "") -> str:
    """
    Revise a saved design using the user's FEEDBACK. Overwrites only the
    fields present in `changes` (e.g. {"outputs": [...], "directives":
    {...}}) and appends the prior version plus the feedback text to the
    design's revision history, so the reasoning trail is preserved.

    Use this whenever the user responds to a proposed design with
    changes -- e.g. "make the output emit a sorted+indexed pair", "use
    the biocontainers samtools image", "it should take the reference as a
    second input". Pass their comment in `feedback`.
    """
    path = _design_path(name)
    if not path.exists():
        return (
            f"No saved design named '{name}'. Use list_module_designs() to "
            f"see existing designs, or save_module_design(...) to create it."
        )
    try:
        design = json.loads(path.read_text())
    except Exception as e:
        return f"Error reading design '{name}': {e}"

    if not isinstance(changes, dict):
        return "Error: 'changes' must be a JSON object of fields to update."

    history = design.get("_history", [])
    snapshot = {k: v for k, v in design.items() if k != "_history"}
    history.append({
        "at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "event": "feedback",
        "feedback": feedback or "(none provided)",
        "changed_fields": sorted(changes.keys()),
        "snapshot": snapshot,
    })

    for k, v in changes.items():
        if k in ("_history", "revision"):
            continue
        design[k] = v

    design["_history"] = history
    design["revision"] = design.get("revision", 1) + 1
    design.setdefault("status", "draft")
    if design["status"] == "approved" and changes:
        design["status"] = "draft"  # changed after approval → back to draft

    ok, err = _validate_design(design)
    if not ok:
        return f"Error: update would make the design invalid ({err}). Not saved."

    path.write_text(json.dumps(design, indent=2))
    preview = _render_nextflow_process(design)
    return (
        f"Updated '{name}' to revision {design['revision']} "
        f"(applied feedback: {feedback or 'n/a'}).\n\n"
        f"DRAFT preview (review only):\n{'-'*50}\n{preview}\n{'-'*50}\n\n"
        f"Show this back to the user and confirm it now matches their intent."
    )


@tool
def get_module_design(name: str) -> str:
    """
    Return the stored JSON for a saved design. Use this to reload a design
    in a later turn (e.g. after the user replies with feedback) before
    updating it. Revision history is omitted for brevity.
    """
    path = _design_path(name)
    if not path.exists():
        return f"No saved design named '{name}'."
    try:
        design = json.loads(path.read_text())
    except Exception as e:
        return f"Error reading design '{name}': {e}"
    design.pop("_history", None)
    return json.dumps(design, indent=2)


@tool
def list_module_designs() -> str:
    """
    List all saved module designs with their status, revision, and
    purpose. Use this to see what has already been designed in this or a
    previous session before proposing new work.
    """
    if not MODULE_DESIGN_DIR.exists():
        return "No designs saved yet."
    files = sorted(MODULE_DESIGN_DIR.glob("*.json"))
    if not files:
        return "No designs saved yet."
    lines = [f"Saved module designs ({len(files)}):\n"]
    for f in files:
        try:
            d = json.loads(f.read_text())
        except Exception:
            lines.append(f"  • {f.stem}  (unreadable)")
            continue
        lines.append(
            f"  • {d.get('name', f.stem)}  "
            f"[{d.get('status', '?')}, rev {d.get('revision', '?')}]"
        )
        if d.get("purpose"):
            lines.append(f"      {d['purpose']}")
    return "\n".join(lines)


@tool
def render_module_design(name: str) -> str:
    """
    Produce a human-readable review of a saved design, including a DRAFT
    preview of the Nextflow process it would generate. Use this to show
    the user a design and invite feedback. This does NOT write anything
    into the pipeline -- it is for review only.
    """
    path = _design_path(name)
    if not path.exists():
        return f"No saved design named '{name}'."
    try:
        d = json.loads(path.read_text())
    except Exception as e:
        return f"Error reading design '{name}': {e}"

    lines = [f"=== Module design: {d.get('name', name)} "
             f"(status: {d.get('status', '?')}, rev {d.get('revision', '?')}) ==="]
    if d.get("purpose"):
        lines.append(f"\nPurpose: {d['purpose']}")
    if d.get("tool"):
        lines.append(f"Tool:    {d['tool']}")
    if d.get("depends_on"):
        lines.append(f"Upstream:   {', '.join(d['depends_on'])}")
    if d.get("feeds_into"):
        lines.append(f"Downstream: {', '.join(d['feeds_into'])}")

    if d.get("inputs"):
        lines.append("\nInputs:")
        for it in d["inputs"]:
            if isinstance(it, str):
                lines.append(f"  • {it}")
            else:
                desc = f" — {it['description']}" if it.get("description") else ""
                lines.append(f"  • {it.get('declaration') or it.get('type','')} "
                             f"{it.get('name', it.get('pattern',''))}{desc}")
    if d.get("outputs"):
        lines.append("\nOutputs:")
        for it in d["outputs"]:
            if isinstance(it, str):
                lines.append(f"  • {it}")
            else:
                desc = f" — {it['description']}" if it.get("description") else ""
                emit = f" (emit: {it['emit']})" if it.get("emit") else ""
                lines.append(f"  • {it.get('declaration') or it.get('type','')} "
                             f"{it.get('name', it.get('pattern',''))}{emit}{desc}")
    if d.get("notes"):
        lines.append(f"\nNotes / open questions: {d['notes']}")

    lines.append("\nDRAFT Nextflow process (review only — not yet implemented):")
    lines.append("-" * 50)
    lines.append(_render_nextflow_process(d))
    lines.append("-" * 50)
    lines.append(
        "\nAsk the user: do the inputs/outputs, resource directives, and "
        "container/conda match their intent? Apply changes with "
        "update_module_design(...)."
    )
    return "\n".join(lines)


@tool
def record_design_feedback(name: str, feedback: str, approved: bool = False) -> str:
    """
    Log the user's feedback (or an approval) against a design WITHOUT
    otherwise changing its fields. Use this when the user comments but you
    aren't changing the spec yet, or when they approve it as-is (pass
    approved=True to mark status "approved"). For feedback that changes
    fields, use update_module_design instead.
    """
    path = _design_path(name)
    if not path.exists():
        return f"No saved design named '{name}'."
    try:
        design = json.loads(path.read_text())
    except Exception as e:
        return f"Error reading design '{name}': {e}"

    history = design.get("_history", [])
    history.append({
        "at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "event": "approval" if approved else "comment",
        "feedback": feedback,
    })
    design["_history"] = history
    if approved:
        design["status"] = "approved"
    path.write_text(json.dumps(design, indent=2))

    if approved:
        return (
            f"Marked '{name}' as APPROVED and logged the feedback. "
            f"It is now ready for the (separate) implementation step; "
            f"nothing has been written to the pipeline yet."
        )
    return f"Logged feedback on '{name}' (status unchanged: {design.get('status')})."
