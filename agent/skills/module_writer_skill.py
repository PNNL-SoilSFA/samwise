"""
module_writer_skill.py
======================

Writes designed Nextflow modules (from module_design_skill.py) to disk
as proper .nf files organized in a standard module structure.

This skill:
  1. Reads module designs from the design store
  2. Renders them as DSL2-compliant Nextflow process definitions
  3. Writes them to a user-specified output directory
  4. Creates appropriate directory structure (following nf-core conventions)
  5. Generates supporting files (meta.yml, includes)

The output follows nf-core local module conventions:
    <output_dir>/
    ├── modules/
    │   └── local/
    │       ├── <module_name>/
    │       │   ├── main.nf        # Process definition
    │       │   └── meta.yml       # Module metadata (optional)
    │       └── ...
    └── modules.nf                 # Optional: includes all modules

References:
  - nf-core module structure: https://nf-co.re/docs/contributing/modules
  - Nextflow DSL2 modules: https://www.nextflow.io/docs/latest/dsl2.html#modules
  - nf-core local modules: https://nf-co.re/docs/contributing/modules#local-modules

Configuration (via .env):
    MODULE_WRITE_ROOTS    - Colon-separated list of directories where modules
                           can be written (default: ./output:./generated)
    MODULE_BACKUP_ENABLED - Whether to backup before overwriting (default: true)

Security:
  - Write operations are restricted to MODULE_WRITE_ROOTS
  - Read operations for designs use module_design_skill's existing store
  - Overwrites require explicit confirmation or create backups
"""

import json
import os
import re
import shutil
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

from langchain_core.tools import tool

# Import from module_design_skill to access designs
from .module_design_skill import (
    MODULE_DESIGN_DIR,
    _safe_design_name,
    _render_nextflow_process,
)

# --- Configuration -----------------------------------------------------------

_ROOTS_SEP = ";" if os.name == "nt" else ":"

# Directories where the agent is allowed to WRITE module files.
# This is separate from ALLOWED_ROOTS (read) for security.
_DEFAULT_WRITE_ROOTS = [
    "./output",
    "./generated",
    "./pipeline_modules",
]

_raw_write_roots = os.getenv("MODULE_WRITE_ROOTS", "")
MODULE_WRITE_ROOTS = [
    Path(p).expanduser().resolve()
    for p in _raw_write_roots.split(_ROOTS_SEP) if p.strip()
] or [Path(p).expanduser().resolve() for p in _DEFAULT_WRITE_ROOTS]

MODULE_BACKUP_ENABLED = os.getenv("MODULE_BACKUP_ENABLED", "true").lower() != "false"

# nf-core style module structure
MODULE_SUBDIR = "modules/local"


# --- Path Safety -------------------------------------------------------------

def _resolve_write_safe(path_str: str) -> Path:
    """
    Resolve a path and ensure it's inside one of MODULE_WRITE_ROOTS.
    
    Similar to fs_skill's _resolve_safe but for write operations.
    Raises ValueError if the path is outside allowed write directories.
    
    Reference: Defense in depth principle for file system access
    https://cheatsheetseries.owasp.org/cheatsheets/File_Upload_Cheat_Sheet.html
    """
    candidate = Path(path_str).expanduser().resolve()
    
    # Check if path is under an allowed write root
    for root in MODULE_WRITE_ROOTS:
        try:
            candidate.relative_to(root)
            return candidate
        except ValueError:
            continue
    
    # Also allow if the path itself would become a write root
    # (for creating new output directories)
    for root in MODULE_WRITE_ROOTS:
        try:
            root.relative_to(candidate)
            # candidate is a parent of a write root - not allowed
            continue
        except ValueError:
            pass
        
        # Check if candidate is a subdirectory we'd create under a write root
        # by checking if any write root is a prefix
        if str(candidate).startswith(str(root)):
            return candidate
    
    allowed = ", ".join(str(r) for r in MODULE_WRITE_ROOTS)
    raise ValueError(
        f"'{candidate}' is outside allowed write directories ({allowed}). "
        f"Set MODULE_WRITE_ROOTS in .env to add write permissions for this path."
    )


def _ensure_write_root_exists(path: Path) -> Path:
    """
    Ensure at least one write root exists, creating the first one if needed.
    Returns the path unchanged after validation.
    """
    # Find which write root this path would be under
    for root in MODULE_WRITE_ROOTS:
        if str(path).startswith(str(root)) or path == root:
            root.mkdir(parents=True, exist_ok=True)
            return path
    
    # If we get here, path wasn't validated - shouldn't happen after _resolve_write_safe
    raise ValueError(f"Cannot determine write root for '{path}'")


# --- Design Access -----------------------------------------------------------

def _load_design(name: str) -> Optional[dict]:
    """
    Load a design from the module_design_skill's store.
    
    Returns None if not found.
    """
    design_path = MODULE_DESIGN_DIR / f"{_safe_design_name(name)}.json"
    if not design_path.exists():
        return None
    
    try:
        return json.loads(design_path.read_text())
    except Exception:
        return None


def _list_designs(status_filter: Optional[str] = None) -> list[dict]:
    """
    List all saved designs, optionally filtering by status.
    
    Args:
        status_filter: If provided, only return designs with this status
                      (e.g., "approved", "draft")
    """
    if not MODULE_DESIGN_DIR.exists():
        return []
    
    designs = []
    for f in MODULE_DESIGN_DIR.glob("*.json"):
        try:
            d = json.loads(f.read_text())
            if status_filter is None or d.get("status") == status_filter:
                designs.append(d)
        except Exception:
            continue
    
    return designs


# --- File Generation ---------------------------------------------------------

def _generate_module_header(design: dict) -> str:
    """
    Generate the header comment for a module file.
    
    Follows nf-core conventions for module documentation.
    Reference: https://nf-co.re/docs/contributing/modules#documentation
    """
    lines = [
        "//",
        "// " + design.get("name", "UNNAMED"),
        "//",
    ]
    
    if design.get("purpose"):
        lines.append("// " + design["purpose"])
        lines.append("//")
    
    if design.get("tool"):
        lines.append("// Tool: " + design["tool"])
    
    if design.get("depends_on"):
        lines.append("// Upstream: " + ", ".join(design["depends_on"]))
    
    if design.get("feeds_into"):
        lines.append("// Downstream: " + ", ".join(design["feeds_into"]))
    
    lines.append("//")
    lines.append("// Generated: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    lines.append("// Design revision: " + str(design.get("revision", "unknown")))
    lines.append("//")
    lines.append("")
    
    return "\n".join(lines)


def _generate_module_file(design: dict) -> str:
    """
    Generate complete module file content including header and process.
    """
    header = _generate_module_header(design)
    process = _render_nextflow_process(design)
    
    return header + process + "\n"


def _generate_meta_yml(design: dict) -> str:
    """
    Generate a meta.yml file for the module (nf-core convention).
    
    Reference: https://nf-co.re/docs/contributing/modules#meta-yml
    """
    lines = [
        "name: " + design.get("name", "unnamed").lower(),
        "description: " + design.get("purpose", "No description"),
    ]
    
    if design.get("tool"):
        lines.append("keywords:")
        lines.append("  - " + design["tool"])
    
    lines.append("tools:")
    if design.get("tool"):
        lines.append("  - " + design["tool"] + ":")
        lines.append("      description: " + design.get("purpose", "Tool used by this module"))
        lines.append("      homepage: # Add tool homepage")
        lines.append("      documentation: # Add documentation URL")
    
    lines.append("")
    lines.append("input:")
    for inp in design.get("inputs", []):
        if isinstance(inp, str):
            lines.append("  - " + inp)
        else:
            name = inp.get("name", inp.get("pattern", "input"))
            desc = inp.get("description", "Input file/value")
            lines.append("  - " + name + ":")
            lines.append("      type: " + inp.get("type", "file"))
            lines.append("      description: " + desc)
    
    lines.append("")
    lines.append("output:")
    for out in design.get("outputs", []):
        if isinstance(out, str):
            lines.append("  - " + out)
        else:
            name = out.get("emit", out.get("name", out.get("pattern", "output")))
            desc = out.get("description", "Output file/value")
            lines.append("  - " + name + ":")
            lines.append("      type: " + out.get("type", "file"))
            lines.append("      description: " + desc)
    
    lines.append("")
    lines.append("authors:")
    lines.append("  - # Add author @github_handle")
    
    return "\n".join(lines)


def _generate_includes_file(module_names: list[str], output_dir: Path) -> str:
    """
    Generate a modules.nf file that includes all local modules.
    
    This makes it easy to import all modules with a single include statement.
    """
    lines = [
        "//",
        "// Local module includes",
        "// Generated: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "//",
        "",
    ]
    
    for name in sorted(module_names):
        module_path = "./" + MODULE_SUBDIR + "/" + name.lower() + "/main.nf"
        lines.append("include { " + name + " } from '" + module_path + "'")
    
    lines.append("")
    return "\n".join(lines)


def _generate_main_nf(pipeline_name: str, safe_name: str) -> str:
    """
    Generate the main.nf file content for a pipeline scaffold.
    
    Reference: Nextflow DSL2 workflow structure
    https://www.nextflow.io/docs/latest/dsl2.html
    """
    lines = [
        "#!/usr/bin/env nextflow",
        "",
        "/*",
        " * " + pipeline_name,
        " * ",
        " * A Nextflow DSL2 pipeline.",
        " * Generated: " + datetime.now().strftime("%Y-%m-%d"),
        " */",
        "",
        "nextflow.enable.dsl = 2",
        "",
        "// Import modules",
        "// include { PROCESS_NAME } from './modules/local/process_name/main.nf'",
        "",
        "// Pipeline parameters",
        "params.input = null",
        "params.outdir = './results'",
        "",
        "// Validate inputs",
        "if (!params.input) {",
        '    error "Please provide input with --input"',
        "}",
        "",
        "// Main workflow",
        "workflow {",
        "    // Channel for input",
        "    ch_input = Channel.fromPath(params.input)",
        "    ",
        "    // TODO: Add process calls here",
        "    // PROCESS_NAME(ch_input)",
        "}",
        "",
        "// On completion",
        "workflow.onComplete {",
        '    println "Pipeline completed at: $workflow.complete"',
        '    println "Execution status: ${ workflow.success ? \'OK\' : \'failed\' }"',
        "}",
    ]
    
    return "\n".join(lines)


def _generate_nextflow_config(pipeline_name: str, safe_name: str) -> str:
    """
    Generate the nextflow.config file content.
    
    Reference: Nextflow configuration documentation
    https://www.nextflow.io/docs/latest/config.html
    """
    lines = [
        "/*",
        " * " + pipeline_name + " configuration",
        " */",
        "",
        "// Pipeline metadata",
        "manifest {",
        "    name            = '" + safe_name + "'",
        "    author          = 'Your Name'",
        "    description     = '" + pipeline_name + "'",
        "    mainScript      = 'main.nf'",
        "    nextflowVersion = '>=21.10.0'",
        "    version         = '0.1.0'",
        "}",
        "",
        "// Default parameters",
        "params {",
        "    input   = null",
        "    outdir  = './results'",
        "    help    = false",
        "}",
        "",
        "// Process defaults",
        "process {",
        "    cpus   = 1",
        "    memory = '2 GB'",
        "    time   = '1h'",
        "    ",
        "    errorStrategy = 'retry'",
        "    maxRetries    = 2",
        "}",
        "",
        "// Execution profiles",
        "profiles {",
        "    standard {",
        "        process.executor = 'local'",
        "    }",
        "    ",
        "    docker {",
        "        docker.enabled = true",
        "        docker.runOptions = '-u $(id -u):$(id -g)'",
        "    }",
        "    ",
        "    singularity {",
        "        singularity.enabled = true",
        "        singularity.autoMounts = true",
        "    }",
        "    ",
        "    conda {",
        "        conda.enabled = true",
        "    }",
        "}",
        "",
        "// Capture execution info",
        "timeline {",
        "    enabled = true",
        '    file    = "${params.outdir}/pipeline_info/timeline.html"',
        "}",
        "",
        "report {",
        "    enabled = true",
        '    file    = "${params.outdir}/pipeline_info/report.html"',
        "}",
    ]
    
    return "\n".join(lines)


def _generate_readme(pipeline_name: str) -> str:
    """
    Generate the README.md file content.
    """
    lines = [
        "# " + pipeline_name,
        "",
        "A Nextflow DSL2 pipeline.",
        "",
        "## Usage",
        "",
        "```",
        "nextflow run main.nf --input <input_path> --outdir <output_dir>",
        "```",
        "",
        "## Parameters",
        "",
        "| Parameter | Description | Default |",
        "|-----------|-------------|---------|",
        "| `--input` | Input file path | (required) |",
        "| `--outdir` | Output directory | `./results` |",
        "",
        "## Profiles",
        "",
        "- `standard`: Local execution",
        "- `docker`: Run with Docker containers",
        "- `singularity`: Run with Singularity containers",
        "- `conda`: Run with Conda environments",
        "",
        "## Modules",
        "",
        "Local modules are in `modules/local/`.",
        "",
        "## Output",
        "",
        "Results are written to the directory specified by `--outdir`.",
    ]
    
    return "\n".join(lines)


# --- File Writing ------------------------------------------------------------

def _backup_file(path: Path) -> Optional[Path]:
    """
    Create a backup of an existing file before overwriting.
    
    Returns the backup path, or None if no backup was needed.
    
    Reference: Standard practice for safe file updates
    https://docs.python.org/3/library/shutil.html#shutil.copy2
    """
    if not path.exists():
        return None
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = path.with_suffix(".backup_" + timestamp + path.suffix)
    
    shutil.copy2(path, backup_path)
    return backup_path


def _write_file_atomic(path: Path, content: str, backup: bool = True) -> dict:
    """
    Write content to a file atomically with optional backup.
    
    Uses write-to-temp-then-rename pattern for atomicity.
    
    Reference: Atomic file writes in Python
    https://docs.python.org/3/library/os.html#os.replace
    
    Returns dict with write status information.
    """
    result = {
        "path": str(path),
        "backup_path": None,
        "created": not path.exists(),
        "bytes_written": 0,
    }
    
    # Create parent directories
    path.parent.mkdir(parents=True, exist_ok=True)
    
    # Backup existing file if requested
    if backup and MODULE_BACKUP_ENABLED and path.exists():
        backup_path = _backup_file(path)
        if backup_path:
            result["backup_path"] = str(backup_path)
    
    # Write to temp file first, then rename (atomic on most filesystems)
    temp_path = path.with_suffix(path.suffix + ".tmp")
    try:
        temp_path.write_text(content, encoding="utf-8")
        result["bytes_written"] = temp_path.stat().st_size
        
        # Atomic rename
        os.replace(temp_path, path)
        
    except Exception as e:
        # Clean up temp file on error
        if temp_path.exists():
            temp_path.unlink()
        raise e
    
    return result


# --- LangChain Tools ---------------------------------------------------------

@tool
def write_module_to_file(
    design_name: str,
    output_dir: str,
    include_meta: bool = True,
    overwrite: bool = False,
) -> str:
    """
    Write a designed Nextflow module to a file in the specified directory.
    
    Creates a proper nf-core style module structure:
        <output_dir>/modules/local/<name>/main.nf
        <output_dir>/modules/local/<name>/meta.yml (optional)
    
    Args:
        design_name: Name of the design to write (from save_module_design)
        output_dir: Directory to write to (must be in MODULE_WRITE_ROOTS)
        include_meta: Whether to generate meta.yml file (default: True)
        overwrite: Whether to overwrite existing files (default: False)
    
    Returns:
        Confirmation with file paths and any warnings.
    
    Note: By default, only writes designs with status "approved". Use
          write_draft_module() to write unapproved designs.
    """
    # Load the design
    design = _load_design(design_name)
    if not design:
        available = [d.get("name") for d in _list_designs()]
        return (
            "Error: No design named '" + design_name + "' found.\n"
            "Available designs: " + (", ".join(available) if available else "none") + "\n"
            "Use list_module_designs() to see all designs."
        )
    
    # Check status
    if design.get("status") != "approved":
        return (
            "Error: Design '" + design_name + "' has status '" + 
            design.get("status", "unknown") + "', not 'approved'.\n"
            "Either:\n"
            "  1. Get user approval: record_design_feedback('" + design_name + "', '...', approved=True)\n"
            "  2. Use write_draft_module('" + design_name + "', '" + output_dir + "') to write anyway"
        )
    
    # Validate output directory
    try:
        base_dir = _resolve_write_safe(output_dir)
        _ensure_write_root_exists(base_dir)
    except ValueError as e:
        return "Error: " + str(e)
    
    # Create module directory structure
    module_name = design.get("name", design_name)
    module_dir = base_dir / MODULE_SUBDIR / module_name.lower()
    main_nf_path = module_dir / "main.nf"
    meta_yml_path = module_dir / "meta.yml"
    
    # Check for existing files
    if main_nf_path.exists() and not overwrite:
        return (
            "Error: Module file already exists: " + str(main_nf_path) + "\n"
            "Use overwrite=True to replace it, or choose a different output_dir."
        )
    
    # Generate content
    main_nf_content = _generate_module_file(design)
    
    # Write files
    results = []
    try:
        write_result = _write_file_atomic(main_nf_path, main_nf_content)
        action = "Created" if write_result["created"] else "Updated"
        results.append("✓ " + action + ": " + str(main_nf_path))
        if write_result["backup_path"]:
            results.append("  Backup: " + write_result["backup_path"])
        
        if include_meta:
            meta_content = _generate_meta_yml(design)
            meta_result = _write_file_atomic(meta_yml_path, meta_content)
            action = "Created" if meta_result["created"] else "Updated"
            results.append("✓ " + action + ": " + str(meta_yml_path))
    
    except Exception as e:
        return "Error writing module files: " + str(e)
    
    # Format response
    lines = ["=== Module Written: " + module_name + " ===", ""]
    lines.extend(results)
    lines.append("")
    lines.append("To use in your pipeline:")
    lines.append("  include { " + module_name + " } from './" + MODULE_SUBDIR + "/" + module_name.lower() + "/main.nf'")
    
    return "\n".join(lines)


@tool
def write_draft_module(
    design_name: str,
    output_dir: str,
    include_meta: bool = False,
) -> str:
    """
    Write a DRAFT (unapproved) module to a file. Use this when you need
    to write a module that hasn't been formally approved yet.
    
    The output file will include a warning comment indicating it's a draft.
    
    Args:
        design_name: Name of the design to write
        output_dir: Directory to write to (must be in MODULE_WRITE_ROOTS)
        include_meta: Whether to generate meta.yml file (default: False for drafts)
    
    Returns:
        Confirmation with file paths and draft warning.
    """
    # Load the design
    design = _load_design(design_name)
    if not design:
        return "Error: No design named '" + design_name + "' found."
    
    # Validate output directory
    try:
        base_dir = _resolve_write_safe(output_dir)
        _ensure_write_root_exists(base_dir)
    except ValueError as e:
        return "Error: " + str(e)
    
    # Create module directory structure
    module_name = design.get("name", design_name)
    module_dir = base_dir / MODULE_SUBDIR / module_name.lower()
    main_nf_path = module_dir / "main.nf"
    
    # Generate content with draft warning
    draft_warning_lines = [
        "// ============================================",
        "// WARNING: This is a DRAFT module",
        "// Status: " + design.get("status", "draft"),
        "// This module has NOT been approved for use.",
        "// Review and test thoroughly before using.",
        "// ============================================",
        "",
    ]
    
    main_nf_content = "\n".join(draft_warning_lines) + _generate_module_file(design)
    
    # Write files
    try:
        write_result = _write_file_atomic(main_nf_path, main_nf_content)
        
        lines = ["=== DRAFT Module Written: " + module_name + " ===", ""]
        lines.append("⚠️  WARNING: This is a DRAFT (status: " + design.get("status", "draft") + ")")
        lines.append("")
        lines.append("✓ Written: " + str(main_nf_path))
        
        if include_meta:
            meta_content = _generate_meta_yml(design)
            meta_yml_path = module_dir / "meta.yml"
            _write_file_atomic(meta_yml_path, meta_content)
            lines.append("✓ Written: " + str(meta_yml_path))
        
        lines.append("")
        lines.append("To approve this design:")
        lines.append("  record_design_feedback('" + design_name + "', 'Approved for use', approved=True)")
        
        return "\n".join(lines)
    
    except Exception as e:
        return "Error writing draft module: " + str(e)


@tool
def write_all_approved_modules(
    output_dir: str,
    generate_includes: bool = True,
    overwrite: bool = False,
) -> str:
    """
    Write ALL approved module designs to the specified directory.
    
    Creates a complete module structure with all approved designs and
    optionally generates a modules.nf file that includes them all.
    
    Args:
        output_dir: Directory to write to (must be in MODULE_WRITE_ROOTS)
        generate_includes: Whether to create modules.nf with all includes (default: True)
        overwrite: Whether to overwrite existing files (default: False)
    
    Returns:
        Summary of all modules written.
    """
    # Get all approved designs
    approved = _list_designs(status_filter="approved")
    
    if not approved:
        draft_count = len(_list_designs(status_filter="draft"))
        return (
            "No approved designs found to write.\n"
            "Draft designs: " + str(draft_count) + "\n"
            "Use record_design_feedback(name, feedback, approved=True) to approve designs."
        )
    
    # Validate output directory
    try:
        base_dir = _resolve_write_safe(output_dir)
        _ensure_write_root_exists(base_dir)
    except ValueError as e:
        return "Error: " + str(e)
    
    # Write each module
    written = []
    errors = []
    
    for design in approved:
        module_name = design.get("name")
        if not module_name:
            continue
        
        module_dir = base_dir / MODULE_SUBDIR / module_name.lower()
        main_nf_path = module_dir / "main.nf"
        
        # Check for existing
        if main_nf_path.exists() and not overwrite:
            errors.append("Skipped " + module_name + ": file exists (use overwrite=True)")
            continue
        
        try:
            content = _generate_module_file(design)
            _write_file_atomic(main_nf_path, content)
            
            # Also write meta.yml
            meta_content = _generate_meta_yml(design)
            _write_file_atomic(module_dir / "meta.yml", meta_content)
            
            written.append(module_name)
        except Exception as e:
            errors.append("Failed " + module_name + ": " + str(e))
    
    # Generate includes file
    includes_path = None
    if generate_includes and written:
        includes_path = base_dir / "modules.nf"
        includes_content = _generate_includes_file(written, base_dir)
        try:
            _write_file_atomic(includes_path, includes_content)
        except Exception as e:
            errors.append("Failed to write modules.nf: " + str(e))
            includes_path = None
    
    # Format response
    lines = ["=== Batch Module Write: " + output_dir + " ===", ""]
    lines.append("Approved designs found: " + str(len(approved)))
    lines.append("Successfully written: " + str(len(written)))
    
    if written:
        lines.append("")
        lines.append("✓ Modules written:")
        for name in written:
            lines.append("    " + name + " → " + MODULE_SUBDIR + "/" + name.lower() + "/main.nf")
    
    if includes_path:
        lines.append("")
        lines.append("✓ Includes file: " + str(includes_path))
        lines.append("  Use: include { ... } from './modules.nf'")
    
    if errors:
        lines.append("")
        lines.append("⚠️ Issues (" + str(len(errors)) + "):")
        for err in errors:
            lines.append("    " + err)
    
    return "\n".join(lines)


@tool
def preview_module_file(design_name: str) -> str:
    """
    Preview what would be written for a module WITHOUT actually writing it.
    
    Use this to review the generated code before committing to disk.
    
    Args:
        design_name: Name of the design to preview
    
    Returns:
        The complete file content that would be written.
    """
    design = _load_design(design_name)
    if not design:
        return "Error: No design named '" + design_name + "' found."
    
    content = _generate_module_file(design)
    
    lines = [
        "=== Preview: " + design_name + " ===",
        "Status: " + design.get("status", "unknown"),
        "Revision: " + str(design.get("revision", "unknown")),
        "",
        "--- main.nf content ---",
        content,
        "--- end of preview ---",
        "",
        "To write this file:",
    ]
    
    if design.get("status") == "approved":
        lines.append("  write_module_to_file('" + design_name + "', '/output/path')")
    else:
        lines.append("  write_draft_module('" + design_name + "', '/output/path')")
        lines.append("  Or approve first: record_design_feedback('" + design_name + "', '...', approved=True)")
    
    return "\n".join(lines)


@tool
def create_pipeline_scaffold(
    output_dir: str,
    pipeline_name: str,
    include_modules: bool = True,
) -> str:
    """
    Create a basic Nextflow pipeline scaffold with standard structure.
    
    Creates:
        <output_dir>/
        ├── main.nf              # Main workflow entry point
        ├── nextflow.config      # Basic configuration
        ├── modules/
        │   └── local/           # Where local modules go
        └── README.md            # Basic documentation
    
    If include_modules is True and there are approved designs, they will
    be written to the modules/local directory.
    
    Args:
        output_dir: Directory to create pipeline in (must be in MODULE_WRITE_ROOTS)
        pipeline_name: Name for the pipeline
        include_modules: Whether to include approved module designs (default: True)
    
    Returns:
        Summary of created files and next steps.
    """
    # Validate output directory
    try:
        base_dir = _resolve_write_safe(output_dir)
        _ensure_write_root_exists(base_dir)
    except ValueError as e:
        return "Error: " + str(e)
    
    # Sanitize pipeline name
    safe_name = re.sub(r"[^A-Za-z0-9_-]", "_", pipeline_name)
    
    created_files = []
    
    # Create main.nf
    main_nf_content = _generate_main_nf(pipeline_name, safe_name)
    
    main_nf_path = base_dir / "main.nf"
    try:
        _write_file_atomic(main_nf_path, main_nf_content, backup=False)
        created_files.append("main.nf")
    except Exception as e:
        return "Error creating main.nf: " + str(e)
    
    # Create nextflow.config
    config_content = _generate_nextflow_config(pipeline_name, safe_name)
    
    config_path = base_dir / "nextflow.config"
    try:
        _write_file_atomic(config_path, config_content, backup=False)
        created_files.append("nextflow.config")
    except Exception as e:
        return "Error creating nextflow.config: " + str(e)
    
    # Create README
    readme_content = _generate_readme(pipeline_name)
    
    readme_path = base_dir / "README.md"
    try:
        _write_file_atomic(readme_path, readme_content, backup=False)
        created_files.append("README.md")
    except Exception:
        pass  # README is optional
    
    # Create modules directory
    modules_dir = base_dir / MODULE_SUBDIR
    modules_dir.mkdir(parents=True, exist_ok=True)
    created_files.append(MODULE_SUBDIR + "/")
    
    # Write approved modules if requested
    modules_written = []
    if include_modules:
        approved = _list_designs(status_filter="approved")
        for design in approved:
            module_name = design.get("name")
            if not module_name:
                continue
            
            module_dir = base_dir / MODULE_SUBDIR / module_name.lower()
            try:
                content = _generate_module_file(design)
                _write_file_atomic(module_dir / "main.nf", content, backup=False)
                modules_written.append(module_name)
            except Exception:
                pass
    
    # Format response
    lines = ["=== Pipeline Scaffold Created: " + pipeline_name + " ===", ""]
    lines.append("Location: " + str(base_dir))
    lines.append("")
    lines.append("Files created:")
    for f in created_files:
        lines.append("  ✓ " + f)
    
    if modules_written:
        lines.append("")
        lines.append("Approved modules included (" + str(len(modules_written)) + "):")
        for m in modules_written:
            lines.append("  ✓ " + m)
    
    lines.append("")
    lines.append("Next steps:")
    lines.append("  1. cd " + str(base_dir))
    lines.append("  2. Edit main.nf to add your workflow logic")
    lines.append("  3. Add modules with write_module_to_file(name, '" + output_dir + "')")
    lines.append("  4. Run: nextflow run main.nf --input <data>")
    
    return "\n".join(lines)


@tool
def list_written_modules(output_dir: str) -> str:
    """
    List all modules that have been written to a directory.
    
    Args:
        output_dir: Directory to check
    
    Returns:
        List of modules found with their paths.
    """
    try:
        base_dir = _resolve_write_safe(output_dir)
    except ValueError as e:
        return "Error: " + str(e)
    
    if not base_dir.exists():
        return "Directory does not exist: " + str(base_dir)
    
    modules_dir = base_dir / MODULE_SUBDIR
    if not modules_dir.exists():
        return "No modules directory found at: " + str(modules_dir)
    
    modules = []
    for item in sorted(modules_dir.iterdir()):
        if item.is_dir():
            main_nf = item / "main.nf"
            if main_nf.exists():
                size = main_nf.stat().st_size
                mtime = datetime.fromtimestamp(main_nf.stat().st_mtime)
                modules.append({
                    "name": item.name.upper(),
                    "path": str(main_nf.relative_to(base_dir)),
                    "size": size,
                    "modified": mtime.strftime("%Y-%m-%d %H:%M"),
                    "has_meta": (item / "meta.yml").exists(),
                })
    
    if not modules:
        return "No modules found in " + str(modules_dir)
    
    lines = ["=== Modules in " + str(base_dir) + " ===", ""]
    lines.append("Found " + str(len(modules)) + " module(s):")
    lines.append("")
    
    for m in modules:
        meta_mark = "📋" if m["has_meta"] else "  "
        lines.append("  " + meta_mark + " " + m["name"])
        lines.append("      Path: " + m["path"])
        lines.append("      Modified: " + m["modified"] + " (" + str(m["size"]) + " bytes)")
    
    lines.append("")
    lines.append("📋 = has meta.yml")
    
    return "\n".join(lines)


@tool
def show_write_permissions() -> str:
    """
    Show which directories the agent is allowed to write modules to.
    
    Returns:
        List of allowed write directories and their status.
    """
    lines = ["=== Module Write Permissions ===", ""]
    lines.append("The agent can write Nextflow modules to these directories:")
    lines.append("")
    
    for root in MODULE_WRITE_ROOTS:
        exists = root.exists()
        status = "✓ exists" if exists else "○ will be created"
        lines.append("  " + status + ": " + str(root))
    
    lines.append("")
    lines.append("To add more directories, set MODULE_WRITE_ROOTS in .env:")
    lines.append("  MODULE_WRITE_ROOTS=/path/one:/path/two")
    
    lines.append("")
    lines.append("Backups enabled: " + str(MODULE_BACKUP_ENABLED))
    
    return "\n".join(lines)
