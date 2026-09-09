"""
nextflow_skill.py
=================

Specialized skill for finding, parsing, and explaining Nextflow pipelines.
Nextflow is a workflow language common in bioinformatics for defining
data-driven computational pipelines.

Tools exposed
-------------
    find_nextflow_pipelines(path)
        Discover all Nextflow files (.nf), configs, and related files
        in a directory tree.

    analyze_nextflow_pipeline(path)
        Parse a .nf file and extract processes, channels, workflows,
        inputs/outputs, and explain what the pipeline does.

    explain_nextflow_config(path)
        Parse nextflow.config and explain execution profiles, parameters,
        resource settings, and container definitions.

    summarize_nextflow_project(path)
        Comprehensive overview of a Nextflow project: main workflow,
        all processes, parameters, and how to run it.

    get_nextflow_process(path, process_name)
        Get detailed information about a specific process including
        its full script, inputs, outputs, and directives.

Security: Uses fs_skill's ALLOWED_ROOTS check.
"""

import os
import re
from pathlib import Path
from typing import Optional

from langchain_core.tools import tool

from .fs_skill import _resolve_safe, MAX_READ_CHARS

MAX_SOURCE_CHARS = int(os.getenv("NF_MAX_SOURCE_CHARS", "15000"))


# --- Nextflow Parsing Helpers -------------------------------------------------

def _extract_processes(content: str) -> list[dict]:
    """
    Extract process definitions from Nextflow DSL2 code.
    Returns list of {name, inputs, outputs, script, directives, docstring}.
    """
    processes = []
    
    # Match process blocks: process NAME { ... }
    # This regex handles nested braces by finding process starts
    process_pattern = re.compile(
        r'process\s+(\w+)\s*\{',
        re.MULTILINE
    )
    
    for match in process_pattern.finditer(content):
        name = match.group(1)
        start = match.end()
        
        # Find matching closing brace (handle nesting)
        brace_count = 1
        pos = start
        while pos < len(content) and brace_count > 0:
            if content[pos] == '{':
                brace_count += 1
            elif content[pos] == '}':
                brace_count -= 1
            pos += 1
        
        process_body = content[start:pos-1]
        
        process_info = {
            "name": name,
            "inputs": _extract_section(process_body, "input"),
            "outputs": _extract_section(process_body, "output"),
            "script": _extract_script(process_body),
            "directives": _extract_directives(process_body),
            "when": _extract_section(process_body, "when"),
            "docstring": _extract_process_comment(content, match.start()),
        }
        processes.append(process_info)
    
    return processes


def _extract_section(body: str, section_name: str) -> str:
    """Extract a section (input:, output:, etc.) from a process body."""
    pattern = re.compile(
        rf'{section_name}\s*:\s*\n(.*?)(?=\n\s*(?:input|output|script|shell|exec|stub|when|directive)\s*:|$)',
        re.DOTALL | re.IGNORECASE
    )
    match = pattern.search(body)
    if match:
        return match.group(1).strip()
    return ""


def _extract_script(body: str) -> str:
    """Extract the script/shell/exec block from a process."""
    # Try script, shell, or exec blocks
    for keyword in ["script", "shell", "exec"]:
        pattern = re.compile(
            rf'{keyword}\s*:\s*\n?\s*(""".*?"""|\'\'\'.*?\'\'\'|".*?"|\'.*?\')',
            re.DOTALL
        )
        match = pattern.search(body)
        if match:
            script = match.group(1)
            # Remove surrounding quotes
            if script.startswith('"""') or script.startswith("'''"):
                script = script[3:-3]
            elif script.startswith('"') or script.startswith("'"):
                script = script[1:-1]
            return script.strip()
    
    # Fallback: look for triple-quoted block
    triple_pattern = re.compile(r'(""".*?"""|\'\'\'.*?\'\'\')', re.DOTALL)
    match = triple_pattern.search(body)
    if match:
        script = match.group(1)[3:-3]
        return script.strip()
    
    return "(script not parsed)"


def _extract_directives(body: str) -> dict:
    """Extract process directives (cpus, memory, container, etc.)."""
    directives = {}
    
    directive_patterns = [
        (r'cpus\s+(\d+|params\.\w+)', 'cpus'),
        (r'memory\s+[\'"]?([^\'"\n]+)[\'"]?', 'memory'),
        (r'time\s+[\'"]?([^\'"\n]+)[\'"]?', 'time'),
        (r'container\s+[\'"]([^"\']+)[\'"]', 'container'),
        (r'conda\s+[\'"]([^"\']+)[\'"]', 'conda'),
        (r'publishDir\s+[\'"]?([^\'",\n]+)', 'publishDir'),
        (r'tag\s+[\'"]?([^\'"\n]+)', 'tag'),
        (r'label\s+[\'"](\w+)[\'"]', 'label'),
        (r'errorStrategy\s+[\'"]?(\w+)', 'errorStrategy'),
        (r'maxRetries\s+(\d+)', 'maxRetries'),
    ]
    
    for pattern, key in directive_patterns:
        match = re.search(pattern, body)
        if match:
            directives[key] = match.group(1).strip()
    
    return directives


def _extract_process_comment(content: str, process_start: int) -> str:
    """Extract comment/docstring immediately before a process definition."""
    # Look backwards for /* */ or // comments
    before = content[:process_start].rstrip()
    
    # Block comment
    if before.endswith("*/"):
        block_start = before.rfind("/*")
        if block_start != -1:
            comment = before[block_start+2:-2].strip()
            # Clean up asterisks from javadoc-style comments
            lines = [re.sub(r'^\s*\*\s?', '', line) for line in comment.split('\n')]
            return '\n'.join(lines).strip()
    
    # Line comments
    lines = before.split('\n')
    comment_lines = []
    for line in reversed(lines):
        stripped = line.strip()
        if stripped.startswith('//'):
            comment_lines.insert(0, stripped[2:].strip())
        elif stripped == '':
            continue
        else:
            break
    
    return '\n'.join(comment_lines)


def _extract_workflows(content: str) -> list[dict]:
    """Extract workflow definitions from Nextflow DSL2 code."""
    workflows = []
    
    # Named workflows: workflow NAME { ... }
    # Entry workflow: workflow { ... }
    workflow_pattern = re.compile(
        r'workflow\s*(\w*)\s*\{',
        re.MULTILINE
    )
    
    for match in workflow_pattern.finditer(content):
        name = match.group(1) or "(entry)"
        start = match.end()
        
        # Find matching closing brace
        brace_count = 1
        pos = start
        while pos < len(content) and brace_count > 0:
            if content[pos] == '{':
                brace_count += 1
            elif content[pos] == '}':
                brace_count -= 1
            pos += 1
        
        body = content[start:pos-1]
        
        # Extract process calls
        process_calls = re.findall(r'(\w+)\s*\(', body)
        # Filter out common functions
        process_calls = [p for p in process_calls if p not in 
                        ('Channel', 'file', 'tuple', 'path', 'val', 'emit', 
                         'take', 'main', 'if', 'else', 'println')]
        
        # Extract channel operations
        channel_ops = re.findall(r'\.(map|filter|flatten|collect|groupTuple|join|combine|mix|concat|first|last|take|splitCsv|splitFasta|splitFastq)\s*[\(\{]', body)
        
        workflows.append({
            "name": name,
            "process_calls": list(dict.fromkeys(process_calls)),  # unique, preserve order
            "channel_operations": list(set(channel_ops)),
            "body_preview": body[:500] + "..." if len(body) > 500 else body,
        })
    
    return workflows


def _extract_params(content: str) -> dict:
    """Extract params definitions from Nextflow code or config."""
    params = {}
    
    # params.name = value
    pattern = re.compile(r'params\.(\w+)\s*=\s*([^\n]+)')
    for match in pattern.finditer(content):
        name = match.group(1)
        value = match.group(2).strip().rstrip(';')
        params[name] = value
    
    # params { name = value } block
    block_pattern = re.compile(r'params\s*\{([^}]+)\}', re.DOTALL)
    for block_match in block_pattern.finditer(content):
        block = block_match.group(1)
        for line in block.split('\n'):
            line = line.strip()
            if '=' in line and not line.startswith('//'):
                parts = line.split('=', 1)
                name = parts[0].strip()
                value = parts[1].strip().rstrip(';')
                params[name] = value
    
    return params


def _extract_includes(content: str) -> list[dict]:
    """Extract include statements from Nextflow code."""
    includes = []
    
    # include { PROCESS } from './module'
    # include { PROCESS as ALIAS } from './module'
    pattern = re.compile(
        r'include\s*\{\s*([^}]+)\s*\}\s*from\s*[\'"]([^\'"]+)[\'"]'
    )
    
    for match in pattern.finditer(content):
        items_str = match.group(1)
        source = match.group(2)
        
        # Parse individual items (may have aliases)
        items = []
        for item in items_str.split(';'):
            item = item.strip()
            if ' as ' in item:
                original, alias = item.split(' as ')
                items.append({"name": original.strip(), "alias": alias.strip()})
            elif item:
                items.append({"name": item.strip()})
        
        includes.append({"source": source, "items": items})
    
    return includes


def _parse_nextflow_config(content: str) -> dict:
    """Parse nextflow.config and extract key settings."""
    config = {
        "params": _extract_params(content),
        "profiles": {},
        "process_defaults": {},
        "docker": None,
        "singularity": None,
        "conda": None,
        "executor": None,
        "manifest": {},
    }
    
    # Extract profiles
    profile_pattern = re.compile(r'profiles\s*\{([^}]+(?:\{[^}]*\}[^}]*)*)\}', re.DOTALL)
    profile_match = profile_pattern.search(content)
    if profile_match:
        profiles_block = profile_match.group(1)
        # Find individual profile names
        profile_names = re.findall(r'(\w+)\s*\{', profiles_block)
        for name in profile_names:
            config["profiles"][name] = f"(profile '{name}' defined)"
    
    # Docker/Singularity enabled
    if re.search(r'docker\s*\{\s*enabled\s*=\s*true', content):
        config["docker"] = "enabled"
    if re.search(r'singularity\s*\{\s*enabled\s*=\s*true', content):
        config["singularity"] = "enabled"
    
    # Executor
    executor_match = re.search(r'executor\s*=\s*[\'"](\w+)[\'"]', content)
    if executor_match:
        config["executor"] = executor_match.group(1)
    
    # Manifest
    manifest_pattern = re.compile(r'manifest\s*\{([^}]+)\}', re.DOTALL)
    manifest_match = manifest_pattern.search(content)
    if manifest_match:
        manifest_block = manifest_match.group(1)
        for line in manifest_block.split('\n'):
            if '=' in line:
                key, val = line.split('=', 1)
                config["manifest"][key.strip()] = val.strip().strip("'\"")
    
    return config


# --- Tools --------------------------------------------------------------------

import os
import re
from pathlib import Path
from typing import Optional, Generator
import time

from .fs_skill import _resolve_safe

# Configuration
MAX_SEARCH_TIME = int(os.getenv("NF_MAX_SEARCH_TIME", "30"))  # seconds
MAX_SEARCH_DEPTH = int(os.getenv("NF_MAX_SEARCH_DEPTH", "5"))  # directory depth
MAX_FILES_TO_CHECK = int(os.getenv("NF_MAX_FILES", "1000"))


def _limited_rglob(base: Path, pattern: str, max_depth: int = MAX_SEARCH_DEPTH, 
                   max_files: int = MAX_FILES_TO_CHECK, 
                   timeout: float = MAX_SEARCH_TIME) -> Generator[Path, None, None]:
    """
    Like Path.rglob but with depth limit, file count limit, and timeout.
    Much faster for large directory trees.
    """
    start_time = time.time()
    files_yielded = 0
    
    def _search(current: Path, depth: int):
        nonlocal files_yielded
        
        if depth > max_depth:
            return
        if time.time() - start_time > timeout:
            return
        if files_yielded >= max_files:
            return
            
        try:
            for item in current.iterdir():
                if time.time() - start_time > timeout:
                    return
                if files_yielded >= max_files:
                    return
                    
                try:
                    if item.is_file():
                        if item.match(pattern):
                            files_yielded += 1
                            yield item
                    elif item.is_dir():
                        # Skip common non-relevant directories
                        skip_dirs = {
                            '.git', '.nextflow', 'work', '.snakemake', 
                            '__pycache__', 'node_modules', '.cache',
                            'results', 'output', 'outputs', 'logs',
                            '.singularity', '.conda'
                        }
                        if item.name not in skip_dirs and not item.name.startswith('.nextflow'):
                            yield from _search(item, depth + 1)
                except PermissionError:
                    continue
                except OSError:
                    continue
        except PermissionError:
            return
        except OSError:
            return
    
    yield from _search(base, 0)


def _quick_find_nextflow(base: Path, timeout: float = MAX_SEARCH_TIME) -> dict:
    """
    Quickly find Nextflow-related files with smart prioritization.
    Checks common locations first before doing a full search.
    """
    start_time = time.time()
    
    results = {
        "nf_files": [],
        "config_files": [],
        "module_dirs": [],
        "other_related": [],
        "search_complete": True,
        "search_time": 0,
    }
    
    # PHASE 1: Check common/expected locations first (instant)
    common_main_files = ["main.nf", "workflow.nf", "pipeline.nf", "nextflow.nf"]
    for name in common_main_files:
        candidate = base / name
        if candidate.exists():
            results["nf_files"].append(str(candidate.relative_to(base)))
    
    # Check for config in root
    for config_name in ["nextflow.config", "nextflow.config.template"]:
        candidate = base / config_name
        if candidate.exists():
            results["config_files"].append(str(candidate.relative_to(base)))
    
    # Check common module directories
    for mod_dir in ["modules", "subworkflows", "lib", "bin", "workflows"]:
        candidate = base / mod_dir
        if candidate.is_dir():
            results["module_dirs"].append(f"{mod_dir}/")
    
    # If we found main files, do a quick targeted search for more
    if results["nf_files"] or results["config_files"]:
        # PHASE 2: Quick search in expected locations only
        search_dirs = [base]
        for mod_dir in ["modules", "subworkflows", "workflows", "lib"]:
            if (base / mod_dir).is_dir():
                search_dirs.append(base / mod_dir)
        
        for search_dir in search_dirs:
            if time.time() - start_time > timeout:
                results["search_complete"] = False
                break
            
            # Only go 3 levels deep in module directories
            for nf_file in _limited_rglob(search_dir, "*.nf", max_depth=3, timeout=timeout/2):
                try:
                    rel = str(nf_file.relative_to(base))
                    if rel not in results["nf_files"]:
                        results["nf_files"].append(rel)
                except ValueError:
                    continue
                
                if len(results["nf_files"]) >= 50:
                    break
    else:
        # PHASE 2b: No obvious Nextflow files found - do broader search but with limits
        for nf_file in _limited_rglob(base, "*.nf", max_depth=MAX_SEARCH_DEPTH, 
                                       max_files=100, timeout=timeout):
            if time.time() - start_time > timeout:
                results["search_complete"] = False
                break
            try:
                rel = str(nf_file.relative_to(base))
                results["nf_files"].append(rel)
            except ValueError:
                continue
        
        # Also look for configs
        for config_file in _limited_rglob(base, "nextflow.config", max_depth=3, 
                                          max_files=20, timeout=timeout/4):
            try:
                rel = str(config_file.relative_to(base))
                if rel not in results["config_files"]:
                    results["config_files"].append(rel)
            except ValueError:
                continue
    
    results["search_time"] = round(time.time() - start_time, 2)
    return results

@tool
def explain_pipeline_directory(path: str, max_files_to_analyze: int = 5) -> str:
    """
    🎯 ONE-SHOT pipeline explanation. Discovers and analyzes Nextflow
    pipelines in a directory in a SINGLE tool call.

    Use this when user asks to "explain" a pipeline directory. This avoids
    multiple round-trips by doing discovery + analysis together.

    Args:
        path: Directory containing the pipeline (must be in ALLOWED_ROOTS)
        max_files_to_analyze: Max .nf files to parse in detail (default: 5)

    Returns comprehensive summary including:
    - All discovered .nf files and configs
    - Parsed processes and workflows from main files
    - Parameters and how to run
    """
    try:
        base = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"

    if not base.is_dir():
        return f"Error: '{base}' is not a directory."

    lines = [f"{'='*60}\n📂 PIPELINE EXPLANATION: {base.name}\n{'='*60}\n"]

    # Step 1: Quick discovery
    discovery = _quick_find_nextflow(base, timeout=15)

    lines.append(f"DISCOVERY (searched in {discovery['search_time']}s):")
    lines.append(f"  • .nf workflow files: {len(discovery['nf_files'])}")
    lines.append(f"  • Config files: {len(discovery['config_files'])}")
    lines.append(f"  • Module directories: {len(discovery['module_dirs'])}")

    if not discovery["nf_files"]:
        # Check subdirectories for pipelines
        lines.append(f"\nNo .nf files in root. Checking subdirectories...")
        subdirs_with_nf = []
        try:
            for subdir in sorted(base.iterdir()):
                if subdir.is_dir() and not subdir.name.startswith('.'):
                    sub_nf = list(subdir.glob("*.nf"))[:5]
                    if sub_nf:
                        subdirs_with_nf.append((subdir.name, [f.name for f in sub_nf]))
        except PermissionError:
            pass

        if subdirs_with_nf:
            lines.append(f"\nSubdirectories containing Nextflow pipelines:")
            for dirname, nf_files in subdirs_with_nf[:10]:
                lines.append(f"\n  📁 {dirname}/")
                for nf in nf_files[:3]:
                    lines.append(f"      • {nf}")
                if len(nf_files) > 3:
                    lines.append(f"      • ... and {len(nf_files) - 3} more")

            # Analyze the first subdirectory with pipelines
            if subdirs_with_nf:
                first_subdir = base / subdirs_with_nf[0][0]
                lines.append(f"\n{'─'*50}")
                lines.append(f"Analyzing: {subdirs_with_nf[0][0]}/")
                lines.append(f"{'─'*50}")

                # Re-run discovery on subdirectory
                discovery = _quick_find_nextflow(first_subdir, timeout=10)
        else:
            lines.append("\n❌ No Nextflow pipelines found in this directory or its subdirectories.")
            lines.append("\nThis might be a data/results directory rather than a pipeline source.")
            return "\n".join(lines)

    # Step 2: List all discovered files
    if discovery["nf_files"]:
        lines.append(f"\n{'─'*50}")
        lines.append("WORKFLOW FILES:")

        # Prioritize main workflow files
        main_files = []
        module_files = []
        for f in discovery["nf_files"]:
            fname = Path(f).name.lower()
            if fname in ("main.nf", "workflow.nf", "pipeline.nf"):
                main_files.insert(0, f)
            elif "module" in fname or "module" in f.lower():
                module_files.append(f)
            else:
                main_files.append(f)

        sorted_files = main_files + module_files

        for f in sorted_files[:15]:
            marker = "⭐" if Path(f).name.lower() in ("main.nf", "workflow.nf") else "📄"
            lines.append(f"  {marker} {f}")
        if len(sorted_files) > 15:
            lines.append(f"  ... and {len(sorted_files) - 15} more files")

    if discovery["config_files"]:
        lines.append(f"\nCONFIG FILES:")
        for f in discovery["config_files"]:
            lines.append(f"  ⚙️  {f}")

    # Step 3: Analyze top workflow files
    lines.append(f"\n{'─'*50}")
    lines.append(f"DETAILED ANALYSIS (top {max_files_to_analyze} files):")
    lines.append(f"{'─'*50}")

    all_processes = []
    all_params = {}
    all_workflows = []
    files_analyzed = 0

    # Determine base path for file resolution
    if discovery["nf_files"]:
        analysis_base = base
        # Check if files are relative to a subdirectory
        first_file = discovery["nf_files"][0]
        if not (base / first_file).exists():
            # Files might be from subdirectory discovery
            for subdir in base.iterdir():
                if subdir.is_dir() and (subdir / Path(first_file).name).exists():
                    analysis_base = subdir
                    break

        # Sort to analyze main files first
        files_to_analyze = sorted(discovery["nf_files"],
            key=lambda f: (0 if Path(f).name.lower() in ("main.nf", "workflow.nf") else 1, f))

        for nf_file in files_to_analyze[:max_files_to_analyze]:
            file_path = analysis_base / nf_file
            if not file_path.exists():
                file_path = base / nf_file
            if not file_path.exists():
                continue

            try:
                content = file_path.read_text(encoding="utf-8", errors="replace")
                if len(content) > 30000:
                    content = content[:30000]

                processes = _extract_processes(content)
                workflows = _extract_workflows(content)
                params = _extract_params(content)

                lines.append(f"\n📄 {nf_file}")

                # File docstring
                if content.strip().startswith("/*"):
                    end = content.find("*/")
                    if end != -1 and end < 500:
                        doc = content[2:end].strip()
                        doc = re.sub(r'^\s*\*\s?', '', doc, flags=re.MULTILINE)
                        if doc:
                            lines.append(f"   \"{doc[:150]}\"")

                if processes:
                    lines.append(f"   Processes: {', '.join(p['name'] for p in processes[:6])}")
                    if len(processes) > 6:
                        lines.append(f"              ... and {len(processes) - 6} more")
                    all_processes.extend(processes)

                if workflows:
                    for wf in workflows:
                        if wf["process_calls"]:
                            flow = " → ".join(wf["process_calls"][:8])
                            lines.append(f"   Flow: {flow}")
                            if len(wf["process_calls"]) > 8:
                                lines.append(f"         ... ({len(wf['process_calls']) - 8} more steps)")
                    all_workflows.extend(workflows)

                all_params.update(params)
                files_analyzed += 1

            except Exception as e:
                lines.append(f"\n📄 {nf_file}: (parse error: {e})")

    # Step 4: Summary
    lines.append(f"\n{'─'*50}")
    lines.append("SUMMARY:")
    lines.append(f"{'─'*50}")
    lines.append(f"  Total .nf files: {len(discovery['nf_files'])}")
    lines.append(f"  Files analyzed: {files_analyzed}")
    lines.append(f"  Total processes found: {len(all_processes)}")
    lines.append(f"  Total parameters: {len(all_params)}")

    if all_processes:
        lines.append(f"\nALL PROCESSES ({len(all_processes)}):")
        for proc in all_processes[:20]:
            lines.append(f"  • {proc['name']}")
        if len(all_processes) > 20:
            lines.append(f"  ... and {len(all_processes) - 20} more")

    if all_params:
        lines.append(f"\nPARAMETERS ({len(all_params)}):")
        for name, val in list(all_params.items())[:10]:
            lines.append(f"  --{name} = {val}")
        if len(all_params) > 10:
            lines.append(f"  ... and {len(all_params) - 10} more")

    # Step 5: How to run
    lines.append(f"\n{'─'*50}")
    lines.append("HOW TO RUN:")
    lines.append(f"{'─'*50}")

    if discovery["nf_files"]:
        main_candidate = None
        for f in discovery["nf_files"]:
            if Path(f).name.lower() in ("main.nf", "workflow.nf", "pipeline.nf"):
                main_candidate = f
                break
        if not main_candidate:
            main_candidate = discovery["nf_files"][0]

        lines.append(f"  cd {base}")
        lines.append(f"  nextflow run {main_candidate}")
        if all_params:
            example_param = list(all_params.keys())[0]
            lines.append(f"  nextflow run {main_candidate} --{example_param} <value>")

    return "\n".join(lines)

@tool
def find_nextflow_pipelines(path: str, thorough: bool = False) -> str:
    """
    🔍 FAST discovery of Nextflow pipelines in a directory.
    
    Quickly finds .nf workflow files, configs, and module directories
    by checking common locations first, then doing a bounded search.
    
    Use this INSTEAD of list_directory when looking for Nextflow/pipeline code.
    
    Args:
        path: Directory to search (must be in ALLOWED_ROOTS)
        thorough: If True, do a more complete (slower) search. Default False
                  uses smart shortcuts for speed.
    
    Works for: DRAM, DRAM2, nf-core pipelines, custom Nextflow projects
    """
    try:
        base = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not base.is_dir():
        return f"Error: '{base}' is not a directory."
    
    timeout = MAX_SEARCH_TIME * 2 if thorough else MAX_SEARCH_TIME
    results = _quick_find_nextflow(base, timeout=timeout)
    
    # Format output
    lines = [f"=== Nextflow Pipeline Discovery: {base.name} ===\n"]
    lines.append(f"(searched in {results['search_time']}s)\n")
    
    if results["nf_files"]:
        lines.append(f"WORKFLOW FILES ({len(results['nf_files'])}):")
        
        # Sort to show main files first
        def sort_key(f):
            name = Path(f).name.lower()
            if name == "main.nf":
                return (0, f)
            if name in ("workflow.nf", "pipeline.nf"):
                return (1, f)
            if "/" not in f:  # root level
                return (2, f)
            return (3, f)
        
        sorted_files = sorted(results["nf_files"], key=sort_key)
        for f in sorted_files[:30]:
            marker = "⭐" if Path(f).name in ("main.nf", "workflow.nf") else "📄"
            lines.append(f"  {marker} {f}")
        if len(results["nf_files"]) > 30:
            lines.append(f"  ... and {len(results['nf_files']) - 30} more ...")
    else:
        lines.append("❌ No .nf workflow files found.")
    
    if results["config_files"]:
        lines.append(f"\nCONFIG FILES:")
        for f in results["config_files"]:
            lines.append(f"  ⚙️  {f}")
    
    if results["module_dirs"]:
        lines.append(f"\nMODULE DIRECTORIES:")
        for d in results["module_dirs"]:
            lines.append(f"  📁 {d}")
    
    if not results["search_complete"]:
        lines.append(f"\n⚠️  Search timed out - results may be incomplete.")
        lines.append(f"   Use find_nextflow_pipelines(path, thorough=True) for deeper search.")
    
    # Suggest next steps
    if results["nf_files"]:
        main_file = None
        for f in results["nf_files"]:
            if Path(f).name in ("main.nf", "workflow.nf", "pipeline.nf"):
                main_file = f
                break
        if not main_file:
            main_file = results["nf_files"][0]
        
        lines.append(f"\n💡 Next: analyze_nextflow_pipeline('{path}/{main_file}')")
        lines.append(f"   Or: summarize_nextflow_project('{path}')")
    
    return "\n".join(lines)


@tool
def summarize_nextflow_project(path: str) -> str:
    """
    📋 FAST comprehensive overview of a Nextflow project.
    
    Quickly identifies the main workflow, processes, parameters, and
    how to run the pipeline. Optimized for large directories.
    
    Use when user asks to "explain" a pipeline like DRAM2, or wants
    to understand what a Nextflow project does.
    """
    try:
        base = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not base.is_dir():
        return f"Error: '{base}' is not a directory."
    
    lines = [f"=== Nextflow Project: {base.name} ===\n"]
    
    # Quick find main workflow (check common names first)
    main_nf = None
    for candidate_name in ["main.nf", "workflow.nf", "pipeline.nf", f"{base.name}.nf"]:
        candidate = base / candidate_name
        if candidate.exists():
            main_nf = candidate
            break
    
    # If not found, quick search root only
    if not main_nf:
        root_nf_files = list(base.glob("*.nf"))
        if root_nf_files:
            main_nf = root_nf_files[0]
    
    # Quick find config
    config_file = None
    if (base / "nextflow.config").exists():
        config_file = base / "nextflow.config"
    
    # Quick count modules (don't traverse deeply)
    module_count = 0
    for mod_dir in ["modules", "subworkflows"]:
        mod_path = base / mod_dir
        if mod_path.is_dir():
            # Just count immediate .nf files, don't recurse deeply
            module_count += len(list(mod_path.glob("*.nf")))
            module_count += len(list(mod_path.glob("*/*.nf")))  # one level deep
    
    # Overview
    lines.append("PROJECT STRUCTURE:")
    lines.append(f"  Main workflow: {main_nf.name if main_nf else 'Not found (searched root)'}")
    lines.append(f"  Config: {config_file.name if config_file else 'Not found'}")
    if module_count:
        lines.append(f"  Module files: ~{module_count}")
    
    # Analyze main workflow if found
    if main_nf:
        try:
            content = main_nf.read_text(encoding="utf-8", errors="replace")
            # Limit content size for parsing
            if len(content) > 50000:
                content = content[:50000]
            
            processes = _extract_processes(content)
            workflows = _extract_workflows(content)
            params = _extract_params(content)
            includes = _extract_includes(content)
            
            lines.append(f"\n--- Main Workflow: {main_nf.name} ---")
            
            # File docstring
            if content.strip().startswith("/*"):
                end = content.find("*/")
                if end != -1 and end < 1000:
                    doc = content[2:end].strip()
                    doc = re.sub(r'^\s*\*\s?', '', doc, flags=re.MULTILINE)
                    lines.append(f"\n{doc[:300]}")
            
            if includes:
                lines.append(f"\nImports from {len(includes)} module(s):")
                for inc in includes[:8]:
                    lines.append(f"  • {inc['source']}")
                if len(includes) > 8:
                    lines.append(f"  • ... and {len(includes) - 8} more")
            
            if processes:
                lines.append(f"\nProcesses ({len(processes)}):")
                for proc in processes[:10]:
                    desc = proc["docstring"][:50] + "..." if proc["docstring"] else ""
                    lines.append(f"  • {proc['name']} {desc}")
                if len(processes) > 10:
                    lines.append(f"  • ... and {len(processes) - 10} more")
            
            if workflows:
                for wf in workflows:
                    if wf["name"] == "(entry)" or wf["name"] == "":
                        if wf["process_calls"]:
                            lines.append(f"\nExecution flow:")
                            flow = " → ".join(wf["process_calls"][:12])
                            lines.append(f"  {flow}")
                            if len(wf["process_calls"]) > 12:
                                lines.append(f"  → ... ({len(wf['process_calls']) - 12} more steps)")
            
            if params:
                lines.append(f"\nParameters ({len(params)}):")
                for name, val in list(params.items())[:8]:
                    lines.append(f"  --{name} = {val}")
                if len(params) > 8:
                    lines.append(f"  ... and {len(params) - 8} more")
                    
        except Exception as e:
            lines.append(f"\n(Could not fully parse workflow: {e})")
    
    # Config summary if found
    if config_file:
        try:
            config_content = config_file.read_text(encoding="utf-8", errors="replace")
            if len(config_content) > 20000:
                config_content = config_content[:20000]
            
            config = _parse_nextflow_config(config_content)
            
            if config["profiles"]:
                lines.append(f"\nProfiles: {', '.join(list(config['profiles'].keys())[:6])}")
            
            if config["params"] and not params:  # Only if not already shown
                lines.append(f"\nConfig parameters ({len(config['params'])}):")
                for name, val in list(config["params"].items())[:5]:
                    lines.append(f"  --{name} = {val}")
        except Exception:
            pass
    
    # How to run
    lines.append(f"\n--- How to Run ---")
    if main_nf:
        lines.append(f"  nextflow run {main_nf.name}")
        lines.append(f"  nextflow run {main_nf.name} --help")
        if config_file:
            lines.append(f"  nextflow run {main_nf.name} -profile <profile>")
    else:
        lines.append("  Main workflow not found - use find_nextflow_pipelines to locate .nf files")
    
    return "\n".join(lines)

@tool
def analyze_nextflow_pipeline(path: str) -> str:
    """
    Parse and analyze a Nextflow workflow file (.nf). Extracts:
    - Pipeline description/comments
    - All processes with their inputs, outputs, scripts, and resources
    - Workflow definitions showing process execution order
    - Channel definitions and data flow
    - Parameters (params.*)
    - Module includes
    
    Use this to understand what a Nextflow pipeline does, what inputs
    it expects, and how processes connect together.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    if target.suffix.lower() != ".nf":
        return f"Error: '{target}' is not a .nf file."
    
    try:
        content = target.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return f"Error reading '{target}': {e}"
    
    if len(content) > MAX_SOURCE_CHARS:
        content = content[:MAX_SOURCE_CHARS]
    
    # Extract components
    processes = _extract_processes(content)
    workflows = _extract_workflows(content)
    params = _extract_params(content)
    includes = _extract_includes(content)
    
    # Extract file-level docstring (top comment)
    file_doc = ""
    if content.strip().startswith("/*"):
        end = content.find("*/")
        if end != -1:
            file_doc = content[2:end].strip()
            # Clean javadoc style
            file_doc = re.sub(r'^\s*\*\s?', '', file_doc, flags=re.MULTILINE)
    elif content.strip().startswith("//"):
        lines = []
        for line in content.split('\n'):
            if line.strip().startswith("//"):
                lines.append(line.strip()[2:].strip())
            elif line.strip() == "":
                continue
            else:
                break
        file_doc = '\n'.join(lines)
    
    # Detect DSL version
    dsl_version = "DSL2" if "DSL2" in content or "workflow {" in content else "DSL1 (legacy)"
    
    # Format output
    lines = [f"=== Nextflow Pipeline: {target.name} ===\n"]
    lines.append(f"DSL Version: {dsl_version}")
    
    if file_doc:
        lines.append(f"\nDESCRIPTION:\n{file_doc[:500]}")
    
    if includes:
        lines.append(f"\nINCLUDES ({len(includes)}):")
        for inc in includes:
            items = ", ".join(i["name"] for i in inc["items"])
            lines.append(f"  from '{inc['source']}': {items}")
    
    if params:
        lines.append(f"\nPARAMETERS ({len(params)}):")
        for name, value in list(params.items())[:15]:
            lines.append(f"  params.{name} = {value}")
        if len(params) > 15:
            lines.append(f"  ... and {len(params) - 15} more ...")
    
    if processes:
        lines.append(f"\nPROCESSES ({len(processes)}):")
        for proc in processes:
            lines.append(f"\n  📦 {proc['name']}")
            if proc["docstring"]:
                lines.append(f"     \"{proc['docstring'][:100]}\"")
            if proc["inputs"]:
                inputs_short = proc["inputs"][:100].replace('\n', ' ')
                lines.append(f"     inputs: {inputs_short}")
            if proc["outputs"]:
                outputs_short = proc["outputs"][:100].replace('\n', ' ')
                lines.append(f"     outputs: {outputs_short}")
            if proc["directives"]:
                dirs = ", ".join(f"{k}={v}" for k, v in list(proc["directives"].items())[:4])
                lines.append(f"     directives: {dirs}")
    
    if workflows:
        lines.append(f"\nWORKFLOWS ({len(workflows)}):")
        for wf in workflows:
            name_display = wf["name"] if wf["name"] != "(entry)" else "(main entry point)"
            lines.append(f"\n  🔀 workflow {name_display}")
            if wf["process_calls"]:
                lines.append(f"     calls: {' → '.join(wf['process_calls'])}")
            if wf["channel_operations"]:
                lines.append(f"     channel ops: {', '.join(wf['channel_operations'])}")
    
    # Suggest how to run
    lines.append(f"\n--- How to Run ---")
    lines.append(f"  nextflow run {target.name}")
    if params:
        example_param = list(params.keys())[0]
        lines.append(f"  nextflow run {target.name} --{example_param} <value>")
    
    return "\n".join(lines)


@tool
def explain_nextflow_config(path: str) -> str:
    """
    Parse and explain a nextflow.config file. Extracts:
    - Default parameters
    - Execution profiles (docker, singularity, conda, slurm, etc.)
    - Process resource defaults
    - Container settings
    - Executor configuration
    
    Use this to understand how to run a pipeline and what options are available.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    
    try:
        content = target.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return f"Error reading '{target}': {e}"
    
    config = _parse_nextflow_config(content)
    
    lines = [f"=== Nextflow Config: {target.name} ===\n"]
    
    if config["manifest"]:
        lines.append("MANIFEST:")
        for key, val in config["manifest"].items():
            lines.append(f"  {key}: {val}")
        lines.append("")
    
    if config["params"]:
        lines.append(f"PARAMETERS ({len(config['params'])}):")
        for name, value in list(config["params"].items())[:20]:
            lines.append(f"  --{name}  (default: {value})")
        if len(config["params"]) > 20:
            lines.append(f"  ... and {len(config['params']) - 20} more ...")
        lines.append("")
    
    if config["profiles"]:
        lines.append(f"PROFILES ({len(config['profiles'])}):")
        for name in config["profiles"]:
            lines.append(f"  -profile {name}")
        lines.append("")
    
    env_info = []
    if config["docker"]:
        env_info.append("Docker: enabled")
    if config["singularity"]:
        env_info.append("Singularity: enabled")
    if config["conda"]:
        env_info.append("Conda: enabled")
    if config["executor"]:
        env_info.append(f"Executor: {config['executor']}")
    
    if env_info:
        lines.append("EXECUTION ENVIRONMENT:")
        for info in env_info:
            lines.append(f"  {info}")
        lines.append("")
    
    # Usage hints
    lines.append("--- Usage Examples ---")
    lines.append(f"  nextflow run main.nf -profile <profile>")
    if config["params"]:
        lines.append(f"  nextflow run main.nf --param_name value")
    
    return "\n".join(lines)


@tool
def summarize_nextflow_project(path: str) -> str:
    """
    Comprehensive overview of a Nextflow project directory. Combines
    find_nextflow_pipelines, analyze_nextflow_pipeline, and
    explain_nextflow_config to give a complete picture of:
    
    - What the pipeline does
    - All processes and their purpose
    - Required inputs and parameters
    - How to run it
    
    Use this for a complete understanding of an unfamiliar Nextflow project.
    """
    try:
        base = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not base.is_dir():
        return f"Error: '{base}' is not a directory."
    
    lines = [f"=== Nextflow Project Summary: {base.name} ===\n"]
    
    # Find main workflow file
    main_nf = None
    for candidate in ["main.nf", "workflow.nf", "pipeline.nf"]:
        if (base / candidate).exists():
            main_nf = base / candidate
            break
    
    if not main_nf:
        # Find any .nf file
        nf_files = list(base.glob("*.nf"))
        if nf_files:
            main_nf = nf_files[0]
    
    # Find config
    config_file = None
    if (base / "nextflow.config").exists():
        config_file = base / "nextflow.config"
    
    # Count modules
    modules_dir = base / "modules"
    module_count = len(list(modules_dir.rglob("*.nf"))) if modules_dir.exists() else 0
    
    subworkflows_dir = base / "subworkflows"
    subworkflow_count = len(list(subworkflows_dir.rglob("*.nf"))) if subworkflows_dir.exists() else 0
    
    # Overview stats
    lines.append("PROJECT STRUCTURE:")
    lines.append(f"  Main workflow: {main_nf.name if main_nf else 'Not found'}")
    lines.append(f"  Config file: {config_file.name if config_file else 'Not found'}")
    lines.append(f"  Modules: {module_count}")
    lines.append(f"  Subworkflows: {subworkflow_count}")
    
    # Parse main workflow
    if main_nf:
        try:
            content = main_nf.read_text(encoding="utf-8", errors="replace")
            processes = _extract_processes(content)
            workflows = _extract_workflows(content)
            params = _extract_params(content)
            includes = _extract_includes(content)
            
            lines.append(f"\n--- Main Workflow: {main_nf.name} ---")
            
            if includes:
                lines.append(f"\nImports {len(includes)} module(s):")
                for inc in includes[:10]:
                    lines.append(f"  • {inc['source']}")
            
            if processes:
                lines.append(f"\nDefines {len(processes)} process(es):")
                for proc in processes[:10]:
                    desc = proc["docstring"][:60] if proc["docstring"] else "(no description)"
                    lines.append(f"  • {proc['name']}: {desc}")
            
            if workflows:
                for wf in workflows:
                    if wf["name"] == "(entry)" or not wf["name"]:
                        lines.append(f"\nExecution flow:")
                        lines.append(f"  {' → '.join(wf['process_calls'][:15])}")
                        if len(wf["process_calls"]) > 15:
                            lines.append(f"  ... and more")
        except Exception as e:
            lines.append(f"\n(Could not parse main workflow: {e})")
    
    # Parse config
    if config_file:
        try:
            content = config_file.read_text(encoding="utf-8", errors="replace")
            config = _parse_nextflow_config(content)
            
            lines.append(f"\n--- Configuration ---")
            
            if config["profiles"]:
                lines.append(f"Available profiles: {', '.join(config['profiles'].keys())}")
            
            if config["params"]:
                lines.append(f"\nKey parameters ({len(config['params'])} total):")
                for name, val in list(config["params"].items())[:10]:
                    lines.append(f"  --{name} = {val}")
        except Exception as e:
            lines.append(f"\n(Could not parse config: {e})")
    
    # How to run
    lines.append(f"\n--- How to Run ---")
    if main_nf:
        lines.append(f"  cd {base}")
        lines.append(f"  nextflow run {main_nf.name}")
        if config_file:
            lines.append(f"  nextflow run {main_nf.name} -profile <profile>")
        lines.append(f"\nView parameters:")
        lines.append(f"  nextflow run {main_nf.name} --help")
    else:
        lines.append("  No main workflow file found - check directory contents.")
    
    return "\n".join(lines)


@tool
def get_nextflow_process(path: str, process_name: str) -> str:
    """
    Get detailed information about a specific Nextflow process, including
    its full script, all inputs/outputs, directives, and any documentation.
    
    Use this after analyze_nextflow_pipeline when you need to understand
    exactly what a particular process does.
    
    Args:
        path: Path to the .nf file containing the process
        process_name: Name of the process to examine
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    
    try:
        content = target.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return f"Error reading: {e}"
    
    processes = _extract_processes(content)
    
    for proc in processes:
        if proc["name"].lower() == process_name.lower():
            lines = [f"=== Process: {proc['name']} ===\n"]
            
            if proc["docstring"]:
                lines.append(f"DESCRIPTION:\n{proc['docstring']}\n")
            
            if proc["inputs"]:
                lines.append(f"INPUTS:\n{proc['inputs']}\n")
            
            if proc["outputs"]:
                lines.append(f"OUTPUTS:\n{proc['outputs']}\n")
            
            if proc["when"]:
                lines.append(f"WHEN (conditional):\n{proc['when']}\n")
            
            if proc["directives"]:
                lines.append("DIRECTIVES:")
                for key, val in proc["directives"].items():
                    lines.append(f"  {key}: {val}")
                lines.append("")
            
            if proc["script"]:
                lines.append(f"SCRIPT:\n```\n{proc['script']}\n```")
            
            return "\n".join(lines)
    
    available = [p["name"] for p in processes]
    return f"Process '{process_name}' not found in '{target}'.\nAvailable processes: {', '.join(available)}"
