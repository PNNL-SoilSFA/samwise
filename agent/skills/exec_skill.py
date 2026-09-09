"""
exec_skill.py
=============

Lets the chat model execute Python scripts in a controlled manner and
collect their output. This is the "do something" counterpart to
code_skill.py's "understand something" tools.

Tools exposed
-------------
    run_python_script(path, args=[], timeout=300)
        Execute a Python script and capture stdout/stderr.

    run_python_function(path, function_name, kwargs={}, timeout=120)
        Import a specific function from a file and call it with given
        arguments. Returns the function's return value (if JSON-serializable)
        or its string representation.

    run_python_snippet(code, timeout=60)
        Execute a short snippet of Python code directly. Good for quick
        data exploration, transformations, or computations.

Security
--------
    - Path restrictions: same ALLOWED_ROOTS as fs_skill
    - Execution environment: runs in a subprocess with configurable timeout
    - Optional: can restrict to a whitelist of "approved" scripts
    - Resource limits: timeout, output size limits

Configure via .env:
    EXEC_TIMEOUT      - default timeout in seconds (default: 300)
    EXEC_MAX_OUTPUT   - max chars of stdout/stderr to return (default: 50000)
    EXEC_PYTHON       - Python interpreter path (default: sys.executable)
    EXEC_ALLOW_SNIPPET - set to "false" to disable run_python_snippet
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from langchain_core.tools import tool

from .fs_skill import _resolve_safe, ALLOWED_ROOTS

EXEC_TIMEOUT = int(os.getenv("EXEC_TIMEOUT", "300"))
EXEC_MAX_OUTPUT = int(os.getenv("EXEC_MAX_OUTPUT", "50000"))
EXEC_PYTHON = os.getenv("EXEC_PYTHON", sys.executable)
EXEC_ALLOW_SNIPPET = os.getenv("EXEC_ALLOW_SNIPPET", "true").lower() != "false"

# Optional: whitelist specific scripts that are known-safe to run.
# If set, only these scripts can be executed. Empty = allow any in ALLOWED_ROOTS.
_APPROVED_SCRIPTS_RAW = os.getenv("EXEC_APPROVED_SCRIPTS", "")
APPROVED_SCRIPTS = [
    Path(p).resolve() for p in _APPROVED_SCRIPTS_RAW.split(":") if p.strip()
] if _APPROVED_SCRIPTS_RAW else []


def _check_script_approved(path: Path) -> bool:
    """If APPROVED_SCRIPTS is set, check that path is in the whitelist."""
    if not APPROVED_SCRIPTS:
        return True  # No whitelist = allow anything in ALLOWED_ROOTS
    return path.resolve() in APPROVED_SCRIPTS


def _truncate_output(text: str, max_chars: int = EXEC_MAX_OUTPUT) -> str:
    if len(text) > max_chars:
        return text[:max_chars] + f"\n\n... truncated at {max_chars} characters ..."
    return text


def _format_run_result(
    returncode: int,
    stdout: str,
    stderr: str,
    timeout: bool = False
) -> str:
    """Format subprocess results for the model."""
    lines = []
    
    if timeout:
        lines.append("⚠️  TIMEOUT: Process exceeded time limit and was killed.\n")
    
    lines.append(f"Exit code: {returncode}")
    
    if stdout.strip():
        lines.append(f"\n--- STDOUT ---\n{_truncate_output(stdout)}")
    else:
        lines.append("\n--- STDOUT ---\n(empty)")
    
    if stderr.strip():
        lines.append(f"\n--- STDERR ---\n{_truncate_output(stderr)}")
    
    return "\n".join(lines)


@tool
def run_python_script(
    path: str,
    args: list[str] | None = None,
    timeout: int = EXEC_TIMEOUT,
    working_dir: str | None = None,
) -> str:
    """
    Execute a Python script and return its stdout/stderr output.
    
    Use this when you need to actually run a script (after understanding
    it with analyze_python_file) to collect results, generate data, or
    perform computations. The script runs in a subprocess with the
    specified timeout.
    
    Args:
        path: Path to the Python script (must be in ALLOWED_ROOTS)
        args: Command-line arguments to pass to the script
        timeout: Max seconds to wait (default from EXEC_TIMEOUT env var)
        working_dir: Directory to run from (defaults to script's directory)
    
    Returns:
        Combined stdout/stderr with exit code. Check exit code 0 for success.
    """
    args = args or []
    
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    if target.suffix.lower() != ".py":
        return f"Error: '{target}' is not a Python file."
    if not _check_script_approved(target):
        return (
            f"Error: '{target}' is not in the approved scripts list. "
            f"Add it to EXEC_APPROVED_SCRIPTS in .env if execution is intended."
        )
    
    # Determine working directory
    if working_dir:
        try:
            cwd = _resolve_safe(working_dir)
        except ValueError as e:
            return f"Error with working_dir: {e}"
    else:
        cwd = target.parent
    
    cmd = [EXEC_PYTHON, str(target)] + [str(a) for a in args]
    
    try:
        result = subprocess.run(
            cmd,
            cwd=str(cwd),
            capture_output=True,
            text=True,
            timeout=timeout,
            env={**os.environ, "PYTHONUNBUFFERED": "1"},
        )
        return _format_run_result(result.returncode, result.stdout, result.stderr)
    
    except subprocess.TimeoutExpired as e:
        stdout = e.stdout or ""
        stderr = e.stderr or ""
        if isinstance(stdout, bytes):
            stdout = stdout.decode(errors="replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode(errors="replace")
        return _format_run_result(-1, stdout, stderr, timeout=True)
    
    except Exception as e:
        return f"Error executing script: {e}"


@tool
def run_python_function(
    path: str,
    function_name: str,
    kwargs: dict | None = None,
    timeout: int = 120,
) -> str:
    """
    Import and call a specific function from a Python file with the
    given keyword arguments. Returns the function's return value.
    
    This is useful when a script has utility functions you want to call
    directly without running the whole script, or when you need to pass
    structured arguments that would be awkward on the command line.
    
    Args:
        path: Path to the Python file containing the function
        function_name: Name of the function to call
        kwargs: Keyword arguments to pass to the function (as a dict)
        timeout: Max seconds to wait for the function to return
    
    Returns:
        The function's return value (JSON-encoded if possible) or error.
    """
    kwargs = kwargs or {}
    
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    if not _check_script_approved(target):
        return f"Error: '{target}' is not in the approved scripts list."
    
    # Build a wrapper script that imports and calls the function
    wrapper_code = f'''
import sys
import json
sys.path.insert(0, {repr(str(target.parent))})

from {target.stem} import {function_name}

kwargs = json.loads({repr(json.dumps(kwargs))})
result = {function_name}(**kwargs)

# Try to JSON-serialize the result; fall back to repr
try:
    print("__RESULT_JSON__")
    print(json.dumps(result, indent=2, default=str))
except Exception:
    print("__RESULT_REPR__")
    print(repr(result))
'''
    
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".py", delete=False
    ) as tmp:
        tmp.write(wrapper_code)
        tmp_path = tmp.name
    
    try:
        result = subprocess.run(
            [EXEC_PYTHON, tmp_path],
            cwd=str(target.parent),
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        
        if result.returncode != 0:
            return f"Function call failed:\n{result.stderr}"
        
        stdout = result.stdout
        if "__RESULT_JSON__" in stdout:
            return stdout.split("__RESULT_JSON__")[1].strip()
        elif "__RESULT_REPR__" in stdout:
            return stdout.split("__RESULT_REPR__")[1].strip()
        else:
            return _truncate_output(stdout)
    
    except subprocess.TimeoutExpired:
        return f"Function call timed out after {timeout} seconds."
    except Exception as e:
        return f"Error calling function: {e}"
    finally:
        Path(tmp_path).unlink(missing_ok=True)


@tool
def run_python_snippet(code: str, timeout: int = 60) -> str:
    """
    Execute a short Python code snippet and return its output.
    
    Good for quick computations, data transformations, or exploring
    data that you've read with other tools. The code runs in an
    isolated subprocess.
    
    Args:
        code: Python code to execute (can be multi-line)
        timeout: Max seconds to wait (default: 60)
    
    Returns:
        stdout/stderr from the code execution.
    
    Example:
        run_python_snippet('''
        import pandas as pd
        df = pd.read_csv("/path/to/data.csv")
        print(df.describe())
        ''')
    """
    if not EXEC_ALLOW_SNIPPET:
        return (
            "Error: run_python_snippet is disabled on this system. "
            "Use run_python_script to run existing scripts instead."
        )
    
    # Basic sanity check on code length
    if len(code) > 10000:
        return "Error: Code snippet too long (max 10000 characters)."
    
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".py", delete=False
    ) as tmp:
        tmp.write(code)
        tmp_path = tmp.name
    
    try:
        result = subprocess.run(
            [EXEC_PYTHON, tmp_path],
            capture_output=True,
            text=True,
            timeout=timeout,
            env={**os.environ, "PYTHONUNBUFFERED": "1"},
        )
        return _format_run_result(result.returncode, result.stdout, result.stderr)
    
    except subprocess.TimeoutExpired as e:
        stdout = (e.stdout or b"").decode(errors="replace")
        stderr = (e.stderr or b"").decode(errors="replace")
        return _format_run_result(-1, stdout, stderr, timeout=True)
    except Exception as e:
        return f"Error executing snippet: {e}"
    finally:
        Path(tmp_path).unlink(missing_ok=True)


@tool
def plan_script_execution(path: str) -> str:
    """
    Analyze a Python script or project directory and suggest how to run
    it: what the entry points are, what arguments might be needed, and
    what order to run things in if there are multiple scripts.
    
    Use this before run_python_script when you're not sure how to
    execute something.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    
    # Import the analysis function from code_skill
    from .code_skill import _analyze_ast, _safe_read_source
    
    lines = [f"=== Execution Plan for {target} ===\n"]
    
    if target.is_file():
        if target.suffix.lower() != ".py":
            return f"'{target}' is not a Python file."
        
        source = _safe_read_source(target)
        analysis = _analyze_ast(source, target.name)
        
        if "error" in analysis:
            return f"Could not parse: {analysis['error']}"
        
        # Determine how to run it
        if analysis.get("has_main_block"):
            lines.append("✓ Has `if __name__ == '__main__'` block")
            lines.append(f"\nRun with:\n  run_python_script('{path}')")
        
        # Check for CLI frameworks
        imports = set(analysis.get("imports", []))
        if any("click" in i for i in imports):
            lines.append("\n✓ Uses Click CLI framework")
            lines.append("  Run with --help to see available commands:")
            lines.append(f"  run_python_script('{path}', args=['--help'])")
        elif any("argparse" in i for i in imports):
            lines.append("\n✓ Uses argparse")
            lines.append("  Run with --help to see available arguments:")
            lines.append(f"  run_python_script('{path}', args=['--help'])")
        elif any("fire" in i for i in imports):
            lines.append("\n✓ Uses Python Fire CLI")
            lines.append("  Run with --help or call functions directly")
        
        # Callable functions
        funcs = analysis.get("functions", [])
        public_funcs = [f for f in funcs if not f["name"].startswith("_")]
        if public_funcs:
            lines.append(f"\nCallable functions ({len(public_funcs)}):")
            for f in public_funcs[:10]:
                lines.append(f"  • {f['signature']}")
            lines.append(f"\nCall directly with:")
            lines.append(f"  run_python_function('{path}', 'function_name', {{'arg': 'value'}})")
    
    elif target.is_dir():
        # Find all Python files and identify entry points
        py_files = list(target.rglob("*.py"))
        
        entry_points = []
        for py_file in py_files[:30]:
            try:
                source = _safe_read_source(py_file)
                analysis = _analyze_ast(source, py_file.name)
                if analysis.get("has_main_block") or analysis.get("entry_points"):
                    entry_points.append((py_file, analysis))
            except Exception:
                continue
        
        if entry_points:
            lines.append("Entry points found:\n")
            for py_file, analysis in entry_points:
                rel = py_file.relative_to(target)
                lines.append(f"  {rel}")
                for ep in analysis.get("entry_points", []):
                    lines.append(f"    → {ep}")
        
        # Check for common patterns
        setup_py = target / "setup.py"
        pyproject = target / "pyproject.toml"
        main_py = target / "main.py"
        run_py = target / "run.py"
        
        if main_py.exists():
            lines.append(f"\n✓ Has main.py - likely entry point")
            lines.append(f"  run_python_script('{main_py}')")
        if run_py.exists():
            lines.append(f"\n✓ Has run.py - likely entry point")
            lines.append(f"  run_python_script('{run_py}')")
        if setup_py.exists() or pyproject.exists():
            lines.append("\n✓ Has package config - may need `pip install -e .` first")
    
    return "\n".join(lines)
