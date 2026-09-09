"""
code_skill.py
=============

Lets the chat model analyze Python source files to understand their
purpose, structure, and how they work -- without executing them.

Tools exposed
-------------
    analyze_python_file(path)
        Extract docstrings, function/class signatures, imports, and
        a structural summary of a Python file.

    explain_python_project(path)
        Analyze all Python files in a directory to understand how
        they relate to each other (entry points, shared imports, etc).

    get_function_details(path, function_name)
        Get detailed info about a specific function including its
        full source code, docstring, and call signature.

Security: Uses fs_skill's ALLOWED_ROOTS check -- same restrictions apply.
"""

import ast
import os
from pathlib import Path
from typing import Optional

from langchain_core.tools import tool

from .fs_skill import _resolve_safe, MAX_READ_CHARS

MAX_SOURCE_CHARS = int(os.getenv("CODE_MAX_SOURCE_CHARS", "12000"))


def _safe_read_source(path: Path) -> str:
    """Read Python source with size limit."""
    text = path.read_text(encoding="utf-8", errors="replace")
    if len(text) > MAX_SOURCE_CHARS:
        return text[:MAX_SOURCE_CHARS] + f"\n# ... truncated at {MAX_SOURCE_CHARS} chars ..."
    return text


def _get_docstring(node) -> str:
    """Extract docstring from an AST node."""
    try:
        return ast.get_docstring(node) or ""
    except Exception:
        return ""


def _format_signature(node: ast.FunctionDef | ast.AsyncFunctionDef) -> str:
    """Build a readable function signature from an AST node."""
    args = []
    
    # Positional args
    defaults_offset = len(node.args.args) - len(node.args.defaults)
    for i, arg in enumerate(node.args.args):
        arg_str = arg.arg
        if arg.annotation:
            arg_str += f": {ast.unparse(arg.annotation)}"
        if i >= defaults_offset:
            default = node.args.defaults[i - defaults_offset]
            arg_str += f" = {ast.unparse(default)}"
        args.append(arg_str)
    
    # *args
    if node.args.vararg:
        arg_str = f"*{node.args.vararg.arg}"
        if node.args.vararg.annotation:
            arg_str += f": {ast.unparse(node.args.vararg.annotation)}"
        args.append(arg_str)
    
    # **kwargs
    if node.args.kwarg:
        arg_str = f"**{node.args.kwarg.arg}"
        if node.args.kwarg.annotation:
            arg_str += f": {ast.unparse(node.args.kwarg.annotation)}"
        args.append(arg_str)
    
    sig = f"def {node.name}({', '.join(args)})"
    if node.returns:
        sig += f" -> {ast.unparse(node.returns)}"
    return sig


def _analyze_ast(source: str, filename: str) -> dict:
    """
    Parse Python source and extract structural information.
    Returns a dict with imports, classes, functions, and module docstring.
    """
    try:
        tree = ast.parse(source, filename=filename)
    except SyntaxError as e:
        return {"error": f"Syntax error: {e}"}
    
    result = {
        "module_docstring": _get_docstring(tree),
        "imports": [],
        "classes": [],
        "functions": [],
        "global_variables": [],
        "has_main_block": False,
        "entry_points": [],
    }
    
    for node in ast.walk(tree):
        # Imports
        if isinstance(node, ast.Import):
            for alias in node.names:
                result["imports"].append(alias.name)
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            for alias in node.names:
                result["imports"].append(f"{module}.{alias.name}")
    
    # Top-level definitions only (not nested)
    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            methods = []
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    methods.append({
                        "name": item.name,
                        "signature": _format_signature(item),
                        "docstring": _get_docstring(item)[:200],
                    })
            result["classes"].append({
                "name": node.name,
                "docstring": _get_docstring(node),
                "methods": methods,
                "bases": [ast.unparse(b) for b in node.bases],
            })
        
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            func_info = {
                "name": node.name,
                "signature": _format_signature(node),
                "docstring": _get_docstring(node),
                "decorators": [ast.unparse(d) for d in node.decorator_list],
            }
            result["functions"].append(func_info)
            
            # Check for common entry point patterns
            if node.name == "main":
                result["entry_points"].append("main()")
            if "click.command" in str(func_info["decorators"]):
                result["entry_points"].append(f"CLI: {node.name}")
            if "app.route" in str(func_info["decorators"]):
                result["entry_points"].append(f"Flask route: {node.name}")
        
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    result["global_variables"].append(target.id)
        
        # Check for if __name__ == "__main__" block
        elif isinstance(node, ast.If):
            try:
                if (isinstance(node.test, ast.Compare) and
                    ast.unparse(node.test) == "__name__ == '__main__'"):
                    result["has_main_block"] = True
                    result["entry_points"].append("if __name__ == '__main__' block")
            except Exception:
                pass
    
    return result


@tool
def analyze_python_file(path: str) -> str:
    """
    Analyze a Python file to understand its purpose and structure.
    Returns the module docstring, imports, classes, functions (with
    signatures and docstrings), and identifies likely entry points.
    Use this to understand what a Python script does before deciding
    whether/how to run it. `path` must be inside an allow-listed
    root directory.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    if target.suffix.lower() != ".py":
        return f"Error: '{target}' is not a Python file."
    
    try:
        source = _safe_read_source(target)
    except Exception as e:
        return f"Error reading '{target}': {e}"
    
    analysis = _analyze_ast(source, target.name)
    
    if "error" in analysis:
        return f"Could not parse '{target}': {analysis['error']}"
    
    # Format output
    lines = [f"=== Analysis of {target.name} ===\n"]
    
    if analysis["module_docstring"]:
        lines.append(f"PURPOSE:\n{analysis['module_docstring']}\n")
    
    if analysis["entry_points"]:
        lines.append(f"ENTRY POINTS: {', '.join(analysis['entry_points'])}")
    
    if analysis["imports"]:
        # Group by top-level package
        std_lib = []
        third_party = []
        local = []
        for imp in analysis["imports"]:
            top = imp.split(".")[0]
            if imp.startswith(".") or top in ("skills", "utils", "lib"):
                local.append(imp)
            else:
                third_party.append(imp)
        
        lines.append(f"\nIMPORTS ({len(analysis['imports'])} total):")
        if third_party:
            lines.append(f"  External: {', '.join(sorted(set(third_party))[:15])}")
        if local:
            lines.append(f"  Local: {', '.join(sorted(set(local)))}")
    
    if analysis["classes"]:
        lines.append(f"\nCLASSES ({len(analysis['classes'])}):")
        for cls in analysis["classes"]:
            bases = f"({', '.join(cls['bases'])})" if cls["bases"] else ""
            lines.append(f"  class {cls['name']}{bases}")
            if cls["docstring"]:
                lines.append(f"    \"\"\"{cls['docstring'][:150]}...\"\"\"")
            for method in cls["methods"][:5]:
                lines.append(f"      {method['signature']}")
            if len(cls["methods"]) > 5:
                lines.append(f"      ... and {len(cls['methods']) - 5} more methods")
    
    if analysis["functions"]:
        lines.append(f"\nFUNCTIONS ({len(analysis['functions'])}):")
        for func in analysis["functions"]:
            dec = f"  @{func['decorators'][0]}\n" if func["decorators"] else ""
            lines.append(f"{dec}  {func['signature']}")
            if func["docstring"]:
                short_doc = func["docstring"].split("\n")[0][:100]
                lines.append(f"    \"\"\"{short_doc}\"\"\"")
    
    return "\n".join(lines)


@tool
def explain_python_project(path: str, max_files: int = 20) -> str:
    """
    Analyze all Python files in a directory to understand the project
    structure: what each file does, how they relate, and where execution
    likely starts. Good for understanding an unfamiliar codebase before
    deciding what to run. `path` must be inside an allow-listed root.
    """
    try:
        base = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not base.is_dir():
        return f"Error: '{base}' is not a directory."
    
    py_files = sorted(base.rglob("*.py"))[:max_files]
    if not py_files:
        return f"No Python files found in '{base}'."
    
    project_info = {
        "files": [],
        "all_imports": set(),
        "entry_points": [],
    }
    
    for py_file in py_files:
        try:
            source = _safe_read_source(py_file)
            analysis = _analyze_ast(source, py_file.name)
            
            rel_path = py_file.relative_to(base)
            file_info = {
                "path": str(rel_path),
                "purpose": (analysis.get("module_docstring") or "No docstring")[:200],
                "functions": len(analysis.get("functions", [])),
                "classes": len(analysis.get("classes", [])),
            }
            project_info["files"].append(file_info)
            project_info["all_imports"].update(analysis.get("imports", []))
            
            if analysis.get("entry_points"):
                for ep in analysis["entry_points"]:
                    project_info["entry_points"].append(f"{rel_path}: {ep}")
        
        except Exception as e:
            project_info["files"].append({
                "path": str(py_file.relative_to(base)),
                "error": str(e),
            })
    
    # Format output
    lines = [f"=== Python Project: {base.name} ===\n"]
    lines.append(f"Found {len(py_files)} Python file(s)")
    
    if project_info["entry_points"]:
        lines.append(f"\nLIKELY ENTRY POINTS:")
        for ep in project_info["entry_points"]:
            lines.append(f"  • {ep}")
    
    # Key dependencies
    external = sorted({
        imp.split(".")[0] for imp in project_info["all_imports"]
        if not imp.startswith(".")
    })[:20]
    if external:
        lines.append(f"\nKEY DEPENDENCIES: {', '.join(external)}")
    
    lines.append(f"\nFILES:")
    for f in project_info["files"]:
        if "error" in f:
            lines.append(f"  {f['path']}: [error: {f['error']}]")
        else:
            lines.append(f"  {f['path']}")
            lines.append(f"    {f['purpose'][:100]}")
            lines.append(f"    ({f['functions']} functions, {f['classes']} classes)")
    
    if len(py_files) >= max_files:
        lines.append(f"\n... truncated at {max_files} files ...")
    
    return "\n".join(lines)


@tool
def get_function_source(path: str, function_name: str) -> str:
    """
    Get the complete source code of a specific function from a Python
    file. Use this after analyze_python_file when you need to see the
    full implementation of a particular function to understand exactly
    what it does or what arguments it needs.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    
    try:
        source = target.read_text(encoding="utf-8", errors="replace")
        tree = ast.parse(source, filename=target.name)
    except Exception as e:
        return f"Error parsing '{target}': {e}"
    
    # Find the function
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name == function_name:
                # Extract source lines
                lines = source.splitlines()
                start = node.lineno - 1
                end = node.end_lineno if hasattr(node, "end_lineno") else start + 50
                func_source = "\n".join(lines[start:end])
                
                if len(func_source) > MAX_SOURCE_CHARS:
                    func_source = func_source[:MAX_SOURCE_CHARS] + "\n# ... truncated ..."
                
                return f"# From {target.name}, line {node.lineno}\n\n{func_source}"
    
    return f"Function '{function_name}' not found in '{target}'."
