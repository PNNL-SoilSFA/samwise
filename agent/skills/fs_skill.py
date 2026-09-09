"""
fs_skill.py
===========

Filesystem-access skill for LangChain: lets the chat model browse and
read files from one or more allow-listed directories on disk, and
(optionally) build a one-off vector index over an arbitrary directory
for semantic search -- without needing to pre-configure DOCS_DIR.

This complements rag_skill.py, which only ever indexes the single
directory pointed to by DOCS_DIR at startup. fs_skill.py is for cases
like:

    "Can you see and read the files in /rcfs/projects/.../some_run?"
    "Can you index the documents in /some/other/path?"

Tools exposed
-------------
    list_directory(path, recursive=True)
        List files/subfolders under `path` with size + dir/file marker.

    read_file(path, max_chars=8000)
        Return the raw text content of a single file (txt, csv, tsv,
        json, log, md, etc). Truncated to max_chars.

    index_directory(path, force_rebuild=False)
        Build (or reuse) a Chroma vector index over every .txt/.md/.pdf
        /.docx file under `path`. Good for prose-style documents; NOT
        useful for structured output like CSV/TSV/JSON -- use
        list_directory + read_file for those instead.

    search_indexed_directory(path, query, k=4)
        Semantic search over a directory previously indexed with
        index_directory.

Security
--------
The model can only reach paths under one or more *allow-listed* root
directories -- never the whole filesystem. Configure this via .env:

    ALLOWED_ROOTS=/rcfs/projects/samwise/microtrait/microtrait_out:./docs

(colon-separated on Linux/macOS; semicolon on Windows -- handled
automatically below.) If ALLOWED_ROOTS is unset, this skill falls back
to DOCS_DIR only, so it's at least as restricted as rag_skill.py by
default -- it never silently opens up broader access.
"""

import hashlib
import json
import os
from pathlib import Path

from langchain_chroma import Chroma
from langchain_core.tools import tool
from langchain_text_splitters import RecursiveCharacterTextSplitter

from .rag_skill import LOADERS_BY_EXT, VECTOR_DB_DIR, build_embeddings

# --- Configuration -----------------------------------------------------------

_ROOTS_SEP = ";" if os.name == "nt" else ":"

# Used only if ALLOWED_ROOTS is not set in .env at all. Add/remove paths
# here for a permanent default, or (preferred) set ALLOWED_ROOTS in .env
# so this list doesn't need a code change + restart every time.
_DEFAULT_ROOTS = [
    os.getenv("DOCS_DIR", "./docs"),
    "/rcfs/projects/samwise/",
]

_raw_roots = os.getenv("ALLOWED_ROOTS", "")
ALLOWED_ROOTS = [
    Path(p).expanduser().resolve() for p in _raw_roots.split(_ROOTS_SEP) if p.strip()
] or [Path(p).expanduser().resolve() for p in _DEFAULT_ROOTS]

MAX_LIST_ENTRIES = int(os.getenv("FS_MAX_LIST_ENTRIES", "200"))
MAX_READ_CHARS = int(os.getenv("FS_MAX_READ_CHARS", "8000"))

# Per-directory indexes built by index_directory() live here, one
# sub-folder per indexed path (named by a short hash of the resolved path).
DYNAMIC_INDEX_ROOT = Path(VECTOR_DB_DIR) / "_dynamic"


# --- Path safety --------------------------------------------------------------

def _resolve_safe(path_str: str) -> Path:
    """
    Resolve a model-supplied path and ensure it falls inside one of
    ALLOWED_ROOTS. Raises ValueError otherwise. Tool functions catch this
    and turn it into a plain error string -- the model never sees a stack
    trace, and a bad path never reaches the actual filesystem call.
    """
    candidate = Path(path_str).expanduser().resolve()
    for root in ALLOWED_ROOTS:
        try:
            candidate.relative_to(root)
            return candidate
        except ValueError:
            continue
    allowed = ", ".join(str(r) for r in ALLOWED_ROOTS)
    raise ValueError(
        f"'{candidate}' is outside the allowed director"
        f"{'y' if len(ALLOWED_ROOTS) == 1 else 'ies'} ({allowed}). "
        f"Add it to ALLOWED_ROOTS in .env if this access is intended."
    )


# --- Browsing tools ------------------------------------------------------------

import fnmatch
from typing import Optional

@tool
def find_files(
    path: str,
    pattern: str = "*.py",
    max_results: int = 100,
    include_size: bool = True
) -> str:
    """
    Search for files matching a glob pattern under `path`. Much more
    efficient than list_directory when you're looking for specific file
    types (like *.py, *.sh, *.R, etc.) in a large directory tree.
    
    Args:
        path: Directory to search (must be in ALLOWED_ROOTS)
        pattern: Glob pattern to match filenames (default: "*.py")
                 Examples: "*.py", "*.R", "run_*.sh", "config*.json"
        max_results: Maximum number of files to return (default: 100)
        include_size: Whether to show file sizes (default: True)
    
    Returns:
        List of matching files with their paths and sizes.
    """
    try:
        base = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not base.exists():
        return f"Error: '{base}' does not exist."
    if not base.is_dir():
        return f"Error: '{base}' is not a directory."
    
    matches = []
    try:
        for p in base.rglob(pattern):
            if p.is_file():
                try:
                    rel = p.relative_to(base)
                    if include_size:
                        size = p.stat().st_size
                        matches.append(f"{rel}  ({size:,} bytes)")
                    else:
                        matches.append(str(rel))
                except (ValueError, OSError):
                    continue
            
            if len(matches) >= max_results:
                break
    except PermissionError as e:
        return f"Error: Permission denied while searching: {e}"
    
    if not matches:
        return f"No files matching '{pattern}' found in '{base}'."
    
    result = [f"Found {len(matches)} file(s) matching '{pattern}' in {base}:\n"]
    result.extend(f"  {m}" for m in sorted(matches))
    
    if len(matches) >= max_results:
        result.append(f"\n... limited to first {max_results} results ...")
    
    return "\n".join(result)


@tool
def find_files_by_content(
    path: str,
    search_text: str,
    file_pattern: str = "*.py",
    max_files: int = 50,
    context_lines: int = 2
) -> str:
    """
    Search for files containing specific text. Useful for finding scripts
    that use particular functions, import certain modules, or contain
    specific configuration.
    
    Args:
        path: Directory to search (must be in ALLOWED_ROOTS)
        search_text: Text to search for within files
        file_pattern: Only search files matching this pattern (default: "*.py")
        max_files: Max files to search through (default: 50)
        context_lines: Lines of context around matches (default: 2)
    
    Returns:
        Files containing the search text with matching line excerpts.
    """
    try:
        base = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not base.is_dir():
        return f"Error: '{base}' is not a directory."
    
    matching_files = []
    files_searched = 0
    
    for p in base.rglob(file_pattern):
        if not p.is_file():
            continue
        if files_searched >= max_files:
            break
        files_searched += 1
        
        try:
            content = p.read_text(encoding="utf-8", errors="replace")
            if search_text.lower() in content.lower():
                # Find matching lines with context
                lines = content.splitlines()
                excerpts = []
                for i, line in enumerate(lines):
                    if search_text.lower() in line.lower():
                        start = max(0, i - context_lines)
                        end = min(len(lines), i + context_lines + 1)
                        excerpt = "\n".join(
                            f"  {j+1}: {lines[j]}" 
                            for j in range(start, end)
                        )
                        excerpts.append(f"  Line {i+1}:\n{excerpt}")
                        if len(excerpts) >= 3:  # Max 3 excerpts per file
                            break
                
                matching_files.append({
                    "path": str(p.relative_to(base)),
                    "excerpts": excerpts[:3]
                })
        except Exception:
            continue
    
    if not matching_files:
        return f"No files matching '{file_pattern}' containing '{search_text}' found in '{base}'."
    
    result = [f"Found '{search_text}' in {len(matching_files)} file(s):\n"]
    for mf in matching_files[:20]:
        result.append(f"\n📄 {mf['path']}")
        for exc in mf["excerpts"]:
            result.append(exc)
    
    if len(matching_files) > 20:
        result.append(f"\n... and {len(matching_files) - 20} more files ...")
    
    return "\n".join(result)

@tool
def list_directory(path: str = ".", recursive: bool = True) -> str:
    """
    List files and subdirectories under `path`. `path` must be inside one
    of the allow-listed root directories (see ALLOWED_ROOTS). Returns
    relative paths and file sizes so the model can decide what's worth
    reading with read_file. Call this before read_file when you don't
    already know exact filenames.
    """
    try:
        base = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"

    if not base.exists():
        return f"Error: '{base}' does not exist."
    if base.is_file():
        return f"'{base}' is a file, not a directory. Use read_file instead."

    entries = []
    walker = base.rglob("*") if recursive else base.glob("*")
    for p in sorted(walker):
        try:
            rel = p.relative_to(base)
        except ValueError:
            continue
        if p.is_dir():
            entries.append(f"[dir]  {rel}/")
        else:
            entries.append(f"[file] {rel}  ({p.stat().st_size:,} bytes)")
        if len(entries) >= MAX_LIST_ENTRIES:
            entries.append(f"... truncated at {MAX_LIST_ENTRIES} entries ...")
            break

    if not entries:
        return f"'{base}' is empty."
    return f"Contents of {base}:\n" + "\n".join(entries)


@tool
def read_file(path: str, max_chars: int = MAX_READ_CHARS) -> str:
    """
    Read and return the text content of a single file at `path` (must be
    inside an allow-listed root directory). Works for plain text, CSV/TSV,
    JSON, logs, markdown, etc. Binary files (images, .npy, compiled
    objects, ...) will likely come back as garbage -- don't use this for
    those. Output is truncated to `max_chars` characters.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"

    if not target.exists():
        return f"Error: '{target}' does not exist."
    if target.is_dir():
        return f"Error: '{target}' is a directory. Use list_directory instead."

    try:
        text = target.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return f"Error reading '{target}': {e}"

    if len(text) > max_chars:
        return (
            text[:max_chars]
            + f"\n\n... truncated, showing first {max_chars} of {len(text)} characters ..."
        )
    return text


# --- Optional on-demand indexing -----------------------------------------------

def _collection_dir_for(path: Path) -> Path:
    """Deterministic, filesystem-safe persist dir for an arbitrary indexed path."""
    digest = hashlib.sha1(str(path).encode()).hexdigest()[:16]
    return DYNAMIC_INDEX_ROOT / digest


@tool
def index_directory(path: str, force_rebuild: bool = False) -> str:
    """
    Build (or reuse) a searchable vector index over every supported file
    (.txt, .md, .pdf, .docx) under `path`, which must be inside an
    allow-listed root directory. Call this once before using
    search_indexed_directory on a new path. Re-run with
    force_rebuild=True if files have changed since the last index. This
    is for prose-style documents -- for structured data (CSV/TSV/JSON/
    logs), use list_directory + read_file instead, since there's nothing
    for embeddings to usefully chunk there.
    """
    try:
        base = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    if not base.is_dir():
        return f"Error: '{base}' is not a directory."

    persist_dir = _collection_dir_for(base)
    persist_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = persist_dir / "_manifest.json"

    current_manifest = {
        str(p.relative_to(base)): p.stat().st_mtime
        for p in sorted(base.rglob("*"))
        if p.is_file() and p.suffix.lower() in LOADERS_BY_EXT
    }

    stale = force_rebuild or not manifest_path.exists()
    if not stale:
        try:
            stale = json.loads(manifest_path.read_text()) != current_manifest
        except Exception:
            stale = True

    if not stale:
        return f"Index for '{base}' is already up to date at {persist_dir}."

    docs = []
    for p in sorted(base.rglob("*")):
        if not p.is_file() or p.suffix.lower() not in LOADERS_BY_EXT:
            continue
        loader_cls = LOADERS_BY_EXT[p.suffix.lower()]
        try:
            docs.extend(loader_cls(str(p)).load())
        except Exception as e:
            print(f"[fs_skill] Warning: failed to load {p}: {e}")

    if not docs:
        return (
            f"No supported documents found under '{base}' "
            f"(supported: {', '.join(LOADERS_BY_EXT)}). For CSV/TSV/JSON/log "
            f"output, use list_directory + read_file instead -- this index "
            f"is for prose-style documents."
        )

    splitter = RecursiveCharacterTextSplitter(chunk_size=1000, chunk_overlap=150)
    chunks = splitter.split_documents(docs)

    embeddings = build_embeddings()
    Chroma.from_documents(
        documents=chunks, embedding=embeddings, persist_directory=str(persist_dir)
    )
    manifest_path.write_text(json.dumps(current_manifest, indent=2))
    return f"Indexed {len(docs)} document(s) -> {len(chunks)} chunk(s) from '{base}'."


@tool
def search_indexed_directory(path: str, query: str, k: int = 4) -> str:
    """
    Search a directory previously indexed with index_directory. `path`
    must match the directory you passed to index_directory. Returns the
    most relevant excerpts with their source files.
    """
    try:
        base = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"

    persist_dir = _collection_dir_for(base)
    if not (persist_dir / "_manifest.json").exists():
        return f"'{base}' hasn't been indexed yet. Call index_directory on it first."

    embeddings = build_embeddings()
    store = Chroma(persist_directory=str(persist_dir), embedding_function=embeddings)
    results = store.similarity_search(query, k=k)
    if not results:
        return "No relevant content found."

    formatted = []
    for i, doc in enumerate(results, start=1):
        source = doc.metadata.get("source", "unknown source")
        formatted.append(f"[{i}] (source: {source})\n{doc.page_content}")
    return "\n\n".join(formatted)
