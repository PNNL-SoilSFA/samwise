"""
large_file_skill.py
===================

Enables the chat agent to read large files completely by breaking them
into manageable chunks. Unlike read_file() in fs_skill.py which truncates
at MAX_READ_CHARS, this skill allows sequential reading of entire files.

The skill maintains reading state in a SQLite database, allowing:
  1. Reading files in chunks that fit within context limits
  2. Resuming reading across multiple tool calls
  3. Persisting progress across chat sessions
  4. Random access to specific portions of files

Use cases:
  - Analyzing large log files
  - Reading complete source code files
  - Processing large CSV/TSV data files
  - Examining full configuration files

Architecture:
  - Files are read in configurable chunks (default: 8000 chars)
  - Reading state (position, file hash) stored in SQLite
  - Supports both sequential reading and random access
  - Memory-mapped I/O for files > 10MB for efficiency

References:
  - Python mmap documentation: https://docs.python.org/3/library/mmap.html
  - SQLite for application state: https://www.sqlite.org/appfileformat.html
  - Chunked file processing patterns: https://realpython.com/read-write-files-python/

Configuration (via .env):
    LARGE_FILE_DB_PATH    - state database path (default: ./large_file_state.db)
    LARGE_FILE_CHUNK_SIZE - chars per chunk (default: 8000)
    LARGE_FILE_MAX_SIZE   - max file size to read in bytes (default: 100MB)
    LARGE_FILE_MMAP_THRESHOLD - use mmap above this size (default: 10MB)

Security: Uses fs_skill's ALLOWED_ROOTS check for path validation.
"""

import hashlib
import json
import mmap
import os
import sqlite3
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

from langchain_core.tools import tool

from .fs_skill import _resolve_safe, ALLOWED_ROOTS

# --- Configuration -----------------------------------------------------------

LARGE_FILE_DB_PATH = Path(os.getenv("LARGE_FILE_DB_PATH", "./large_file_state.db"))
LARGE_FILE_CHUNK_SIZE = int(os.getenv("LARGE_FILE_CHUNK_SIZE", "8000"))
LARGE_FILE_MAX_SIZE = int(os.getenv("LARGE_FILE_MAX_SIZE", str(100 * 1024 * 1024)))  # 100MB
LARGE_FILE_MMAP_THRESHOLD = int(os.getenv("LARGE_FILE_MMAP_THRESHOLD", str(10 * 1024 * 1024)))  # 10MB

# File extensions that are safe to read as text
# Reference: Common text file extensions
TEXT_EXTENSIONS = {
    # Code
    ".py", ".js", ".ts", ".java", ".c", ".cpp", ".h", ".hpp", ".cs", ".go",
    ".rs", ".rb", ".php", ".pl", ".pm", ".r", ".R", ".scala", ".kt", ".swift",
    ".nf", ".groovy", ".gradle", ".sh", ".bash", ".zsh", ".fish", ".ps1",
    # Data
    ".json", ".xml", ".yaml", ".yml", ".toml", ".ini", ".cfg", ".conf",
    ".csv", ".tsv", ".txt", ".log", ".md", ".rst", ".tex",
    # Web
    ".html", ".htm", ".css", ".scss", ".sass", ".less", ".vue", ".jsx", ".tsx",
    # Config
    ".env", ".gitignore", ".dockerignore", ".editorconfig",
    ".config", ".properties", ".plist",
    # Documentation
    ".markdown", ".adoc", ".org", ".wiki",
    # Other
    ".sql", ".graphql", ".proto", ".thrift", ".avsc",
    ".makefile", ".cmake", ".dockerfile",
}


# --- Database Setup ----------------------------------------------------------

def _init_database() -> sqlite3.Connection:
    """
    Initialize SQLite database for tracking file reading state.
    
    Schema stores:
      - File path and content hash (to detect changes)
      - Current reading position
      - Total file size
      - Session metadata
    
    Reference: https://www.sqlite.org/lang_createtable.html
    """
    LARGE_FILE_DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    
    conn = sqlite3.connect(str(LARGE_FILE_DB_PATH))
    conn.row_factory = sqlite3.Row
    
    conn.execute("""
        CREATE TABLE IF NOT EXISTS file_sessions (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            file_path TEXT NOT NULL,
            file_hash TEXT NOT NULL,
            file_size INTEGER NOT NULL,
            current_position INTEGER DEFAULT 0,
            chunk_size INTEGER NOT NULL,
            total_chunks INTEGER NOT NULL,
            chunks_read INTEGER DEFAULT 0,
            created_at TEXT NOT NULL,
            updated_at TEXT NOT NULL,
            completed INTEGER DEFAULT 0,
            metadata TEXT,
            UNIQUE(file_path, file_hash)
        )
    """)
    
    conn.execute("""
        CREATE INDEX IF NOT EXISTS idx_file_sessions_path 
        ON file_sessions(file_path)
    """)
    
    conn.commit()
    return conn


def _get_connection() -> sqlite3.Connection:
    """Get a database connection."""
    return _init_database()


# --- File Utilities ----------------------------------------------------------

def _compute_file_hash(path: Path, sample_size: int = 8192) -> str:
    """
    Compute a hash of the file for change detection.
    
    Uses first + last bytes plus file size for efficiency on large files.
    This is a common pattern for quick file identity checks.
    
    Reference: Similar approach used in rsync and git for quick comparisons
    https://git-scm.com/book/en/v2/Git-Internals-Git-Objects
    """
    file_size = path.stat().st_size
    hasher = hashlib.md5()
    
    # Include file size in hash
    hasher.update(str(file_size).encode())
    
    with open(path, "rb") as f:
        # Read first chunk
        hasher.update(f.read(sample_size))
        
        # Read last chunk if file is large enough
        if file_size > sample_size * 2:
            f.seek(-sample_size, 2)  # Seek from end
            hasher.update(f.read(sample_size))
    
    return hasher.hexdigest()


def _is_text_file(path: Path) -> bool:
    """
    Check if a file is likely a text file.
    
    Uses extension check first (fast), then samples content if needed.
    
    Reference: Standard heuristic used by tools like `file` command
    https://linux.die.net/man/1/file
    """
    # Check extension
    if path.suffix.lower() in TEXT_EXTENSIONS:
        return True
    
    # For unknown extensions, sample the file
    try:
        with open(path, "rb") as f:
            sample = f.read(8192)
            
        # Check for null bytes (binary indicator)
        if b"\x00" in sample:
            return False
        
        # Try to decode as UTF-8
        try:
            sample.decode("utf-8")
            return True
        except UnicodeDecodeError:
            # Try common encodings
            for encoding in ["latin-1", "cp1252"]:
                try:
                    sample.decode(encoding)
                    return True
                except UnicodeDecodeError:
                    continue
            return False
            
    except Exception:
        return False


def _read_chunk_standard(path: Path, start: int, size: int, encoding: str = "utf-8") -> tuple[str, int]:
    """
    Read a chunk from a file using standard file I/O.
    
    Returns (content, bytes_read).
    
    Reference: Python file I/O documentation
    https://docs.python.org/3/tutorial/inputoutput.html#reading-and-writing-files
    """
    with open(path, "r", encoding=encoding, errors="replace") as f:
        f.seek(start)
        content = f.read(size)
        actual_pos = f.tell()
    
    return content, actual_pos - start


def _read_chunk_mmap(path: Path, start: int, size: int, encoding: str = "utf-8") -> tuple[str, int]:
    """
    Read a chunk using memory-mapped I/O for large files.
    
    Memory mapping is more efficient for large files as it lets the OS
    handle paging and caching.
    
    Reference: Python mmap documentation
    https://docs.python.org/3/library/mmap.html
    """
    file_size = path.stat().st_size
    
    with open(path, "rb") as f:
        # Memory-map the file (read-only)
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            # Calculate byte positions (approximate for text)
            byte_start = min(start, file_size)
            byte_end = min(start + size * 4, file_size)  # *4 for UTF-8 worst case
            
            raw_bytes = mm[byte_start:byte_end]
            
    # Decode and trim to requested char count
    content = raw_bytes.decode(encoding, errors="replace")[:size]
    return content, len(content.encode(encoding, errors="replace"))


def _read_chunk(path: Path, start: int, size: int) -> tuple[str, int]:
    """
    Read a chunk from a file, choosing the appropriate method.
    
    Uses mmap for files above LARGE_FILE_MMAP_THRESHOLD.
    """
    file_size = path.stat().st_size
    
    if file_size > LARGE_FILE_MMAP_THRESHOLD:
        return _read_chunk_mmap(path, start, size)
    else:
        return _read_chunk_standard(path, start, size)


# --- Session Management ------------------------------------------------------

def _get_or_create_session(path: Path, chunk_size: int = LARGE_FILE_CHUNK_SIZE) -> dict:
    """
    Get existing session or create new one for a file.
    
    If the file has changed (different hash), creates a new session.
    """
    file_hash = _compute_file_hash(path)
    file_size = path.stat().st_size
    total_chunks = (file_size + chunk_size - 1) // chunk_size  # Ceiling division
    
    conn = _get_connection()
    
    # Check for existing session
    row = conn.execute(
        "SELECT * FROM file_sessions WHERE file_path = ? AND file_hash = ?",
        (str(path), file_hash)
    ).fetchone()
    
    if row:
        conn.close()
        return dict(row)
    
    # Create new session
    now = datetime.now().isoformat()
    cursor = conn.execute(
        """
        INSERT INTO file_sessions 
        (file_path, file_hash, file_size, current_position, chunk_size, 
         total_chunks, chunks_read, created_at, updated_at, completed)
        VALUES (?, ?, ?, 0, ?, ?, 0, ?, ?, 0)
        """,
        (str(path), file_hash, file_size, chunk_size, total_chunks, now, now)
    )
    session_id = cursor.lastrowid
    conn.commit()
    
    # Fetch and return the new session
    row = conn.execute(
        "SELECT * FROM file_sessions WHERE id = ?", (session_id,)
    ).fetchone()
    conn.close()
    
    return dict(row)


def _update_session(session_id: int, position: int, chunks_read: int, completed: bool = False) -> None:
    """Update session with new reading position."""
    conn = _get_connection()
    conn.execute(
        """
        UPDATE file_sessions 
        SET current_position = ?, chunks_read = ?, completed = ?, updated_at = ?
        WHERE id = ?
        """,
        (position, chunks_read, 1 if completed else 0, datetime.now().isoformat(), session_id)
    )
    conn.commit()
    conn.close()


def _reset_session(session_id: int) -> None:
    """Reset a session to start reading from the beginning."""
    conn = _get_connection()
    conn.execute(
        """
        UPDATE file_sessions 
        SET current_position = 0, chunks_read = 0, completed = 0, updated_at = ?
        WHERE id = ?
        """,
        (datetime.now().isoformat(), session_id)
    )
    conn.commit()
    conn.close()


# --- LangChain Tools ---------------------------------------------------------

@tool
def open_large_file(path: str, chunk_size: int = LARGE_FILE_CHUNK_SIZE) -> str:
    """
    Open a large file for reading and get its metadata. Call this FIRST
    before using read_next_chunk() to read the file contents.
    
    This creates a reading session that tracks your progress through the
    file, allowing you to read it in chunks across multiple tool calls.
    
    Args:
        path: Path to the file (must be in ALLOWED_ROOTS)
        chunk_size: Characters per chunk (default: 8000). Smaller chunks
                   are easier to process but require more calls.
    
    Returns:
        File metadata including size, number of chunks, and session info.
        If a previous incomplete session exists, it will be resumed.
    
    Example workflow:
        1. open_large_file("/path/to/large.log")
        2. read_next_chunk("/path/to/large.log")  # repeat until done
        3. Or: read_file_chunk("/path/to/large.log", chunk_number=5)  # random access
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    if target.is_dir():
        return f"Error: '{target}' is a directory. Use list_directory instead."
    
    # Check file size
    file_size = target.stat().st_size
    if file_size > LARGE_FILE_MAX_SIZE:
        return (
            f"Error: File is {file_size / (1024*1024):.1f}MB, exceeds maximum of "
            f"{LARGE_FILE_MAX_SIZE / (1024*1024):.0f}MB. Consider processing it externally."
        )
    
    # Check if it's a text file
    if not _is_text_file(target):
        return (
            f"Error: '{target}' appears to be a binary file. "
            f"This tool only supports text files."
        )
    
    # Get or create session
    try:
        session = _get_or_create_session(target, chunk_size)
    except Exception as e:
        return f"Error creating reading session: {e}"
    
    # Format response
    lines = [f"=== Opened: {target.name} ===\n"]
    lines.append(f"File size: {file_size:,} bytes ({file_size / 1024:.1f} KB)")
    lines.append(f"Chunk size: {session['chunk_size']:,} characters")
    lines.append(f"Total chunks: {session['total_chunks']}")
    
    if session['chunks_read'] > 0 and not session['completed']:
        lines.append(f"\n📍 Resuming from chunk {session['chunks_read'] + 1} of {session['total_chunks']}")
        lines.append(f"   Position: {session['current_position']:,} bytes")
        lines.append(f"   Progress: {session['chunks_read'] / session['total_chunks'] * 100:.1f}%")
    elif session['completed']:
        lines.append(f"\n✓ File was previously read completely.")
        lines.append(f"  Use reset_file_reading('{path}') to start over,")
        lines.append(f"  or read_file_chunk('{path}', chunk_number=N) for random access.")
    else:
        lines.append(f"\n📖 Ready to read. Use read_next_chunk('{path}') to begin.")
    
    lines.append(f"\nSession ID: {session['id']}")
    
    return "\n".join(lines)


@tool
def read_next_chunk(path: str) -> str:
    """
    Read the next chunk of a file that was opened with open_large_file().
    
    Call this repeatedly to read through the entire file. The reading
    position is automatically tracked, so you can stop and resume later.
    
    Args:
        path: Path to the file (must match what was passed to open_large_file)
    
    Returns:
        The next chunk of file content, plus progress information.
        When the file is complete, indicates that reading is finished.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    
    # Get session
    try:
        session = _get_or_create_session(target)
    except Exception as e:
        return f"Error: {e}"
    
    if session['completed']:
        return (
            f"✓ File '{target.name}' has been completely read.\n"
            f"  Total chunks: {session['total_chunks']}\n"
            f"  Use reset_file_reading('{path}') to read again, or\n"
            f"  read_file_chunk('{path}', chunk_number=N) for specific sections."
        )
    
    # Read the chunk
    try:
        content, bytes_read = _read_chunk(
            target, 
            session['current_position'], 
            session['chunk_size']
        )
    except Exception as e:
        return f"Error reading file: {e}"
    
    # Update session
    new_position = session['current_position'] + bytes_read
    new_chunks_read = session['chunks_read'] + 1
    is_complete = new_position >= session['file_size'] or len(content) < session['chunk_size']
    
    _update_session(session['id'], new_position, new_chunks_read, is_complete)
    
    # Format response
    chunk_num = new_chunks_read
    total = session['total_chunks']
    progress = new_chunks_read / total * 100
    
    header = f"--- Chunk {chunk_num}/{total} ({progress:.1f}%) ---"
    
    if is_complete:
        footer = f"\n--- END OF FILE ({target.name}) ---"
    else:
        footer = f"\n--- [Use read_next_chunk('{path}') for next chunk] ---"
    
    return f"{header}\n\n{content}\n{footer}"


@tool
def read_file_chunk(path: str, chunk_number: int) -> str:
    """
    Read a specific chunk from a file by chunk number (random access).
    
    Use this when you need to jump to a specific part of the file rather
    than reading sequentially. Chunk numbers start at 1.
    
    Args:
        path: Path to the file
        chunk_number: Which chunk to read (1-based)
    
    Returns:
        The requested chunk content.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    
    # Get session
    try:
        session = _get_or_create_session(target)
    except Exception as e:
        return f"Error: {e}"
    
    if chunk_number < 1 or chunk_number > session['total_chunks']:
        return (
            f"Error: Invalid chunk number {chunk_number}. "
            f"Valid range: 1-{session['total_chunks']}"
        )
    
    # Calculate position for this chunk
    start_position = (chunk_number - 1) * session['chunk_size']
    
    # Read the chunk
    try:
        content, bytes_read = _read_chunk(target, start_position, session['chunk_size'])
    except Exception as e:
        return f"Error reading file: {e}"
    
    # Format response
    header = f"--- Chunk {chunk_number}/{session['total_chunks']} (random access) ---"
    is_last = chunk_number == session['total_chunks']
    
    if is_last:
        footer = f"\n--- END OF FILE ({target.name}) ---"
    else:
        footer = f"\n--- [Chunk {chunk_number} of {session['total_chunks']}] ---"
    
    return f"{header}\n\n{content}\n{footer}"


@tool
def get_file_reading_status(path: str) -> str:
    """
    Check the reading status of a file - how much has been read, current
    position, etc.
    
    Args:
        path: Path to the file
    
    Returns:
        Reading progress and session information.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    
    # Check for existing session
    file_hash = _compute_file_hash(target)
    conn = _get_connection()
    row = conn.execute(
        "SELECT * FROM file_sessions WHERE file_path = ? AND file_hash = ?",
        (str(target), file_hash)
    ).fetchone()
    conn.close()
    
    if not row:
        file_size = target.stat().st_size
        return (
            f"File: {target.name}\n"
            f"Size: {file_size:,} bytes\n"
            f"Status: Not opened yet. Use open_large_file('{path}') to begin reading."
        )
    
    session = dict(row)
    progress = session['chunks_read'] / session['total_chunks'] * 100
    
    lines = [f"=== Reading Status: {target.name} ===\n"]
    lines.append(f"File size: {session['file_size']:,} bytes")
    lines.append(f"Chunk size: {session['chunk_size']:,} characters")
    lines.append(f"Total chunks: {session['total_chunks']}")
    lines.append(f"Chunks read: {session['chunks_read']}")
    lines.append(f"Progress: {progress:.1f}%")
    lines.append(f"Current position: {session['current_position']:,} bytes")
    
    if session['completed']:
        lines.append(f"\n✓ Status: COMPLETED")
    else:
        lines.append(f"\n📖 Status: IN PROGRESS")
        lines.append(f"   Next: read_next_chunk('{path}')")
    
    lines.append(f"\nSession started: {session['created_at'][:19]}")
    lines.append(f"Last updated: {session['updated_at'][:19]}")
    
    return "\n".join(lines)


@tool
def reset_file_reading(path: str) -> str:
    """
    Reset reading progress for a file, starting over from the beginning.
    
    Use this when you want to re-read a file from the start, or if the
    file content has changed and you need to start fresh.
    
    Args:
        path: Path to the file
    
    Returns:
        Confirmation that the reading session was reset.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    
    # Find and reset session
    file_hash = _compute_file_hash(target)
    conn = _get_connection()
    row = conn.execute(
        "SELECT id FROM file_sessions WHERE file_path = ? AND file_hash = ?",
        (str(target), file_hash)
    ).fetchone()
    conn.close()
    
    if not row:
        return f"No reading session found for '{target.name}'. Use open_large_file() to start."
    
    _reset_session(row['id'])
    
    return (
        f"✓ Reset reading session for '{target.name}'.\n"
        f"  Position: 0 bytes\n"
        f"  Use read_next_chunk('{path}') to begin reading from the start."
    )


@tool
def read_file_range(path: str, start_line: int = 1, end_line: int = 100) -> str:
    """
    Read a specific range of lines from a file. Useful when you know
    approximately where the content you need is located.
    
    Args:
        path: Path to the file
        start_line: First line to read (1-based, default: 1)
        end_line: Last line to read (inclusive, default: 100)
    
    Returns:
        The specified lines from the file.
    
    Note: For very large files, this may be slow as it needs to count
    lines from the beginning. For large files, consider using chunk-based
    reading instead.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    
    if start_line < 1:
        start_line = 1
    if end_line < start_line:
        return f"Error: end_line ({end_line}) must be >= start_line ({start_line})"
    
    max_lines = 500  # Safety limit
    if end_line - start_line + 1 > max_lines:
        return f"Error: Cannot read more than {max_lines} lines at once. Requested: {end_line - start_line + 1}"
    
    try:
        lines_to_return = []
        current_line = 0
        
        with open(target, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                current_line += 1
                
                if current_line < start_line:
                    continue
                if current_line > end_line:
                    break
                
                lines_to_return.append(f"{current_line:6d}: {line.rstrip()}")
        
        if not lines_to_return:
            return f"No lines found in range {start_line}-{end_line}. File has {current_line} lines."
        
        content = "\n".join(lines_to_return)
        
        header = f"--- Lines {start_line}-{min(end_line, current_line)} of {target.name} ---"
        
        return f"{header}\n\n{content}"
        
    except Exception as e:
        return f"Error reading file: {e}"


@tool
def search_in_large_file(path: str, search_term: str, context_lines: int = 2, max_matches: int = 20) -> str:
    """
    Search for a term in a large file and return matching lines with context.
    
    More efficient than reading the entire file when you're looking for
    specific content.
    
    Args:
        path: Path to the file
        search_term: Text to search for (case-insensitive)
        context_lines: Number of lines before/after each match to show (default: 2)
        max_matches: Maximum number of matches to return (default: 20)
    
    Returns:
        Matching lines with context and line numbers.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    
    if not search_term:
        return "Error: search_term is required."
    
    search_lower = search_term.lower()
    
    try:
        # Read all lines (with limit for very large files)
        all_lines = []
        with open(target, "r", encoding="utf-8", errors="replace") as f:
            for i, line in enumerate(f):
                all_lines.append(line.rstrip())
                if i > 100000:  # Safety limit: 100k lines
                    break
        
        # Find matches
        matches = []
        for i, line in enumerate(all_lines):
            if search_lower in line.lower():
                matches.append(i)
                if len(matches) >= max_matches:
                    break
        
        if not matches:
            return f"No matches found for '{search_term}' in '{target.name}'."
        
        # Format output with context
        lines = [f"Found {len(matches)} match(es) for '{search_term}' in {target.name}:\n"]
        
        for match_idx in matches:
            start = max(0, match_idx - context_lines)
            end = min(len(all_lines), match_idx + context_lines + 1)
            
            lines.append(f"--- Match at line {match_idx + 1} ---")
            for i in range(start, end):
                marker = ">>>" if i == match_idx else "   "
                lines.append(f"{marker} {i + 1:6d}: {all_lines[i]}")
            lines.append("")
        
        if len(matches) >= max_matches:
            lines.append(f"... (limited to first {max_matches} matches)")
        
        return "\n".join(lines)
        
    except Exception as e:
        return f"Error searching file: {e}"


@tool
def get_file_summary(path: str) -> str:
    """
    Get a quick summary of a large file without reading it completely:
    line count, file size, first few lines, last few lines.
    
    Use this to understand what's in a file before deciding how to read it.
    
    Args:
        path: Path to the file
    
    Returns:
        File summary including size, line count, and preview.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"
    
    if not target.exists():
        return f"Error: '{target}' does not exist."
    
    file_size = target.stat().st_size
    
    if not _is_text_file(target):
        return (
            f"File: {target.name}\n"
            f"Size: {file_size:,} bytes ({file_size / 1024:.1f} KB)\n"
            f"Type: Binary file (not readable as text)"
        )
    
    try:
        # Count lines and get first/last lines
        first_lines = []
        last_lines = []
        line_count = 0
        
        with open(target, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line_count += 1
                
                if line_count <= 10:
                    first_lines.append(line.rstrip())
                
                # Keep last 5 lines in a rotating buffer
                last_lines.append(line.rstrip())
                if len(last_lines) > 5:
                    last_lines.pop(0)
        
        lines = [f"=== File Summary: {target.name} ===\n"]
        lines.append(f"Size: {file_size:,} bytes ({file_size / 1024:.1f} KB)")
        lines.append(f"Lines: {line_count:,}")
        lines.append(f"Estimated chunks: {(file_size + LARGE_FILE_CHUNK_SIZE - 1) // LARGE_FILE_CHUNK_SIZE}")
        
        lines.append(f"\n--- First {len(first_lines)} lines ---")
        for i, line in enumerate(first_lines, 1):
            lines.append(f"  {i}: {line[:100]}{'...' if len(line) > 100 else ''}")
        
        if line_count > 15:
            lines.append(f"\n--- Last {len(last_lines)} lines ---")
            start_num = line_count - len(last_lines) + 1
            for i, line in enumerate(last_lines):
                lines.append(f"  {start_num + i}: {line[:100]}{'...' if len(line) > 100 else ''}")
        
        lines.append(f"\nTo read: open_large_file('{path}') then read_next_chunk('{path}')")
        
        return "\n".join(lines)
        
    except Exception as e:
        return f"Error analyzing file: {e}"


@tool
def list_reading_sessions() -> str:
    """
    List all active file reading sessions - files that have been opened
    but not completely read, or completed sessions.
    
    Returns:
        List of all reading sessions with their status.
    """
    conn = _get_connection()
    rows = conn.execute("""
        SELECT file_path, file_size, chunks_read, total_chunks, completed, updated_at
        FROM file_sessions
        ORDER BY updated_at DESC
        LIMIT 20
    """).fetchall()
    conn.close()
    
    if not rows:
        return "No reading sessions found. Use open_large_file() to start reading a file."
    
    lines = [f"File Reading Sessions ({len(rows)}):\n"]
    
    for row in rows:
        path = Path(row['file_path']).name
        progress = row['chunks_read'] / row['total_chunks'] * 100
        status = "✓ Complete" if row['completed'] else f"📖 {progress:.0f}%"
        size_kb = row['file_size'] / 1024
        
        lines.append(f"  {status} {path}")
        lines.append(f"       {size_kb:.1f}KB, {row['chunks_read']}/{row['total_chunks']} chunks")
        lines.append(f"       Last: {row['updated_at'][:16]}")
        lines.append("")
    
    return "\n".join(lines)
