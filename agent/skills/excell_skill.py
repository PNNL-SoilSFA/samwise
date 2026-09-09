"""
excel_skill.py
==============

A standalone "Excel skill" for LangChain, built to handle workbooks of
any size — from a handful of rows up to files with millions of them.

It exposes three LangChain `@tool`s:

  * read_excel(path, sheet=None, start_row=1, max_rows=None)
        Overview of every sheet (call with just `path`), or a windowed,
        row-by-row dump of one sheet. Supports pagination via `start_row`
        so you can walk a large sheet in bounded pages.

  * find_in_excel(path, query, sheet=None, max_hits=None)
        Single streaming pass over the workbook that returns only the
        rows containing `query`. Constant memory, so it works on huge
        sheets where dumping everything is impossible.

  * export_excel_sheet(path, sheet, out_path=None, fmt="csv")
        Stream a whole sheet to a .csv/.tsv file on disk, one row at a
        time (constant memory). Use this when a sheet is too big to read
        into context — then point your normal text tools (or the RAG
        skill) at the exported file.

Design notes for large files:
  * Reading uses openpyxl in read-only/streaming mode, so rows are
    processed one at a time rather than loaded whole into memory.
  * Opening a large workbook is the expensive step (parsing shared
    strings), so the opened workbook is cached and reused across calls.
  * openpyxl re-parses a sheet from the top on each iteration, so
    jumping to a very large `start_row` gets progressively slower. For
    whole-sheet work on huge files, prefer find_in_excel /
    export_excel_sheet (one pass) over paging with big offsets.

Supported formats: .xlsx, .xlsm, .xltx, .xltm. The legacy .xls format is
NOT supported by openpyxl and needs a different reader (xlrd/pandas).

Dependencies:
    pip install openpyxl langchain-core

Configuration (via .env / environment):
    EXCEL_MAX_ROWS     - default rows per read_excel page      (default: 500)
    EXCEL_MAX_COLS     - max columns rendered per row          (default: 50)
    EXCEL_PREVIEW_ROWS - rows shown per sheet in overview mode (default: 10)
    EXCEL_MAX_HITS     - max matches returned by find_in_excel (default: 50)
    EXCEL_DATA_ONLY    - "1"/"true" to return cached formula VALUES,
                          "0"/"false" for formula strings       (default: true)
    EXCEL_DEEP_ROW_WARN- warn when start_row exceeds this       (default: 50000)
"""

import csv
import os
import re
from pathlib import Path
from typing import Optional

from langchain_core.tools import tool
from openpyxl import load_workbook

EXCEL_MAX_ROWS = int(os.getenv("EXCEL_MAX_ROWS", "1000000"))
EXCEL_MAX_COLS = int(os.getenv("EXCEL_MAX_COLS", "200"))
EXCEL_PREVIEW_ROWS = int(os.getenv("EXCEL_PREVIEW_ROWS", "50"))
EXCEL_MAX_HITS = int(os.getenv("EXCEL_MAX_HITS", "5000000"))
EXCEL_DATA_ONLY = os.getenv("EXCEL_DATA_ONLY", "true").strip().lower() in {"1", "true", "yes"}
EXCEL_DEEP_ROW_WARN = int(os.getenv("EXCEL_DEEP_ROW_WARN", "50000000"))

# openpyxl can read these; .xls (legacy) needs a different library.
SUPPORTED_EXTS = {".xlsx", ".xlsm", ".xltx", ".xltm"}


# --------------------------------------------------------------------------
# Workbook cache: opening a large workbook is slow (~seconds), so keep the
# most-recently-used one open and reuse it across tool calls. Size-1 cache;
# re-opens automatically if the file changes on disk.
# --------------------------------------------------------------------------
_WB_CACHE: dict = {}


def _get_workbook(path: Path):
    key = str(path.resolve())
    mtime = path.stat().st_mtime
    cached = _WB_CACHE.get(key)
    if cached is not None and cached[1] == mtime:
        return cached[0]
    # Evict (and close) anything else before opening the new one.
    for old_wb, _ in _WB_CACHE.values():
        try:
            old_wb.close()
        except Exception:
            pass
    _WB_CACHE.clear()
    wb = load_workbook(filename=str(path), read_only=True, data_only=EXCEL_DATA_ONLY)
    _WB_CACHE[key] = (wb, mtime)
    return wb


def _validate(path: str):
    """Return (Path, error_message). error_message is None when OK."""
    p = Path(path)
    if not p.exists():
        return p, f"Error: file not found: {path}"
    if p.suffix.lower() not in SUPPORTED_EXTS:
        return p, (
            f"Error: unsupported file type '{p.suffix}'. "
            f"Supported: {', '.join(sorted(SUPPORTED_EXTS))}. "
            f"(The legacy .xls format needs a different reader.)"
        )
    return p, None


def _fmt_cell(value) -> str:
    """Render a single cell as a one-line string (no tabs/newlines)."""
    if value is None:
        return ""
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def _dims(ws):
    """Stored (max_row, max_col); cheap in read-only mode, None if unsized."""
    try:
        return ws.max_row, ws.max_column
    except Exception:
        return None, None


def _resolve_sheet(wb, sheet: str):
    """Resolve a sheet by exact name, or by 0-based index given as a string."""
    if sheet in wb.sheetnames:
        return wb[sheet]
    s = str(sheet).strip()
    if s.lstrip("-").isdigit():
        idx = int(s)
        if 0 <= idx < len(wb.sheetnames):
            return wb[wb.sheetnames[idx]]
    return None


def _render_window(ws, start_row: int, max_rows: int, max_cols: int):
    """
    Render rows [start_row, start_row+max_rows) of a worksheet as text.
    Returns (text, rows_emitted, reached_end).
    """
    lines = []
    emitted = 0
    end_row = start_row + max_rows - 1
    for i, row in enumerate(
        ws.iter_rows(min_row=start_row, max_row=end_row, values_only=True)
    ):
        r_idx = start_row + i
        cells = list(row)
        cols_truncated = len(cells) > max_cols
        rendered = "\t".join(_fmt_cell(c) for c in cells[:max_cols])
        suffix = "  …(cols truncated)" if cols_truncated else ""
        lines.append(f"{r_idx:>7} | {rendered}{suffix}")
        emitted += 1
    reached_end = emitted < max_rows
    return "\n".join(lines), emitted, reached_end


def _overview(wb) -> str:
    """List every sheet with dimensions and a short preview of each."""
    parts = [
        f"Workbook contains {len(wb.sheetnames)} sheet(s): {', '.join(wb.sheetnames)}",
        "",
    ]
    for name in wb.sheetnames:
        ws = wb[name]
        max_row, max_col = _dims(ws)
        dims = (
            f"rows: {max_row if max_row is not None else 'unknown'}, "
            f"cols: {max_col if max_col is not None else 'unknown'}"
        )
        parts.append(f"=== Sheet: {name}  ({dims}) ===")
        text, emitted, reached_end = _render_window(
            ws, 1, EXCEL_PREVIEW_ROWS, EXCEL_MAX_COLS
        )
        parts.append(text if text else "(empty)")
        if not reached_end:
            more = f"{max_row - EXCEL_PREVIEW_ROWS}" if isinstance(max_row, int) else "more"
            parts.append(
                f"... {more} row(s) not shown. "
                f"Call read_excel(path, sheet='{name}') to page through them, "
                f"or find_in_excel / export_excel_sheet for the whole sheet."
            )
        parts.append("")
    return "\n".join(parts).rstrip()


@tool
def read_excel(
    path: str,
    sheet: Optional[str] = None,
    start_row: int = 1,
    max_rows: Optional[int] = None,
) -> str:
    """
    Read an Excel workbook (.xlsx/.xlsm/.xltx/.xltm) as plain text.

    Call with just `path` for an overview of every sheet (names,
    dimensions, and a short preview of each). Call with `sheet` (a sheet
    name, or a 0-based index like "0") to extract that sheet's contents
    row by row: one row per line, prefixed with its 1-based Excel row
    number, cells tab-separated.

    For large sheets, page through the data with `start_row` (1-based)
    and `max_rows` (page size). The response ends with a footer telling
    you the next start_row to request. Note: very large `start_row`
    values get progressively slower — to scan or dump an entire large
    sheet, use find_in_excel or export_excel_sheet instead.
    """
    p, err = _validate(path)
    if err:
        return err
    try:
        wb = _get_workbook(p)
    except Exception as e:
        return f"Error: failed to open {path}: {e}"

    if sheet is None:
        return _overview(wb)

    ws = _resolve_sheet(wb, sheet)
    if ws is None:
        return f"Error: sheet '{sheet}' not found. Available sheets: {', '.join(wb.sheetnames)}"

    if start_row < 1:
        start_row = 1
    page = max_rows if (max_rows and max_rows > 0) else EXCEL_MAX_ROWS

    text, emitted, reached_end = _render_window(ws, start_row, page, EXCEL_MAX_COLS)
    header = f"=== Sheet: {ws.title} (rows {start_row}–{start_row + emitted - 1}) ===" \
        if emitted else f"=== Sheet: {ws.title} ==="
    if not emitted:
        return f"{header}\n(no rows at or after row {start_row})"

    footer_bits = []
    if reached_end:
        footer_bits.append("[reached end of sheet]")
    else:
        next_row = start_row + emitted
        footer_bits.append(
            f"[more rows follow — next page: read_excel(path, sheet='{ws.title}', "
            f"start_row={next_row})]"
        )
    if start_row > EXCEL_DEEP_ROW_WARN:
        footer_bits.append(
            "[note: deep offsets are slow; prefer find_in_excel/export_excel_sheet "
            "for whole-sheet work]"
        )
    return f"{header}\n{text}\n" + "\n".join(footer_bits)


@tool
def find_in_excel(
    path: str,
    query: str,
    sheet: Optional[str] = None,
    max_hits: Optional[int] = None,
) -> str:
    """
    Search a workbook for rows containing `query` (case-insensitive
    substring match across all cells) and return just those rows. Makes a
    single streaming pass, so it works on very large sheets without
    loading them into memory.

    Searches every sheet by default, or restrict to one with `sheet` (a
    name or 0-based index). Results are capped (EXCEL_MAX_HITS); refine
    the query if you hit the cap.
    """
    p, err = _validate(path)
    if err:
        return err
    try:
        wb = _get_workbook(p)
    except Exception as e:
        return f"Error: failed to open {path}: {e}"

    if sheet is not None:
        ws = _resolve_sheet(wb, sheet)
        if ws is None:
            return f"Error: sheet '{sheet}' not found. Available sheets: {', '.join(wb.sheetnames)}"
        targets = [ws]
    else:
        targets = [wb[name] for name in wb.sheetnames]

    cap = max_hits if (max_hits and max_hits > 0) else EXCEL_MAX_HITS
    needle = query.lower()
    hits = []
    capped = False
    for ws in targets:
        for r_idx, row in enumerate(ws.iter_rows(values_only=True), start=1):
            cells = [_fmt_cell(c) for c in row]
            if needle in "\t".join(cells).lower():
                rendered = "\t".join(cells[:EXCEL_MAX_COLS])
                hits.append(f"[{ws.title}] row {r_idx}: {rendered}")
                if len(hits) >= cap:
                    capped = True
                    break
        if capped:
            break

    if not hits:
        return f"No rows containing '{query}' were found."
    header = f"Found {len(hits)} row(s) containing '{query}'" + (
        f" (capped at {cap}; refine your query for more):" if capped else ":"
    )
    return header + "\n" + "\n".join(hits)


@tool
def export_excel_sheet(
    path: str,
    sheet: str,
    out_path: Optional[str] = None,
    fmt: str = "csv",
) -> str:
    """
    Stream an entire sheet to a delimited text file on disk, one row at a
    time (constant memory) — the right way to handle a sheet too large to
    read into context. `sheet` is a name or 0-based index. `fmt` is "csv"
    or "tsv". Returns the output path and the number of rows written.

    After exporting, read the resulting file with your normal text tools
    or index it with the RAG skill.
    """
    p, err = _validate(path)
    if err:
        return err
    try:
        wb = _get_workbook(p)
    except Exception as e:
        return f"Error: failed to open {path}: {e}"

    ws = _resolve_sheet(wb, sheet)
    if ws is None:
        return f"Error: sheet '{sheet}' not found. Available sheets: {', '.join(wb.sheetnames)}"

    fmt = fmt.strip().lower()
    if fmt not in {"csv", "tsv"}:
        return f"Error: fmt must be 'csv' or 'tsv', got '{fmt}'."
    delimiter = "," if fmt == "csv" else "\t"

    if out_path:
        out = Path(out_path)
    else:
        safe = re.sub(r"[^\w.-]+", "_", ws.title).strip("_") or "sheet"
        out = p.parent / f"{p.stem}.{safe}.{fmt}"
    out.parent.mkdir(parents=True, exist_ok=True)

    rows_written = 0
    try:
        with out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter=delimiter)
            for row in ws.iter_rows(values_only=True):
                writer.writerow(["" if c is None else c for c in row])
                rows_written += 1
    except Exception as e:
        return f"Error: failed while exporting to {out}: {e}"

    return (
        f"Exported sheet '{ws.title}' -> {out} "
        f"({rows_written} row(s), {fmt.upper()}). "
        f"Read it with your text tools or index it with the RAG skill."
    )


if __name__ == "__main__":
    # Standalone usage:
    #   python excel_skill.py <workbook> [sheet] [start_row]
    # No sheet -> overview; sheet given -> paged dump from start_row.
    import sys

    target = sys.argv[1] if len(sys.argv) > 1 else os.getenv("EXCEL_TEST_FILE", "")
    if not target:
        print("Usage: python excel_skill.py <path-to-xlsx> [sheet] [start_row]")
        raise SystemExit(1)

    sheet_arg = sys.argv[2] if len(sys.argv) > 2 else None
    start_arg = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    print(read_excel.invoke({"path": target, "sheet": sheet_arg, "start_row": start_arg}))
