"""
rds_skill.py
============

Lets the chat model read R serialized data (.rds) files -- e.g. output
from R packages like MicroTrait -- without requiring you to convert them
by hand first.

Two read strategies, tried in order:
  1. pyreadr (pure Python, no R installation needed) -- handles the
     common case: a single data.frame, or a handful of named ones.
  2. Rscript fallback, if R is installed and on PATH -- handles anything
     pyreadr can't (S4 objects, nested lists, etc.) by dumping either a
     CSV (if it can coerce the object to a data.frame) or an R str()
     summary otherwise.

Path safety: this reuses fs_skill's `_resolve_safe()` / ALLOWED_ROOTS
check directly rather than re-implementing it, so there is exactly one
place that decides what's reachable on disk.

Extra dependency: `pip install pyreadr`  (pandas comes along with it)
"""

import shutil
import subprocess
import tempfile
from pathlib import Path

from langchain_core.tools import tool

from .fs_skill import _resolve_safe  # shared allow-list check -- see module docstring

try:
    import pyreadr
    _HAVE_PYREADR = True
except ImportError:
    _HAVE_PYREADR = False

MAX_PREVIEW_ROWS = 50
MAX_PREVIEW_CHARS = 8000

# Minimal R script: readRDS the input, write a CSV if it's (or can become)
# a data.frame, otherwise dump an R str() summary so there's still
# something useful to read.
_R_CONVERTER = r"""
args <- commandArgs(trailingOnly = TRUE)
obj <- readRDS(args[1])
out <- args[2]
if (is.data.frame(obj)) {
  write.csv(obj, out, row.names = FALSE)
} else {
  ok <- tryCatch({
    write.csv(as.data.frame(obj), out, row.names = FALSE)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    sink(out)
    str(obj)
    sink()
  }
}
"""


def _format_dataframe_preview(df, name: str, max_rows: int, max_chars: int) -> str:
    lines = [
        f"--- {name} ---",
        f"shape: {df.shape[0]} rows x {df.shape[1]} columns",
        f"columns: {', '.join(str(c) for c in df.columns)}",
        "",
        df.head(max_rows).to_csv(index=False),
    ]
    text = "\n".join(lines)
    if len(text) > max_chars:
        text = text[:max_chars] + f"\n\n... truncated at {max_chars} characters ..."
    return text


@tool
def read_rds_file(
    path: str, max_rows: int = MAX_PREVIEW_ROWS, max_chars: int = MAX_PREVIEW_CHARS
) -> str:
    """
    Read an R .rds file and return a text summary: shape, column names,
    and a preview of the first `max_rows` rows (as CSV text). Works for
    .rds files containing one or more data.frames -- the common case for
    R analysis output such as MicroTrait results. `path` must be inside
    an allow-listed root directory (same rule as list_directory /
    read_file). Falls back to Rscript if the object isn't a plain
    data.frame and R is available on the host; otherwise reports what
    went wrong so you know whether to install pyreadr or check for R.
    """
    try:
        target = _resolve_safe(path)
    except ValueError as e:
        return f"Error: {e}"

    if not target.exists():
        return f"Error: '{target}' does not exist."
    if target.suffix.lower() != ".rds":
        return f"Error: '{target}' doesn't look like an .rds file."

    # --- Strategy 1: pyreadr, no R installation required ---
    if _HAVE_PYREADR:
        try:
            result = pyreadr.read_r(str(target))
            if not result:
                return f"'{target}' loaded but contained no readable objects."
            parts = []
            for name, df in result.items():
                label = name if name else target.name
                parts.append(_format_dataframe_preview(df, label, max_rows, max_chars))
            return "\n\n".join(parts)
        except Exception as e:
            pyreadr_error = str(e)
    else:
        pyreadr_error = "pyreadr is not installed (pip install pyreadr)."

    # --- Strategy 2: Rscript fallback, for objects pyreadr can't handle ---
    rscript = shutil.which("Rscript")
    if not rscript:
        return (
            f"Could not read '{target}' with pyreadr ({pyreadr_error}), and no "
            f"Rscript was found on PATH to fall back to. Install pyreadr "
            f"(`pip install pyreadr`) for plain data.frames, or make sure R/"
            f"Rscript is available on this host for everything else."
        )

    with tempfile.TemporaryDirectory() as tmp:
        script_path = Path(tmp) / "convert.R"
        out_path = Path(tmp) / "out.csv"
        script_path.write_text(_R_CONVERTER)
        try:
            subprocess.run(
                [rscript, str(script_path), str(target), str(out_path)],
                capture_output=True,
                text=True,
                timeout=120,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            return (
                f"pyreadr failed ({pyreadr_error}) and the Rscript fallback also "
                f"failed:\n{e.stderr.strip()}"
            )
        except subprocess.TimeoutExpired:
            return f"Rscript fallback timed out reading '{target}'."

        text = out_path.read_text(encoding="utf-8", errors="replace")
        if len(text) > max_chars:
            text = text[:max_chars] + f"\n\n... truncated at {max_chars} characters ..."
        return f"(read via Rscript fallback -- pyreadr could not parse this object)\n\n{text}"
