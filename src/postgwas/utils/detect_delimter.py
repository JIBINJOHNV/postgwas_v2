from dataclasses import dataclass
from typing import Optional, List, Literal
import os, gzip, bz2, lzma, zipfile, csv
from pathlib import Path
import subprocess, shlex, os
from pathlib import Path 

# =========================================================
# RESULT OBJECT (safe + informative)
# =========================================================
@dataclass
class DelimiterDetectionResult:
    kind: Literal["delimiter", "whitespace", "unknown"]
    value: Optional[str] = None
    confidence: float = 0.0
    method: str = ""


# =========================================================
# FILE OPENER (compressed-safe)
# =========================================================
def _open_file(filepath: str):
    ext = Path(filepath).suffix.lower()
    if ext == ".gz":
        return gzip.open(filepath, "rt", errors="ignore")
    if ext == ".bz2":
        return bz2.open(filepath, "rt", errors="ignore")
    if ext in {".xz", ".lzma"}:
        return lzma.open(filepath, "rt", errors="ignore")
    if ext == ".zip":
        with zipfile.ZipFile(filepath) as zf:
            members = [f for f in zf.namelist() if not f.startswith("__MACOSX")]
            if len(members) != 1:
                raise ValueError(f"Zip must contain exactly one data file. Found: {members}")
            return (line.decode(errors="ignore") for line in zf.open(members[0]))
    return open(filepath, "rt", errors="ignore")


# =========================================================
# SAMPLE READER
# =========================================================
def _read_sample_lines(filepath, comment_char, sample_lines):
    lines = []
    with _open_file(filepath) as f:
        for line in f:
            if comment_char and line.lstrip().startswith(comment_char):
                continue
            if line.strip():
                lines.append(line.strip("\n"))
            if len(lines) >= sample_lines:
                break
    return lines


# =========================================================
# CONSISTENCY SCORER
# =========================================================
def _score_delimiter(lines, delim):
    counts = [len(l.split(delim)) for l in lines if delim in l]
    if not counts:
        return 0.0
    # consistency = low variance + multi-column
    unique = len(set(counts))
    avg = sum(counts) / len(counts)
    return avg / (1 + unique)


# =========================================================
# WHITESPACE DETECTOR (CRITICAL)
# =========================================================
def _detect_whitespace(lines):
    counts = [len(l.split()) for l in lines]
    if not counts:
        return False, 0.0
    unique = len(set(counts))
    avg = sum(counts) / len(counts)
    # strong signal: consistent multi-column split
    if avg > 2 and unique <= 2:
        return True, avg / (1 + unique)
    return False, 0.0


# =========================================================
# MAIN FUNCTION
# =========================================================
def detect_delimiter(
    filepath: str,
    comment_char: Optional[str] = "#",
    sample_lines: int = 50,
    candidates: Optional[List[str]] = None,
    allow_whitespace: bool = True,
    fail_policy: Literal["raise", "unknown"] = "raise",
) -> DelimiterDetectionResult:
    """
    Robust delimiter detection for general-purpose tabular files.
    Returns:
        DelimiterDetectionResult
    """
    if candidates is None:
        candidates = [",", "\t", ";", "|"]  # safer defaults
    if not os.path.exists(filepath):
        raise FileNotFoundError(filepath)
    lines = _read_sample_lines(filepath, comment_char, sample_lines)
    if not lines:
        if fail_policy == "raise":
            raise ValueError("No valid data lines found")
        return DelimiterDetectionResult("unknown")
    # -----------------------------------------------------
    # 1. Try CSV Sniffer (opportunistic)
    # -----------------------------------------------------
    try:
        sniffer = csv.Sniffer()
        dialect = sniffer.sniff("\n".join(lines), delimiters="".join(candidates))
        delim = dialect.delimiter
        score = _score_delimiter(lines, delim)
        if score > 1.5:
            return DelimiterDetectionResult(
                kind="delimiter",
                value=delim,
                confidence=score,
                method="sniffer",
            )
    except csv.Error:
        pass
    # -----------------------------------------------------
    # 2. Heuristic scoring for candidates
    # -----------------------------------------------------
    scores = {d: _score_delimiter(lines, d) for d in candidates}
    best_delim = max(scores, key=scores.get)
    if scores[best_delim] > 1.5:
        return DelimiterDetectionResult(
            kind="delimiter",
            value=best_delim,
            confidence=scores[best_delim],
            method="heuristic",
        )
    # -----------------------------------------------------
    # 3. Whitespace detection (IMPORTANT)
    # -----------------------------------------------------
    if allow_whitespace:
        is_ws, score = _detect_whitespace(lines)
        if is_ws:
            return DelimiterDetectionResult(
                kind="whitespace",
                value=None,
                confidence=score,
                method="whitespace",
            )
    # -----------------------------------------------------
    # 4. Fail
    # -----------------------------------------------------
    if fail_policy == "raise":
        raise ValueError("Could not determine delimiter")
    return DelimiterDetectionResult("unknown")



def normalize_whitespace_file(input_path: str, output_dir: str) -> str:
    base_name = Path(input_path).name
    while Path(base_name).suffix:
        base_name = Path(base_name).stem
    clean_path = os.path.join(output_dir, f"{base_name}_cleaned.tsv")
    # choose correct decompression
    if input_path.endswith(".gz"):
        cat_cmd = "gzip -dc"
    elif input_path.endswith(".xz") or input_path.endswith(".lzma"):
        cat_cmd = "xz -dc"
    else:
        cat_cmd = "cat"
    # normalize whitespace → tab
    cmd = (
        f"{cat_cmd} {shlex.quote(input_path)} | "
        f"grep -v '^##' | "
        f"tr -s '[:blank:]' '\\t' > {shlex.quote(clean_path)}"
    )
    subprocess.run(cmd, shell=True, check=True)
    return clean_path


def get_polars_separator(file, output_dir):
    res = detect_delimiter(file)
    # simple cases → direct
    if res.kind == "delimiter":
        return file, res.value
    # whitespace → normalize → tab
    elif res.kind == "whitespace":
        clean_file = normalize_whitespace_file(file, output_dir)
        return clean_file, "\t"
    else:
        raise ValueError(f"Could not detect delimiter: {file}")


# import csv
# import gzip
# import zipfile
# import bz2
# import lzma
# import io
# import os
# from collections import Counter
# from typing import Optional, List


# def detect_delimiter(
#     filepath: str,
#     comment_char: str = "#",
#     sample_lines: int = 50,
#     candidates: Optional[List[str]] = None,
# ) -> Optional[str]:
#     """
#     Detect the delimiter of a flat file (csv, tsv, txt), including
#     compressed variants (gzip, zip, bz2, xz/lzma).

#     Args:
#         filepath:     Path to the file.
#         comment_char: Lines starting with this character are ignored.
#                       Pass None to disable comment stripping entirely.
#         sample_lines: Number of non-comment lines to sample.
#         candidates:   Delimiters to consider if csv.Sniffer fails.

#     Returns:
#         The detected delimiter string, or None if it cannot be determined.

#     Raises:
#         FileNotFoundError: If the file does not exist.
#         ValueError:        If the file is a zip archive with multiple members.
#     """
#     if candidates is None:
#         candidates = [",", "\t", "|", ";", ":"]

#     if not os.path.exists(filepath):
#         raise FileNotFoundError(f"File not found: {filepath}")

#     # --- 1. Open the file (handles compression transparently) ---
#     raw_lines = _read_sample_lines(filepath, comment_char, sample_lines)

#     if not raw_lines:
#         return None

#     sample_text = "\n".join(raw_lines)

#     # --- 2. Try csv.Sniffer first (fast, handles quoted fields correctly) ---
#     sniffer = csv.Sniffer()
#     try:
#         dialect = sniffer.sniff(sample_text, delimiters="".join(candidates))
#         # Sniffer can guess wrongly on simple single-column files; verify it
#         # actually splits lines into more than one field.
#         if _sniffer_is_plausible(dialect.delimiter, raw_lines):
#             return dialect.delimiter
#     except csv.Error:
#         pass  # Sniffer failed — fall through to voting

#     # --- 3. Candidate voting: count occurrences per line, pick the winner ---
#     return _vote_delimiter(raw_lines, candidates)


# # ---------------------------------------------------------------------------
# # Internal helpers
# # ---------------------------------------------------------------------------

# def _open_file(filepath: str):
#     """
#     Return a text-mode file-like object regardless of compression format.
#     Supports: plain, .gz, .bz2, .xz/.lzma, .zip (single-member).
#     """
#     lower = filepath.lower()

#     if lower.endswith(".gz"):
#         return gzip.open(filepath, "rt", encoding="utf-8", errors="replace")

#     if lower.endswith(".bz2"):
#         return bz2.open(filepath, "rt", encoding="utf-8", errors="replace")

#     if lower.endswith((".xz", ".lzma")):
#         return lzma.open(filepath, "rt", encoding="utf-8", errors="replace")

#     if lower.endswith(".zip"):
#         zf = zipfile.ZipFile(filepath, "r")
#         members = zf.namelist()
#         # Filter out macOS metadata, directory entries, etc.
#         data_members = [
#             m for m in members
#             if not m.startswith("__MACOSX") and not m.endswith("/")
#         ]
#         if len(data_members) != 1:
#             raise ValueError(
#                 f"ZIP archive has {len(data_members)} data members; "
#                 "expected exactly 1. Pass the member path explicitly."
#             )
#         raw_bytes = zf.read(data_members[0])
#         zf.close()
#         return io.TextIOWrapper(
#             io.BytesIO(raw_bytes), encoding="utf-8", errors="replace"
#         )

#     # Plain text (csv / tsv / txt / no extension / anything else)
#     return open(filepath, "r", encoding="utf-8", errors="replace")


# def _read_sample_lines(
#     filepath: str, comment_char: Optional[str], max_lines: int
# ) -> List[str]:
#     """
#     Open the file, strip comment lines and blank lines, return up to
#     `max_lines` representative data rows.
#     """
#     lines: List[str] = []
#     fh = _open_file(filepath)
#     try:
#         for raw in fh:
#             stripped = raw.rstrip("\r\n")
#             if stripped.strip() == "":
#                 continue
#             if comment_char and stripped.lstrip().startswith(comment_char):
#                 continue
#             lines.append(stripped)
#             if len(lines) >= max_lines:
#                 break
#     finally:
#         fh.close()
#     return lines


# def _sniffer_is_plausible(delimiter: str, lines: List[str]) -> bool:
#     """
#     Return True when the sniffed delimiter actually produces >1 field
#     in the majority of sampled lines (guards against false positives).
#     """
#     if not lines:
#         return False
#     multi_field = sum(1 for ln in lines if delimiter in ln)
#     return multi_field / len(lines) >= 0.6


# def _vote_delimiter(lines: List[str], candidates: List[str]) -> Optional[str]:
#     """
#     For each candidate delimiter count how many lines it splits into
#     multiple fields *and* how consistent the field count is across lines.
#     The delimiter with the highest (consistent) score wins.
#     """
#     best_delim: Optional[str] = None
#     best_score: float = 0.0

#     for delim in candidates:
#         counts = [ln.count(delim) for ln in lines if delim in ln]
#         if not counts:
#             continue

#         coverage = len(counts) / len(lines)                    # fraction of lines that contain it
#         mode_count = Counter(counts).most_common(1)[0][0]
#         consistency = counts.count(mode_count) / len(counts)   # how regular the field count is
#         score = coverage * consistency * mode_count             # reward consistency + prevalence

#         if score > best_score:
#             best_score = score
#             best_delim = delim

#     return best_delim