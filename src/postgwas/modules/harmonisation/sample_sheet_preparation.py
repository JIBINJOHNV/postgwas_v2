"""Streaming, validated whitespace copies for sample-sheet preparation."""

from contextlib import contextmanager
from dataclasses import dataclass
import os
from pathlib import Path
import re
import tempfile

from postgwas.core.errors import ConfigurationError
from postgwas.core.io.delimiters import open_text
from postgwas.core.resource_preparation import sha256


class DelimiterPreparationError(ConfigurationError):
    """Normal delimiter detection failed; a configured rescue may be attempted."""


@dataclass(frozen=True)
class PreparedWhitespace:
    temporary: Path
    destination: Path
    data_rows: int
    changed_lines: int

    def publish(self) -> Path:
        """Reuse identical bytes or publish exclusively; never overwrite user data."""
        if self.destination.is_symlink():
            raise ConfigurationError("Converted destination must not be a symlink: %s" % self.destination)
        if self.destination.exists():
            if not self.destination.is_file() or sha256(self.destination) != sha256(self.temporary):
                raise ConfigurationError(
                    "Converted destination already exists with different content: %s. "
                    "Choose another preparation directory or move the existing file before rerunning."
                    % self.destination
                )
        else:
            os.link(self.temporary, self.destination)
        return self.destination.resolve()


@contextmanager
def whitespace_copy(source: Path, directory: Path, suffix: str, *, strip_double_hash: bool):
    """Apply tr '\\t' ' ' | tr -s ' ', validating every row before publication.

    Quote and empty-tab rejection are loss-prevention invariants: squeezing
    separators cannot safely preserve quoted text or positional empty fields.
    """
    if directory.is_symlink():
        raise ConfigurationError("Preparation directory must not be a symlink: %s" % directory)
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / (source.name + suffix)
    before = source.stat()
    descriptor, name = tempfile.mkstemp(prefix=".preparing-", dir=directory)
    temporary = Path(name)
    data_rows = changed_lines = 0
    width = None
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as output, open_text(source) as input_handle:
            for line_number, raw in enumerate(input_handle, 1):
                line = raw.rstrip("\r\n")
                converted = re.sub(" +", " ", line.replace("\t", " "))
                if line and not line.strip():
                    raise ConfigurationError("Line %d: whitespace-only row makes conversion ambiguous." % line_number)
                if line.strip() and not (strip_double_hash and line.startswith("##")):
                    if '"' in line or "\ufffd" in line or "\x00" in line:
                        raise ConfigurationError("Line %d: quoted or invalid text makes whitespace conversion ambiguous." % line_number)
                    if re.search(r"\t *\t", line) or not converted or converted.startswith(" ") or converted.endswith(" "):
                        raise ConfigurationError("Line %d: empty or boundary fields make whitespace conversion ambiguous." % line_number)
                    fields = converted.split(" ")
                    if width is None:
                        width = len(fields)
                    elif len(fields) != width:
                        raise ConfigurationError("Line %d: converted row has %d fields; expected %d." % (line_number, len(fields), width))
                    else:
                        data_rows += 1
                changed_lines += converted != line
                output.write(converted + "\n")
        if width is None or not data_rows:
            raise ConfigurationError("Whitespace conversion requires a header and at least one data row.")
        after = source.stat()
        if (before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_ino, after.st_size, after.st_mtime_ns):
            raise ConfigurationError("Input changed during whitespace conversion: %s" % source)
        yield PreparedWhitespace(temporary, destination, data_rows, changed_lines)
    finally:
        temporary.unlink(missing_ok=True)
