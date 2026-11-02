"""Expected input-read failures shared by harmonisation validation and reading."""

import lzma
import os
import shlex
import zipfile
import zlib


# Python's supported compression readers signal incomplete/corrupt streams
# through these exception types. In particular, gzip can raise zlib.error and
# LZMAError is not an OSError. Configuration ValueError must remain separate.
SUMMARY_STATISTICS_READ_ERRORS = (
    OSError, EOFError, lzma.LZMAError, zlib.error, zipfile.BadZipFile,
)


class SummaryStatisticsReadError(RuntimeError):
    """A supplied study stream cannot be read completely; retrying cannot fix it.

    The native cause is retained without inferring whether a transfer failed or
    the publisher supplied binary parts. The shared logger consumes these same
    presentation attributes; the batch runner marks this dataset failed.
    """

    def __init__(self, file_path, error):
        self.file_path = os.path.abspath(os.fspath(file_path))
        self.plain_english = (
            "Input file could not be read completely: %s. Read error (%s): %s"
            % (self.file_path, type(error).__name__, str(error) or "no details provided")
        )
        self.what_to_do = (
            "Check that the file is readable and the download/copy is complete. "
            "Compressed input may be incomplete, damaged, or only one part of a "
            "binary-split file. If the publisher supplied parts, combine all parts "
            "in the documented order before using the complete file."
        )
        if self.file_path.lower().endswith((".gz", ".bgz")):
            self.what_to_do += " Verify the complete gzip file with: gzip -t %s" % shlex.quote(self.file_path)
        super().__init__("%s %s" % (self.plain_english, self.what_to_do))
