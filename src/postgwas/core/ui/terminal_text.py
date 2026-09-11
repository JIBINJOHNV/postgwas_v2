"""Bounded terminal-control decoding shared by display and durable capture."""

import codecs
import re


# ASCII terminal controls except the printable record separators TAB and LF.
_TERMINAL_CONTROLS = re.compile(r"[\x00-\x08\x0b-\x1f\x7f]")


class TerminalTextDecoder:
    """Incremental UTF-8/ANSI decoding, including split escape sequences.

    Native diagnostics retain their printable text. SGR, cursor movement and
    OSC strings (for example terminal hyperlinks) must not enter durable logs
    or override the PostGWAS theme. CR progress records become separate lines;
    CRLF remains a single newline. State is bounded even for malformed escapes.
    """

    def __init__(self, encoding):
        self.decoder = codecs.getincrementaldecoder(encoding)(errors="replace")
        self.state = "text"
        self.after_cr = False

    def decode(self, data, *, final=False):
        return self.clean(self.decoder.decode(data, final=final))

    def clean(self, decoded):
        """Strip terminal controls from already-decoded literal text."""
        if (self.state == "text" and not self.after_cr
                and _TERMINAL_CONTROLS.search(decoded) is None):
            return decoded
        output = []
        for char in decoded:
            # A malformed unterminated escape must not hide later diagnostic
            # lines. Terminal control strings do not own the next text record.
            if char in "\r\n":
                self.state = "text"
            if self.state == "escape":
                self.state = {
                    "[": "csi", "]": "string", "P": "string",
                    "^": "string", "_": "string",
                }.get(char, "text")
                continue
            if self.state == "csi":
                if "@" <= char <= "~":
                    self.state = "text"
                continue
            if self.state == "string":
                if char == "\x07":
                    self.state = "text"
                elif char == "\x1b":
                    self.state = "string_escape"
                continue
            if self.state == "string_escape":
                self.state = "text" if char == "\\" else "string"
                continue
            if char == "\x1b":
                self.state = "escape"
                continue
            if char == "\r":
                output.append("\n")
                self.after_cr = True
                continue
            if char == "\n" and self.after_cr:
                self.after_cr = False
                continue
            self.after_cr = False
            if char in "\n\t" or ord(char) >= 32 and char != "\x7f":
                output.append(char)
        return "".join(output)


