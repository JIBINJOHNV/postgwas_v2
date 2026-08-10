#!/usr/bin/env python3
"""Developer wrapper for ``postgwas resources prepare magma``."""

import sys

from postgwas.resources.cli import main


if __name__ == "__main__":
    raise SystemExit(main(["prepare", "magma", *sys.argv[1:]]))
