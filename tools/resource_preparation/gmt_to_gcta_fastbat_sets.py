#!/usr/bin/env python3
"""Standalone entry point for GCTA fastBAT GMT pathway preparation."""

from postgwas.modules.gcta_gene.pathway_sets import main


if __name__ == "__main__":
    raise SystemExit(main())
