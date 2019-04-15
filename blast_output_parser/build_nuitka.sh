#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

python -m nuitka --standalone --experimental=use_pefile --experimental=use_pefile_recurse blast_parser.py
