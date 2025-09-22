# lox-py

[![PyPI version](https://img.shields.io/pypi/v/lox-py)](https://pypi.org/project/lox-py/)

**lox-py** is a Python utility for processing paired-end sequencing LoxCode data.  
It extracts and counts LoxCode barcodes from FASTQ files, producing CSV summaries per sample.

---

## Features

- Align paired-end sequencing reads and extract LoxCode barcodes.
- Handles both `.fastq` and `.fastq.gz` files.
- Optional Biopython acceleration for FASTQ parsing.
- Approximate matching with edit distance (via `python-Levenshtein`).
- Fully Python refactor of Tom Weber’s original C++ scripts.

---

## Installation

```bash
pip install lox-py
