#!/usr/bin/env python
"""Check that the R and Python tracks wrote the same tables (numbers equal to a tolerance).

Usage: python tests/compare_tables.py DIR_R DIR_PY [file ...]
With no file names, every .tsv / .tsv.gz in DIR_R is compared.

Variant ID columns (lead_id, second_id) may differ where several variants are in perfect LD
and tie to rounding error; those rows pass when every numeric column still agrees.
"""
import os
import sys

import numpy as np
import pandas as pd

TOL = 1e-6


TIE_COLS = {"lead_id", "second_id"}


def compare(a, b):
    x = pd.read_csv(a, sep="\t")
    y = pd.read_csv(b, sep="\t")
    if list(x.columns) != list(y.columns) or x.shape != y.shape:
        return f"shape or columns differ: {x.shape} vs {y.shape}"
    worst = 0.0
    for c in x.columns:
        if pd.api.types.is_numeric_dtype(x[c]) and pd.api.types.is_numeric_dtype(y[c]):
            xv, yv = x[c].to_numpy(float), y[c].to_numpy(float)
            if c == "lead_r":  # a tied variant in perfect negative LD flips the sign
                xv, yv = np.abs(xv), np.abs(yv)
            d = np.nanmax(np.abs(xv - yv) / np.maximum(1, np.abs(xv))) if len(x) else 0
            if (x[c].isna() != y[c].isna()).any():
                return f"missing values differ in {c}"
            worst = max(worst, d)
        elif not (x[c].astype(str) == y[c].astype(str)).all():
            if c in TIE_COLS:
                n = (x[c].astype(str) != y[c].astype(str)).sum()
                print(f"     note: {c} differs in {n} of {len(x)} rows (tied variants)")
                continue
            return f"text column {c} differs"
    return None if worst <= TOL else f"max relative difference {worst:.2e}"


dir_r, dir_py = sys.argv[1:3]
files = sys.argv[3:] or sorted(f for f in os.listdir(dir_r) if f.endswith((".tsv", ".tsv.gz")))
bad = 0
for f in files:
    msg = compare(os.path.join(dir_r, f), os.path.join(dir_py, f))
    print(f"{'OK  ' if msg is None else 'FAIL'} {f}{'' if msg is None else ': ' + msg}")
    bad += msg is not None
sys.exit(1 if bad else 0)
