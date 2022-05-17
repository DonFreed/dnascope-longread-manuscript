#!/usr/bin/env python3

"""
Pull selected subsets, subtypes, filters, and columns from hap.py extended CSV output
"""

# Copyright (c) 2022 Sentieon Inc. All rights reserved

import argparse
import sys

import pandas as pd

def process_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--infiles", nargs="+", help="The input CSV")
    parser.add_argument("--sample_names", nargs="+", help="A sample name for each CSV")
    parser.add_argument("--subtypes", nargs="*", help="Subtypes to select")
    parser.add_argument("--subsets", nargs="*", help="Subsets to select")
    parser.add_argument("--filters", nargs="*", help="Filters to select")
    parser.add_argument("--columns", nargs="*", help="Columns to output")
    parser.add_argument("--outfile", default=sys.stdout, type=argparse.FileType('w'))
    return parser.parse_args(argv)

def main(args):
    assert len(args.infiles) == len(args.sample_names)

    # col_list
    col_list = []
    if hasattr(args, "columns") and args.columns:
        col_list = list(args.columns)

    # Read all input data
    dfs = []
    for infile, sample in zip(args.infiles, args.sample_names):
        sample_df = pd.read_csv(infile)
        sample_df["Sample"] = sample
        dfs.append(sample_df)
    df = pd.concat(dfs)

    # Add additional rows
    if col_list and ("TOTAL.ERRORS" in col_list or "ERRORS.PER.MB" in col_list):
        df["TOTAL.ERRORS"] = df["TRUTH.FN"] + df["QUERY.FP"]
        if "ERRORS.PER.MB" in col_list:
            df["ERRORS.PER.MB"] = df["TOTAL.ERRORS"] / df["Subset.IS_CONF.Size"] * 1e6

    # Filter the selected rows
    row_filters = []
    for arg_list, col_name in zip((args.subtypes, args.subsets, args.filters), ("Subtype", "Subset", "Filter")):
        if arg_list:
            row_filters.append(getattr(df, col_name).isin(arg_list))
        else:
            row_filters.append(pd.Series([True] * df[[col_name]].size, index=list(df.index)))
    row_filters = row_filters[0] & row_filters[1] & row_filters[2]
    df = df.loc[row_filters]

    # Filter columns
    if col_list:
        if "Sample" not in col_list:
            col_list.append("Sample")
        df = df[col_list]

    
    df = df.sort_values(by=["Subset", "Type", "Sample"])
    df.to_csv(args.outfile, sep='\t', index=False)

if __name__ == "__main__":
    args = process_args()
    main(args)
