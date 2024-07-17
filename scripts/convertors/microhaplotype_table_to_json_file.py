#!/usr/bin/env python3
import argparse
from sys import argv
import json
import os

from pmotools.microhaplotype_table_to_pmo_dict import microhaplotype_table_to_pmo_dict


def parse_args_excel_meta_to_json_meta():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, required=True,
                        help='Input excel file path')
    parser.add_argument('--bioinfo_id', type=str, required=True,
                        help='Identifier of bioinformatics processing run')
    parser.add_argument('--sampleID_col', type=str,
                        default='sampleID', help='Column name containing sampleIDs')
    parser.add_argument('--locus_col', type=str, default='locus',
                        help='Column name containing locus information')
    parser.add_argument('--mhap_col', type=str, default='asv',
                        help='Column name containing microhaplotypes')
    parser.add_argument('--reads_col', type=str, default='reads',
                        help='Column name containing reads per microhaplotype')
    parser.add_argument('--delim', type=str, default='\t',
                        help='Delimiter of input file')
    parser.add_argument('--output', type=str, required=True,
                        help='Output json file path')
    parser.add_argument('--overwrite', action='store_true',
                        help='If output file exists, overwrite it')
    return parser.parse_args()


def microhaplotype_table_to_json_file():
    args = parse_args_excel_meta_to_json_meta()
    if not args.output.endswith('.json'):
        args.output += '.json'
    # make sure file exists
    if not os.path.exists(args.file):
        raise FileNotFoundError(args.file)
    # only overwrite an existing file if --overwrite is on
    if os.path.exists(args.output) and not args.overwrite:
        raise Exception(
            "Output file already exists, use --overWrite to overwrite it")
    output_data = microhaplotype_table_to_pmo_dict(args.file, args.bioinfo_id, args.sampleID_col, args.locus_col,
                                                   args.mhap_col, args.reads_col, args.delim, args.output, args.overwrite)
    # Write output as json
    json_str = json.dumps(output_data, indent=4)
    with open(args.output, 'w') as json_file:
        json_file.write(json_str)


if __name__ == "__main__":
    microhaplotype_table_to_json_file()
