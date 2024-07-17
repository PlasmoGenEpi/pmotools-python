#!/usr/bin/env python3
import os, argparse, json
import pandas as pd



def parse_args_extractor_microhapseq_with_selected_meta():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, required=True, help='PMO file')
    parser.add_argument('--delim', default="tab", type=str, required=False, help='the delimiter of the input text file, examples tab,comma')
    parser.add_argument('--output', type=str, required=True, help='Output json file path')
    parser.add_argument('--overwrite', action = 'store_true', help='If output file exists, overwrite it')
    parser.add_argument('--metaFields', type=str, required=False, help='Meta Fields to include')

    return parser.parse_args()

def extractor_microhapseq_with_selected_meta():
    args = parse_args_extractor_microhapseq_with_selected_meta()
    outputExtension = ".tsv"
    if args.delim == "tab":
        args.delim = "\t"
        outputExtension = ".tsv"
    elif args.delim == "comma":
        args.delim = ","
        outputExtension = ".csv"

    if not args.output.endswith(outputExtension):
        args.output += outputExtension
    # make sure file exists
    if not os.path.exists(args.file):
        raise FileNotFoundError(args.file)
    # only overwrite an existing file if --overwrite is on
    if os.path.exists(args.output) and not args.overwrite:
        raise Exception("Output file already exists, use --overWrite to overwrite it")
    pmodata = json.load(open(args.file))

    print(pmodata.keys())


if __name__ == "__main__":
    extractor_microhapseq_with_selected_meta()

