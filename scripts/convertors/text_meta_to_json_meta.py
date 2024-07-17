#!/usr/bin/env python3
import os, argparse, json
import pandas as pd


def parse_args_text_meta_to_json_meta():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, required=True, help='Input text file path')
    parser.add_argument('--delim', default="tab", type=str, required=False, help='the delimiter of the input text file, examples tab,comma')
    parser.add_argument('--index_col_name', type=str, required=False, help='by default output is a list, if an index column name is supplied it will be a dict with this column as index')
    parser.add_argument('--output', type=str, required=True, help='Output json file path')
    parser.add_argument('--overwrite', action = 'store_true', help='If output file exists, overwrite it')
    return parser.parse_args()

def text_meta_to_json_meta():
    args = parse_args_text_meta_to_json_meta()
    if args.delim == "tab":
        args.delim = "\t"
    elif args.delim == "comma":
        args.delim = ","
    if not args.output.endswith('.json'):
        args.output += '.json'
    # make sure file exists
    if not os.path.exists(args.file):
        raise FileNotFoundError(args.file)
    # only overwrite an existing file if --overwrite is on
    if os.path.exists(args.output) and not args.overwrite:
        raise Exception("Output file already exists, use --overWrite to overwrite it")
    index_col_name = None
    if args.index_col_name is not None:
        index_col_name = args.index_col_name
    contents = pd.read_csv(args.file, sep = args.delim, index_col=index_col_name)

    # Custom object_hook to replace None with an empty string
    def custom_object_hook(d):
        return {k: ("" if v is None else v) for k, v in d.items()}
    if args.index_col_name is not None:
        contents_json = json.loads(contents.to_json(orient="index",
                                                   index= True,
                                                    date_format="iso"),
                                   object_hook=custom_object_hook)
    else:
        contents_json = json.loads(contents.to_json(orient="records", date_format = "iso"), object_hook=custom_object_hook )

    json.dump(contents_json, open(args.output, 'w'), indent=4)

if __name__ == "__main__":
    text_meta_to_json_meta()

