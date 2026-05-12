#!/usr/bin/env python3
import argparse
import json
import pandas as pd

from pmotools.pmo_builder.mhap_table_to_pmo import mhap_table_to_pmo
from pmotools.utils.small_utils import Utils


def parse_args_microhaplotype_table_to_json_file():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True, help="Input excel file path")
    parser.add_argument(
        "--bioinfo_name",
        type=str,
        required=False,
        help="Identifier of bioinformatics processing run",
    )
    parser.add_argument(
        "--library_sample_name_col",
        type=str,
        default="library_sample_name",
        help="Column name containing library_sample_name",
    )
    parser.add_argument(
        "--target_name_col",
        type=str,
        default="target_name",
        help="Column name containing target_name information",
    )
    parser.add_argument(
        "--seq_col",
        type=str,
        default="asv",
        help="Column name containing microhaplotypes",
    )
    parser.add_argument(
        "--reads_col",
        type=str,
        default="reads",
        help="Column name containing reads per microhaplotype",
    )
    parser.add_argument(
        "--additional_cols",
        type=str,
        help="Additional column name to add to detected haplotypes table, comma separated e.g. --additional_cols addCol,adddCol2",
    )
    parser.add_argument(
        "--delim", type=str, default="\t", help="Delimiter of input file"
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Output json file path"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="If output file exists, overwrite it"
    )
    return parser.parse_args()


def microhaplotype_table_to_json_file():
    args = parse_args_microhaplotype_table_to_json_file()

    ext = ".json.gz" if args.output.endswith(".json.gz") else ".json"
    args.output = Utils.appendStrAsNeeded(args.output, ext)

    addCols = None
    if args.additional_cols is not None:
        addCols = {}
        addColsToks = args.additional_cols.split(",")
        for addCol in addColsToks:
            if ":" in addCol:
                addColTok = addCol.split(":")
                if len(addColTok) == 2:
                    addCols[addColTok[0]] = addColTok[1]
                else:
                    raise Exception(
                        "should have only 1 :, found more than 1 while parsing: "
                        + addCol
                    )
            else:
                addCols[addCol] = addCol

    # check if input file exists and if output file exists check if --overwrite flag is set
    Utils.inputOutputFileCheckFromArgParse(args)

    contents = pd.read_csv(args.file, sep=args.delim)
    output_data = mhap_table_to_pmo(
        contents,
        args.bioinfo_name,
        args.library_sample_name_col,
        args.target_name_col,
        args.seq_col,
        args.reads_col,
        addCols,
    )
    # Write output as json
    json_str = json.dumps(output_data, indent=4)
    with Utils.smart_open_write(args.output) as json_file:
        json_file.write(json_str)


if __name__ == "__main__":
    microhaplotype_table_to_json_file()
