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
        "--bioinfo_id",
        type=str,
        required=True,
        help="Identifier of bioinformatics processing run",
    )
    parser.add_argument(
        "--sampleID_col",
        type=str,
        default="sampleID",
        help="Column name containing sampleIDs",
    )
    parser.add_argument(
        "--locus_col",
        type=str,
        default="locus",
        help="Column name containing locus information",
    )
    parser.add_argument(
        "--mhap_col",
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

    args.output = Utils.appendStrAsNeeded(args.output, ".json")

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
        args.bioinfo_id,
        args.sampleID_col,
        args.locus_col,
        args.mhap_col,
        args.reads_col,
        addCols,
    )
    # Write output as json
    json_str = json.dumps(output_data, indent=4)
    with open(args.output, "w") as json_file:
        json_file.write(json_str)


if __name__ == "__main__":
    microhaplotype_table_to_json_file()
