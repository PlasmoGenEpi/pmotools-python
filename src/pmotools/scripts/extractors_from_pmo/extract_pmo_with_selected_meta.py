#!/usr/bin/env python3
import argparse
import sys


from pmotools.pmo_engine.pmo_processor import PMOProcessor
from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.pmo_engine.pmo_writer import PMOWriter
from pmotools.utils.small_utils import Utils


def parse_args_extract_pmo_with_selected_meta():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True, help="PMO file")
    parser.add_argument(
        "--output", type=str, required=True, help="Output json file path"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="If output file exists, overwrite it"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="write out various messages about extraction",
    )
    parser.add_argument(
        "--metaFieldsValues",
        type=str,
        required=True,
        help="Meta Fields to include, should either be a table with columns field, values (and optionally group) or supplied command line as field1=value1,value2,value3:field2=value1,value2",
    )
    return parser.parse_args()


def extract_pmo_with_selected_meta():
    args = parse_args_extract_pmo_with_selected_meta()

    # check files
    Utils.inputOutputFileCheck(args.file, args.output, args.overwrite)

    # read in pmo
    pmo = PMOReader.read_in_pmo(args.file)

    # extract out of PMO
    pmo_out, group_counts = PMOProcessor.extract_from_pmo_samples_with_meta_groupings(
        pmo, args.metaFieldsValues
    )

    # write out the extracted
    args.output = PMOWriter.add_pmo_extension_as_needed(
        args.output, args.file.endswith(".gz") or args.output.endswith(".gz")
    )
    PMOWriter.write_out_pmo(pmo_out, args.output, args.overwrite)

    if args.verbose:
        sys.stdout.write(
            "Extracted the following number of specimens per group:" + "\n"
        )
        group_counts.to_csv(sys.stdout, sep="\t", index=True)


if __name__ == "__main__":
    extract_pmo_with_selected_meta()
