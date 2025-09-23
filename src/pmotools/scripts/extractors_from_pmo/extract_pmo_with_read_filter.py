#!/usr/bin/env python3
import argparse


from pmotools.pmo_engine.pmo_processor import PMOProcessor
from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.pmo_engine.pmo_writer import PMOWriter
from pmotools.utils.small_utils import Utils


def parse_args_extract_pmo_with_read_filter():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True, help="PMO file")
    parser.add_argument(
        "--output", type=str, required=True, help="Output json file path"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="If output file exists, overwrite it"
    )
    parser.add_argument(
        "--read_count_minimum",
        default=0.0,
        type=float,
        required=True,
        help="the minimum read count (inclusive) for detected haplotypes to be kept",
    )
    return parser.parse_args()


def extract_pmo_with_read_filter():
    args = parse_args_extract_pmo_with_read_filter()

    # check files
    Utils.inputOutputFileCheck(args.file, args.output, args.overwrite)

    # read in pmo
    pmo = PMOReader.read_in_pmo(args.file)

    # extract
    pmo_out = PMOProcessor.extract_from_pmo_with_read_filter(
        pmo, args.read_count_minimum
    )

    # write out the extracted
    args.output = PMOWriter.add_pmo_extension_as_needed(
        args.output, args.file.endswith(".gz") or args.output.endswith(".gz")
    )
    PMOWriter.write_out_pmo(pmo_out, args.output, args.overwrite)


if __name__ == "__main__":
    extract_pmo_with_read_filter()
