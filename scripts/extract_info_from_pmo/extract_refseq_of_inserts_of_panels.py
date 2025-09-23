#!/usr/bin/env python3
import argparse


from pmotools.pmo_engine.pmo_processor import PMOProcessor
from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.utils.small_utils import Utils


def parse_args_extract_refseq_of_inserts_of_panels():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True, help="PMO file")
    parser.add_argument(
        "--output", type=str, default="STDOUT", required=False, help="output file"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="If output file exists, overwrite it"
    )
    parser.description = "extract ref_seq of inserts of panels, but if no ref_seq is save in the PMO will just be blank"
    return parser.parse_args()


def extract_refseq_of_inserts_of_panels():
    args = parse_args_extract_refseq_of_inserts_of_panels()

    # check files
    Utils.inputOutputFileCheck(args.file, args.output, args.overwrite)

    # read in PMO
    pmo = PMOReader.read_in_pmo(args.file)

    # get panel insert locations
    panel_bed_locs = PMOProcessor.extract_panels_insert_bed_loc(pmo)

    # write
    with Utils.smart_open_write(args.output) as f:
        f.write("\t".join(["target_id", "ref_seq"]) + "\n")
        for loc in panel_bed_locs:
            f.write("\t".join([loc.name, loc.ref_seq]) + "\n")


if __name__ == "__main__":
    extract_refseq_of_inserts_of_panels()
