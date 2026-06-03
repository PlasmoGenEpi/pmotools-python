#!/usr/bin/env python3
import argparse

from pmotools.pmo_engine.pmo_exporter import PMOExporter
from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.utils.small_utils import Utils


def get_parser_export_pmo_into_xlsx() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pmotools-python export_pmo_into_xlsx",
        description="export all parts of a PMO into a .xlsx file",
    )
    parser.add_argument("--pmo", type=str, required=True, help="PMO file")
    parser.add_argument("--output", type=str, required=True, help="output file")
    parser.add_argument(
        "--overwrite", action="store_true", help="If output file exists, overwrite it"
    )
    return parser


def parse_args_export_pmo_into_xlsx():
    parser = get_parser_export_pmo_into_xlsx()
    return parser.parse_args()


def export_pmo_into_xlsx():
    args = parse_args_export_pmo_into_xlsx()

    # check files

    args.output = Utils.appendStrAsNeeded(args.output, ".xlsx")
    Utils.inputOutputFileCheck(args.pmo, args.output, args.overwrite)

    # read in PMO
    pmo = PMOReader.read_in_pmo(args.pmo)

    # export
    PMOExporter.export_to_excel(pmo, args.output)


if __name__ == "__main__":
    export_pmo_into_xlsx()
