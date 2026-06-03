#!/usr/bin/env python3
import os
import argparse
import json


from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.pmo_engine.pmo_checker import PMOChecker
from pmotools.utils.small_utils import Utils
from pmotools import __version__ as __pmotools_version__


def get_parser_validate_pmo() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pmotools-python validate_pmo",
        description="Validate a PMO file against a JSON Schema",
    )
    parser.add_argument("--pmo", type=str, required=True, help="a pmo file to validate")
    parser.add_argument(
        "--jsonschema_version",
        default=__pmotools_version__,
        type=str,
        required=False,
        help="version of the jsonschema to validate against (default: %(default)s)",
    )
    parser.add_argument(
        "--jsonschema_file",
        default=None,
        type=str,
        required=False,
        help="explicit jsonschema file path to validate against (overrides --jsonschema_version)",
    )
    return parser


def parse_args_validate_pmo():
    parser = get_parser_validate_pmo()
    args = parser.parse_args()

    # post-parse resolution stays here, NOT in get_parser
    if args.jsonschema_file is None:
        args.jsonschema_file = os.path.join(
            os.path.dirname(
                os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            ),
            "schemas/",
            f"portable_microhaplotype_object_v{args.jsonschema_version}.schema.json",
        )
    return args


def validate_pmo():
    args = parse_args_validate_pmo()

    # read in the PMO
    pmo = PMOReader.read_in_pmo(args.pmo)

    # create checker
    with Utils.smart_open_read_by_ext(args.jsonschema_file) as f:
        checker = PMOChecker(json.load(f))
        # validate
        checker.validate_pmo_json(pmo)


if __name__ == "__main__":
    validate_pmo()
