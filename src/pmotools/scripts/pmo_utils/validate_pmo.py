#!/usr/bin/env python3
import os
import argparse
import json


from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.pmo_engine.pmo_checker import PMOChecker
from pmotools.utils.small_utils import Utils
from pmotools import __version__ as __pmotools_version__


def parse_args_validate_pmo():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pmo", type=str, required=True, help="a pmo file to validate")
    parser.add_argument(
        "--jsonschema_file",
        default=os.path.join(
            os.path.dirname(
                os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            ),
            "schemas/",
            f"portable_microhaplotype_object_v{__pmotools_version__}.schema.json",
        ),
        type=str,
        required=False,
        help="jsonschema to validate against",
    )

    return parser.parse_args()


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
