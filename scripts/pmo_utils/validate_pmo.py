#!/usr/bin/env python3
import gzip
import os, argparse, json
import jsonschema


from pmotools.pmo_engine.PMOReader import PMOReader
from pmotools.pmo_engine.pmo_checker import PMOChecker
from pmotools.utils.small_utils import Utils

def parse_args_validate_pmo():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pmo', type=str, required=True, help='a pmo file to validate')
    parser.add_argument('--jsconschema_file', default= os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                       "etc/portable_microhaplotype_object.schema.json"), type=str, required=False, help='jsonschema to validate against')

    return parser.parse_args()

def validate_pmo():
    args = parse_args_validate_pmo()

    # read in the PMO
    pmo = PMOReader.read_in_pmo(args.pmo)

    # create checker
    checker = PMOChecker(args.jsconschema_file)

    # validate
    checker.validate_pmo_json(pmo)


if __name__ == "__main__":
    validate_pmo()

