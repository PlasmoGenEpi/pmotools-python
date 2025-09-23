#!/usr/bin/env python3

import argparse
import json
import os

from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.utils.small_utils import Utils
from pmotools.pmo_engine.pmo_checker import PMOChecker
from pmotools.pmo_engine.pmo_processor import PMOProcessor
from pmotools import __version__ as __pmotools_version__


def parse_args_extract_for_allele_table():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True, help="PMO file")
    parser.add_argument(
        "--jsonschema",
        default=os.path.join(
            os.path.dirname(
                os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            ),
            "schemas/",
            f"portable_microhaplotype_object_v{__pmotools_version__}.schema.json",
        ),
        type=str,
        required=False,
        help="the jsonschema to check the PMO against",
    )

    parser.add_argument(
        "--delim",
        default="tab",
        type=str,
        required=False,
        help="the delimiter of the input text file, examples tab,comma",
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Output allele table file name path"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="If output file exists, overwrite it"
    )
    parser.add_argument(
        "--allele_freqs_output",
        type=str,
        help="if also writing out allele frequencies, write to this file",
    )

    parser.add_argument(
        "--specimen_info_meta_fields",
        type=str,
        required=False,
        help="Meta Fields if any to include from the specimen table",
    )
    parser.add_argument(
        "--library_sample_info_meta_fields",
        type=str,
        required=False,
        help="Meta Fields if any to include from the library sample table",
    )
    parser.add_argument(
        "--microhap_fields",
        type=str,
        required=False,
        help="additional optional fields from the detected microhaplotype object to include",
    )
    parser.add_argument(
        "--representative_haps_fields",
        type=str,
        required=False,
        help="additional optional fields from the detected representative object to include",
    )
    parser.add_argument(
        "--default_base_col_names",
        type=str,
        required=False,
        default="library_sample_name,target_name,mhap_id",
        help="default base column names, must be length 3",
    )

    return parser.parse_args()


def extract_for_allele_table():
    args = parse_args_extract_for_allele_table()
    compressed_output = (
        "." not in args.output and args.file.endswith(".gz")
    ) or args.output.endswith(".gz")

    output_delim, output_extension = Utils.process_delimiter_and_output_extension(
        args.delim, gzip=compressed_output
    )

    allele_per_sample_table_out_fnp = (
        args.output
        if "STDOUT" == args.output
        else Utils.appendStrAsNeeded(args.output, output_extension)
    )
    Utils.inputOutputFileCheck(
        args.file, allele_per_sample_table_out_fnp, args.overwrite
    )

    allele_freq_output = ""
    if args.allele_freqs_output is not None:
        allele_freq_output = Utils.appendStrAsNeeded(
            args.allele_freqs_output, output_extension
        )
        Utils.inputOutputFileCheck(args.file, allele_freq_output, args.overwrite)

    pmodata = PMOReader.read_in_pmo(args.file)
    with open(args.jsonschema, "r") as f:
        schema_dict = json.load(f)
        checker = PMOChecker(schema_dict)
        # make sure PMO is valid
        checker.validate_pmo_json(pmodata)

    if args.specimen_info_meta_fields is not None:
        args.specimen_info_meta_fields = Utils.parse_delimited_input_or_file(
            args.specimen_info_meta_fields, ","
        )
    if args.microhap_fields is not None:
        args.microhap_fields = Utils.parse_delimited_input_or_file(
            args.microhap_fields, ","
        )
    if args.library_sample_info_meta_fields is not None:
        args.library_sample_info_meta_fields = Utils.parse_delimited_input_or_file(
            args.library_sample_info_meta_fields, ","
        )
    if args.representative_haps_fields is not None:
        args.representative_haps_fields = Utils.parse_delimited_input_or_file(
            args.representative_haps_fields, ","
        )

    allele_table = PMOProcessor.extract_alleles_per_sample_table(
        pmodata,
        additional_specimen_info_fields=args.specimen_info_meta_fields,
        additional_library_sample_info_fields=args.library_sample_info_meta_fields,
        additional_microhap_fields=args.microhap_fields,
        additional_representative_info_fields=args.representative_haps_fields,
        default_base_col_names=args.default_base_col_names.split(","),
    )
    with Utils.smart_open_write(allele_per_sample_table_out_fnp) as f:
        allele_table.to_csv(f, sep=output_delim, index=False)

    if args.allele_freqs_output is not None:
        allele_freqs = PMOProcessor.extract_allele_counts_freq_from_pmo(pmodata)
        with Utils.smart_open_write(allele_freq_output) as f:
            allele_freqs.to_csv(f, sep=output_delim, index=False)


if __name__ == "__main__":
    extract_for_allele_table()
