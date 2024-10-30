#!/usr/bin/env python3
import argparse

from pmotools.pmo_utils.PMOReader import PMOReader
from pmotools.utils.small_utils import Utils
from pmotools.pmo_utils.PMOChecker import PMOChecker
from pmotools.pmo_utils.PMOExtractor import PMOExtractor



def parse_args_extract_for_allele_table():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bioid', type=str, required=True, help='bio ID to extract for')
    parser.add_argument('--file', type=str, required=True, help='PMO file')
    parser.add_argument('--delim', default="tab", type=str, required=False, help='the delimiter of the input text file, examples tab,comma')
    parser.add_argument('--output', type=str, required=True, help='Output allele table file name path')
    parser.add_argument('--overwrite', action = 'store_true', help='If output file exists, overwrite it')
    parser.add_argument('--allele_freqs_output',type=str, help='if also writing out allele frequencies, write to this file')

    parser.add_argument('--specimen_info_meta_fields', type=str, required=False, help='Meta Fields if any to include from the specimen table')
    parser.add_argument('--experiment_info_meta_fields', type=str, required=False, help='Meta Fields if any to include from the experiment table')
    parser.add_argument('--microhap_fields', type=str, required=False, help='additional optional fields from the detected microhaplotype object to include')
    parser.add_argument('--representative_haps_fields', type=str, required=False, help='additional optional fields from the detected representative object to include')
    parser.add_argument('--default_base_col_names', type=str, required=False, default="sampleID,locus,allele", help='default base column names, must be length 3')


    return parser.parse_args()

def extract_for_allele_table():
    args = parse_args_extract_for_allele_table()

    output_delim, output_extension = Utils.process_delimiter_and_output_extension(args.delim, gzip=args.output.endswith('.gz'))

    allele_per_sample_table_out_fnp = args.output if "STDOUT" == args.output else Utils.appendStrAsNeeded(args.output, output_extension)
    Utils.inputOutputFileCheck(args.file, allele_per_sample_table_out_fnp, args.overwrite)

    allele_freq_output = ""
    if args.allele_freqs_output is not None:
        allele_freq_output = Utils.appendStrAsNeeded(args.allele_freqs_output, output_extension)
        Utils.inputOutputFileCheck(args.file, allele_freq_output, args.overwrite)

    checker = PMOChecker()
    pmodata = PMOReader.read_in_pmo(args.file)


    checker.check_for_required_base_fields(pmodata)
    checker.check_bioinformatics_ids_consistency(pmodata)
    checker.check_for_bioinformatics_id(pmodata, args.bioid)

    if args.specimen_info_meta_fields is not None:
        args.specimen_info_meta_fields = Utils.parse_delimited_input_or_file(args.specimen_info_meta_fields, ",")
    if args.microhap_fields is not None:
        args.microhap_fields = Utils.parse_delimited_input_or_file(args.microhap_fields, ",")
    if args.experiment_info_meta_fields is not None:
        args.experiment_info_meta_fields = Utils.parse_delimited_input_or_file(args.experiment_info_meta_fields, ",")
    if args.representative_haps_fields is not None:
        args.representative_haps_fields = Utils.parse_delimited_input_or_file(args.representative_haps_fields, ",")

    PMOExtractor.write_alleles_per_sample_table(pmodata, args.bioid, allele_per_sample_table_out_fnp,
                                                 additional_specimen_infos_fields = args.specimen_info_meta_fields,
                                                 additional_experiment_infos_fields = args.experiment_info_meta_fields,
                                                 additional_microhap_fields = args.microhap_fields,
                                                 additional_representative_infos_fields = args.representative_haps_fields,
                                                 output_delimiter=output_delim,
                                                 default_base_col_names=args.default_base_col_names.split(","))
    if args.allele_freqs_output is not None:
        allele_counts, allele_freqs, target_totals = PMOExtractor.extract_allele_counts_freq_from_pmo(pmodata, args.bioid)
        PMOExtractor.write_allele_freq_output(allele_freq_output, allele_freqs, output_delimiter=output_delim)

if __name__ == "__main__":
    extract_for_allele_table()

