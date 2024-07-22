#!/usr/bin/env python3
import os, argparse, json
from collections import defaultdict

import pandas as pd
from pmotools.utils.small_utils import Utils




def parse_args_extract_for_dcifer():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bioid', type=str, required=True, help='bio ID to extract for')
    parser.add_argument('--file', type=str, required=True, help='PMO file')
    parser.add_argument('--delim', default="tab", type=str, required=False, help='the delimiter of the input text file, examples tab,comma')
    parser.add_argument('--output_stub', type=str, required=True, help='Output json file path')
    parser.add_argument('--overwrite', action = 'store_true', help='If output file exists, overwrite it')
    parser.add_argument('--metaFields', type=str, required=False, help='Meta Fields to include')

    return parser.parse_args()

def extract_for_dcifer():
    args = parse_args_extract_for_dcifer()
    outputExtension = ".tsv"
    if args.delim == "tab":
        args.delim = "\t"
        outputExtension = ".tsv"
    elif args.delim == "comma":
        args.delim = ","
        outputExtension = ".csv"
    allele_freq_output = args.output_stub + "_allele_freq" + outputExtension
    allele_per_sample_table = args.output_stub + "_allele_table" + outputExtension

    # make sure file exists
    Utils.inputOutputFileCheck(args.file, allele_freq_output, args.overwrite)
    Utils.inputOutputFileCheck(args.file, allele_per_sample_table, args.overwrite)

    pmodata = json.load(open(args.file))

    # print(pmodata.keys())
    # calculate allele freqs
    allele_freqs = {}
    allele_counts = {}
    locus_totals = defaultdict(int)

    for locus in pmodata["panel_info"]["targets"]:
        allele_counts[locus] = defaultdict(int)
        allele_freqs[locus] = defaultdict(float)
        locus_totals[locus] = 0





    with open(allele_per_sample_table, mode='w') as allele_per_sample_table_file:
        allele_per_sample_table_file.write("\t".join(["sampleID", "locus", "allele"]) + "\n")

        for sample in pmodata["microhaplotypes_detected"][args.bioid]["experiment_samples"]:
            for locus in pmodata["microhaplotypes_detected"][args.bioid]["experiment_samples"][sample]["target_results"]:
                for microhapid in pmodata["microhaplotypes_detected"][args.bioid]["experiment_samples"][sample]["target_results"][locus]["microhaplotypes"]:
                    # print(microhapid["microhaplotype_id"])
                    allele_counts[locus][microhapid["microhaplotype_id"]] = allele_counts[locus][microhapid["microhaplotype_id"]] + 1
                    locus_totals[locus] += 1
                    allele_per_sample_table_file.write("\t".join([sample, locus, microhapid["microhaplotype_id"]]) + "\n")

    for locus in allele_counts:
        for microhapid in allele_counts[locus]:
            allele_freqs[locus][microhapid] = allele_counts[locus][microhapid]/locus_totals[locus]

    with open(allele_freq_output, mode='w') as allele_freqs_file:
        allele_freqs_file.write("\t".join(["locus", "allele", "freq"]) + "\n")
        for locus in allele_freqs:
            for microhapid in allele_freqs[locus]:
                allele_freqs_file.write("\t".join([locus, microhapid, str(allele_freqs[locus][microhapid])]) + "\n")


if __name__ == "__main__":
    extract_for_dcifer()

