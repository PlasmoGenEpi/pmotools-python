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
    parser.add_argument('--meta_fields', type=str, required=False, help='Meta Fields if any to include')
    parser.add_argument('--microhap_fields', type=str, required=False, help='additional optional fields from the microhaplotype object to include')


    return parser.parse_args()

def extract_for_dcifer():
    args = parse_args_extract_for_dcifer()
    output_extension = ".tsv"
    if args.delim == "tab":
        args.delim = "\t"
    elif args.delim == "comma" or args.delim == ",":
        args.delim = ","
        output_extension = ".csv"

    allele_freq_output = f"{args.output_stub}_allele_freq{output_extension}"
    allele_per_sample_table = f"{args.output_stub}_allele_table{output_extension}"

    Utils.inputOutputFileCheck(args.file, allele_freq_output, args.overwrite)
    Utils.inputOutputFileCheck(args.file, allele_per_sample_table, args.overwrite)

    with open(args.file) as f:
        pmodata = json.load(f)

    allele_counts, allele_freqs, target_totals = extract_allele_counts_freq_from_pmo(pmodata, args.bioid)
    if args.meta_fields is not None:
        args.meta_fields = args.meta_fields.split(",")
    if args.microhap_fields is not None:
        args.microhap_fields = args.microhap_fields.split(",")

    write_allele_per_sample_table(pmodata, args.bioid, allele_per_sample_table, args.microhap_fields, args.meta_fields)
    write_allele_freq_output(allele_freq_output, allele_freqs)

def extract_allele_counts_freq_from_pmo(pmodata, bioid, select_samples : list[str] = None, select_targets : list[str] = None):
    """
    Extract the allele counts from pmo for the bioinfomratics id. Can optionally only count for a subset of samples and targets

    :param pmodata: the pmo data structure
    :param bioid: the bioinfomratics id
    :param select_samples: an optional list of samples to get the info from
    :param select_targets: an optional list of targets to only get the info from
    :return: 3 dictionaries, 1) allele counts dict of dict (key1=target, key2=microhap, value=allele counts), 2) allele freq dict of dict (key1=target, key2=microhap, value=allele freq), 3) target_totals dict  (key1=target, value=target total alleles)
    """
    # prepare output data with defaultdicts so values default to 0 when adding to them
    allele_counts = defaultdict(lambda: defaultdict(int))
    allele_freqs = defaultdict(lambda: defaultdict(float))
    target_totals = defaultdict(int)

    for sample, sample_data in pmodata["microhaplotypes_detected"][bioid]["experiment_samples"].items():
        # only count for a subset of samples, could be only samples from a specific country or meta subfield
        if select_samples is not None and sample in select_samples:
            continue
        for target, target_data in sample_data["target_results"].items():
            # only count for specific sub set of targets
            if select_targets is not None and target not in select_targets:
                continue
            for microhapid in target_data["microhaplotypes"]:
                allele_id = microhapid["microhaplotype_id"]
                allele_counts[target][allele_id] += 1
                target_totals[target] += 1

    for target in allele_counts:
        for microhapid in allele_counts[target]:
            allele_freqs[target][microhapid] = allele_counts[target][microhapid] / target_totals[target]

    return allele_counts, allele_freqs, target_totals

def write_allele_per_sample_table(pmodata, bioid : str, output_file: str, additional_microhap_fields : list[str]  = None, meta_fields: list[str] = None):
    # check to see if at least 1 sample has supplied meta field
    # samples with missing field will get NA output
    # Check to see if at least 1 sample has supplied meta field
    if meta_fields is not None:
        # Find meta fields that have at least some data
        meta_fields_with_data = {
            metafield
            for metafield in meta_fields
            for specimen_data in pmodata["specimen_infos"].values()
            if metafield in specimen_data
        }

        # Determine meta fields with no samples having data
        meta_fields_with_no_samples = set(meta_fields) - meta_fields_with_data

        if meta_fields_with_no_samples:
            raise Exception(f"No samples have data for fields: {', '.join(meta_fields_with_no_samples)}")
    if additional_microhap_fields is not None:
        # Find meta fields that have at least some data
        additional_microhap_fields_with_data = {
            additional_microhap_field
            for additional_microhap_field in additional_microhap_fields
            for experiment_samples_data in pmodata["microhaplotypes_detected"][bioid]["experiment_samples"].values()
            for target_data in experiment_samples_data["target_results"].values()
            for microhap_data in target_data["microhaplotypes"]
            if additional_microhap_field in microhap_data
        }


        # Determine meta fields with no samples having data
        additional_microhap_fields_with_no_samples = set(additional_microhap_fields) - additional_microhap_fields_with_data

        if additional_microhap_fields_with_no_samples:
            raise Exception(f"No haplotypes have data for fields: {', '.join(additional_microhap_fields_with_no_samples)}")


    with open(output_file, mode='w') as f:
        # Writing target as locus as that is the default expectation for dcifer
        f.write("sampleID\tlocus\tallele")
        if additional_microhap_fields is not None:
            for additional_microhap_field in additional_microhap_fields:
                f.write(f"\t{additional_microhap_field}")
        if meta_fields is not None:
            for metafield in meta_fields:
                f.write(f"\t{metafield}")
        f.write("\n")

        specimen_infos = pmodata["specimen_infos"]
        experiment_infos = pmodata["experiment_infos"]
        detected_samples = pmodata["microhaplotypes_detected"][bioid]["experiment_samples"]

        for sample, sample_data in detected_samples.items():
            specimen_id = experiment_infos[sample]["specimen_id"]
            specimen_meta = specimen_infos[specimen_id]

            for target, target_data in sample_data["target_results"].items():
                for microhapid in target_data["microhaplotypes"]:
                    allele_id = microhapid["microhaplotype_id"]
                    f.write(f"{sample}\t{target}\t{allele_id}")
                    if additional_microhap_fields is not None:
                        for additional_microhap_field in additional_microhap_fields:
                            f.write(f"\t{microhapid.get(additional_microhap_field, 'NA')}")
                    if meta_fields is not None:
                        for metafield in meta_fields:
                            f.write(f"\t{specimen_meta.get(metafield, 'NA')}")
                    f.write("\n")

def write_allele_freq_output(output_file, allele_freqs):
    with open(output_file, mode='w') as f:
        # writing target as locus as that is the default expectation for dcifer
        f.write("locus\tallele\tfreq\n")
        for target, freq_data in allele_freqs.items():
            for microhapid, freq in freq_data.items():
                f.write(f"{target}\t{microhapid}\t{freq}\n")

if __name__ == "__main__":
    extract_for_dcifer()

