#!/usr/bin/env python3
import gzip
import os, argparse, json
from collections import defaultdict

import pandas as pd

from pmotools.extract_from_pmo.PMOWriter import PMOWriter
from pmotools.utils.small_utils import Utils
from pmotools.extract_from_pmo.PMOReader import PMOReader



def parse_args_combine_pmos():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pmo_files', type=str, required=True, help='a list of PMO files to combine into 1 PMO file, must be from same amplicon panel')
    parser.add_argument('--output', type=str, required=True, help='Output new combined PMO file')
    parser.add_argument('--overwrite', action = 'store_true', help='If output file exists, overwrite it')

    return parser.parse_args()

def combine_pmos():
    args = parse_args_combine_pmos()

    if args.output.endswith('.gz'):
        args.output = Utils.appendStrAsNeeded(args.output, ".json.gz")
    else:
        args.output = Utils.appendStrAsNeeded(args.output, ".json")
    Utils.outputfile_check(args.output, args.overwrite)

    pmo_files_list = Utils.parse_delimited_input_or_file(args.pmo_files, ",")
    if len(pmo_files_list) < 2:
        raise Exception("Only supplied " + str(len(pmo_files_list)) + " but multiple PMO files were expected")
    pmos = PMOReader.read_in_pmos(pmo_files_list)

    pmo_out = {}
    # create new output name
    analysis_name = ""
    for pmo in pmos:
        if "" != analysis_name:
            analysis_name += ";"
        analysis_name += pmo["analysis_name"]
    pmo_out["analysis_name"] = analysis_name
    # add in blank objects for default bases
    pmo_out["panel_info"] = {}
    pmo_out["experiment_infos"] = {}
    pmo_out["specimen_infos"] = {}
    pmo_out["sequencing_infos"] = {}
    pmo_out["microhaplotypes_detected"] = {}
    pmo_out["representative_microhaplotype_sequences"] = {}
    pmo_out["taramp_bioinformatics_infos"] = {}

    # have to check to make sure the specimen/experiment tables don't have repeated samples (or maybe add an option to check to see if they are the same)
    all_specimen_names = []
    duplicate_specimen_names = []
    for pmo in pmos:
        for specimen in pmo["specimen_infos"].keys():
            if specimen in all_specimen_names:
                duplicate_specimen_names.append(specimen)
            else:
                all_specimen_names.append(specimen)
    all_experiment_names = []
    duplicate_experiment_names = []
    for pmo in pmos:
        for experiment in pmo["experiment_infos"].keys():
            if experiment in all_experiment_names:
                duplicate_experiment_names.append(experiment)
            else:
                all_experiment_names.append(experiment)
    warnings = []
    if len(duplicate_specimen_names) > 0:
        warnings.append("Duplicate specimen names were supplied for the following specimens: " + ",".join(
            duplicate_specimen_names))
    if len(duplicate_experiment_names) > 0:
        warnings.append(
            "Duplicate experiment names were supplied for the following experiments: " + ",".join(duplicate_experiment_names))

    # check if the panels are the same
    panel_info_names = []
    for pmo in pmos:
        if len(pmo["panel_info"]) > 1:
            warnings.append("currently only supporting combining PMOs with 1 panel info")
        panel_id = next(iter(pmo["panel_info"].values()))["panel_id"]
        if panel_id not in panel_info_names:
            panel_info_names.append(panel_id)
    if len(panel_info_names) > 1:
        warnings.append("currently only supporting combining PMOs with 1 panel info")
    if len(warnings) > 0:
        raise Exception("\n".join(warnings))

    # combine specimen and experiment infos
    for pmo in pmos:
        pmo_out["specimen_infos"].update(pmo["specimen_infos"])
        pmo_out["experiment_infos"].update(pmo["experiment_infos"])
    # add panel info
    pmo_out["panel_info"].update(pmos[0]["panel_info"])

    # combine the taramp_bioinformatics_infos and sequencing_infos
    for pmo in pmos:
        pmo_out["taramp_bioinformatics_infos"].update(pmo["taramp_bioinformatics_infos"])
        pmo_out["sequencing_infos"].update(pmo["sequencing_infos"])

    # combine representative microhaplotypes
    def default_int_value():
        return defaultdict(int)
    def default_str_value():
        return defaultdict(str)
    rep_seq_counts = defaultdict(default_int_value)
    all_representative_haps_ids = []

    for pmo in pmos:
        for rep_haps in pmo["representative_microhaplotype_sequences"].values():
            all_representative_haps_ids.append(rep_haps["representative_microhaplotype_id"])
            for tar in rep_haps["targets"].values():
                for seq in tar["seqs"].values():
                    rep_seq_counts[tar["target_id"]][seq["seq"]] +=1

    representative_haps_id = ";".join([str(key) for key in all_representative_haps_ids] )
    # create a key for the old IDs
    new_rep_seq_ids = defaultdict(default_str_value)
    rep_targets = dict()
    for tar_name, tar_seqs in rep_seq_counts.items():
        microhaplotypes = dict()
        tar_seq_count = 0
        for seq in tar_seqs.keys():
            new_name = tar_name + "." + str(tar_seq_count).zfill(len(str(len(tar_seqs))))
            new_rep_seq_ids[tar_name][seq] = new_name
            tar_seq_count += 1
            microhaplotypes[new_name] =  {
                "microhaplotype_id": new_name,
                "seq": seq
            }
        rep_targets[tar_name] = microhaplotypes

    pmo_out["representative_microhaplotype_sequences"] = {
        representative_haps_id : {
            "representative_microhaplotype_id" : representative_haps_id,
            "targets" : rep_targets
        }
    }
    # rename the microhaplotypes detected
    pmo_out["microhaplotypes_detected"] = {}


    for pmo in pmos:
        pmo_out["microhaplotypes_detected"].update(pmo["microhaplotypes_detected"])
        for microhaplotypes_detected_id,microhaplotypes_detected_val in pmo["microhaplotypes_detected"].items():
            old_rep_haps_id = microhaplotypes_detected_val["representative_microhaplotype_id"]
            microhaplotypes_detected_val["representative_microhaplotype_id"] = representative_haps_id
            for exp_samples in pmo_out["microhaplotypes_detected"][microhaplotypes_detected_id]["experiment_samples"].values():
                for tar_name, tar_haps in exp_samples["target_results"].items():
                    for micro_hap in tar_haps["microhaplotypes"]:
                        old_micro_id = micro_hap["microhaplotype_id"]
                        new_micro_id = new_rep_seq_ids[tar_name][pmo["representative_microhaplotype_sequences"][old_rep_haps_id]["targets"][tar_name]["seqs"][old_micro_id]["seq"]]
                        micro_hap["microhaplotype_id"] = new_micro_id

    # write
    PMOWriter.write_out_pmo(pmo_out, args.output, args.overwrite)


if __name__ == "__main__":
    combine_pmos()

