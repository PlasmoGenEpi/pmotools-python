#!/usr/bin/env python3


import json
import gzip
import os
import sys
from collections import defaultdict

from pmotools.pmo_engine.pmo_checker import PMOChecker


class PMOReader:
    """
    A class for reading in PMO from files
    """
    @staticmethod
    def read_in_pmo(fnp : str | os.PathLike[str]):
        """
        Read in a PMO file, can either be compressed(.gz) or uncompressed
        :param fnp: the file name path of the PMO file to read in
        :return: a PMO like object
        """
        # check input
        checker = PMOChecker()
        if "STDIN" == fnp:
            pmo_data = json.load(sys.stdin)
        else:
            if fnp.endswith(".gz"):
                with gzip.open(fnp) as f:
                    pmo_data = json.load(f)
            else:
                with open(fnp) as f:
                    pmo_data = json.load(f)

        checker.check_for_required_base_fields(pmo_data)
        return pmo_data

    @staticmethod
    def read_in_pmos(fnps : list[str] | list[os.PathLike[str]]):
        """
        Read in a PMO file, can either be compressed(.gz) or uncompressed
        :param fnps: the file name path of the PMO file to read in
        :return: a list of PMO like object
        """
        ret = []
        for fnp in fnps:
            ret.append(PMOReader.read_in_pmo(fnp))
        return ret

    @staticmethod
    def combine_multiple_pmos(pmos : list[dict]):
        """
        Combine multiple PMOs into one pmo
        :param pmos: a list of PMO objects
        :return: a combined PMO
        """
        # create new pmo out
        pmo_out = {}
        # create new output name
        pmo_name = ""
        for pmo in pmos:
            if "" != pmo_name:
                pmo_name += ";"
            pmo_name += pmo["pmo_name"]
        pmo_out["pmo_name"] = pmo_name
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
                "Duplicate experiment names were supplied for the following experiments: " + ",".join(
                    duplicate_experiment_names))

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
            for tar in pmo["representative_microhaplotype_sequences"]["targets"].values():
                for seq in tar["seqs"].values():
                    rep_seq_counts[tar["target_id"]][seq["seq"]] += 1

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
                microhaplotypes[new_name] = {
                    "microhaplotype_id": new_name,
                    "seq": seq
                }
            rep_targets[tar_name] = microhaplotypes

        pmo_out["representative_microhaplotype_sequences"] = {
            "targets": rep_targets
        }
        # rename the microhaplotypes detected
        pmo_out["microhaplotypes_detected"] = {}

        for pmo in pmos:
            pmo_out["microhaplotypes_detected"].update(pmo["microhaplotypes_detected"])
            for microhaplotypes_detected_id, microhaplotypes_detected_val in pmo["microhaplotypes_detected"].items():

                for exp_samples in pmo_out["microhaplotypes_detected"][microhaplotypes_detected_id][
                    "experiment_samples"].values():
                    for tar_name, tar_haps in exp_samples["target_results"].items():
                        for micro_hap in tar_haps["microhaplotypes"]:
                            old_micro_id = micro_hap["microhaplotype_id"]
                            new_micro_id = new_rep_seq_ids[tar_name][
                                pmo["representative_microhaplotype_sequences"]["targets"][tar_name]["seqs"][old_micro_id]["seq"]]
                            micro_hap["microhaplotype_id"] = new_micro_id
        return pmo_out


