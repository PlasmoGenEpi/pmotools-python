#!/usr/bin/env python3


import json
from turtledemo.penrose import start


class PMOChecker:
    """
    A class to house utilites to help check the formatting of read in PMO files.
    """

    def __init__(self):
        self.all_required_base_fields = [
            "analysis_name",
            "panel_info",
            "experiment_infos",
            "specimen_infos",
            "sequencing_infos",
            "microhaplotypes_detected",
            "representative_microhaplotype_sequences",
            "taramp_bioinformatics_infos"
        ]

        self.optional_base_fields = [
            "target_demultiplexed_experiment_samples",
            "postprocessing_bioinformatics_infos"
        ]


    def check_for_required_base_fields(self, pmo_object):
        missing_base_fields = []
        for base_field in self.all_required_base_fields:
            if base_field not in pmo_object:
                missing_base_fields.append(base_field)
        if len(missing_base_fields) > 0:
            raise Exception("Missing required base fields: {}".format(missing_base_fields))

    @staticmethod
    def check_bioinformatics_ids(pmo_object):
        warnings = []
        for bioid in pmo_object["taramp_bioinformatics_infos"].keys():
            if bioid not in pmo_object["microhaplotypes_detected"]:
                warnings.append("Missing " + bioid + " from " + "microhaplotypes_detected")
            if bioid not in pmo_object["representative_microhaplotype_sequences"]:
                warnings.append("Missing " + bioid + " from " + "representative_microhaplotype_sequences")
        if len(warnings) > 0:
            raise Exception("\n".join(warnings))

    @staticmethod
    def check_for_bioinformatics_id(pmo_object, bioid):
        if bioid not in pmo_object["taramp_bioinformatics_infos"]:
            raise Exception(
                "Bioid ID {bioid} not found in PMO file, options are {avail_bioids}".format(bioid=bioid,
                                                                                            avail_bioids=",".join(
                                                                                                pmo_object[
                                                                                                    "taramp_bioinformatics_infos"].keys())))


