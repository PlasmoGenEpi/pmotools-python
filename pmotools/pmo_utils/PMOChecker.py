#!/usr/bin/env python3


import json


class PMOChecker:
    """
    A class to house utilities to help check the formatting of read in PMO files.
    """

    def __init__(self):
        self.all_required_base_fields = [
            "pmo_name",
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

        self.required_specimen_infos_fields = [
            "specimen_id",
            "samp_taxon_id",
            "collection_date",
            "collection_country",
            "collector",
            "samp_store_loc"
            "samp_collect_device",
            "project_name"
        ]

        self.required_experiment_infos_fields = [
            "experiment_sample_id",
            "sequencing_info_id",
            "panel_id",
            "specimen_id"
        ]

        self.required_sequencing_infos_fields = [
            "sequencing_info_id",
            "seq_instrument",
            "seq_date",
            "nucl_acid_ext",
            "nucl_acid_amp",
            "nucl_acid_ext_date",
            "nucl_acid_amp_date",
            "pcr_cond",
            "lib_screen",
            "lib_layout",
            "lib_kit"
        ]


    def check_for_required_base_fields(self, pmo_object):
        """
        Check that all required base fields are present in object

        :param pmo_object: the pmo object to check
        :return: return void if passes, otherwise raises an exception
        """
        missing_base_fields = []
        for base_field in self.all_required_base_fields:
            if base_field not in pmo_object:
                missing_base_fields.append(base_field)
        if len(missing_base_fields) > 0:
            raise Exception("Missing required base fields: {}".format(missing_base_fields))

    @staticmethod
    def check_bioinformatics_ids_consistency(pmo_object):
        """
        Check that all bio ids match

        :param pmo_object: the file to check
        :return: none, will raise exception if not consistent
        """
        warnings = []
        for bioid in pmo_object["taramp_bioinformatics_infos"].keys():
            if bioid not in pmo_object["microhaplotypes_detected"]:
                warnings.append("Missing " + bioid + " from " + "microhaplotypes_detected")
        if len(warnings) > 0:
            warnings.append("Available bioid(s) are {}".format(", ".join(pmo_object["microhaplotypes_detected"].keys())))
            raise Exception("\n".join(warnings))

    @staticmethod
    def check_for_bioinformatics_id(pmo_object, bioid):
        if bioid not in pmo_object["taramp_bioinformatics_infos"]:
            raise Exception(
                "Bioid ID {bioid} not found in PMO file, options are {avail_bioids}".format(bioid=bioid,
                                                                                            avail_bioids=",".join(
                                                                                                pmo_object[
                                                                                                    "taramp_bioinformatics_infos"].keys())))


