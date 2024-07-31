#!/usr/bin/env python3

import pandas as pd
import json
from collections import defaultdict
from pmotools.utils.PMOChecker import PMOChecker


class PMOExtractor:

    @staticmethod
    def extract_allele_counts_freq_from_pmo(pmodata, bioid, select_samples : list[str] = None, select_targets : list[str] = None):
        """
        Extract the allele counts from pmo for the bioinfomratics id. Can optionally only count for a subset of samples and targets

        :param pmodata: the pmo data structure
        :param bioid: the bioinfomratics id
        :param select_samples: an optional list of samples to get the info from
        :param select_targets: an optional list of targets to only get the info from
        :return: 3 dictionaries, 1) allele counts dict of dict (key1=target, key2=microhap, value=allele counts), 2) allele freq dict of dict (key1=target, key2=microhap, value=allele freq), 3) target_totals dict  (key1=target, value=target total alleles)
        """
        # check input
        checker = PMOChecker()
        checker.check_for_required_base_fields(pmodata)
        checker.check_bioinformatics_ids_consistency(pmodata)
        checker.check_for_bioinformatics_id(pmodata, bioid)

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

    @staticmethod
    def write_alleles_per_sample_table(pmodata,
                                       bioid: str,
                                       output_fnp: str,
                                       additional_specimen_infos_fields: list[str] = None,
                                       additional_experiment_infos_fields: list[str] = None,
                                       additional_microhap_fields: list[str] = None,
                                       additional_representative_infos_fields: list[str] = None,
                                       output_delimiter: str = '\t',
                                       default_base_col_names: list[str] = ["sampleID", "locus", "allele"]):
        """
        Write out sample, target and allele. Can optionally add on any other additional fields

        :param output_delimiter: the delimiter used to write the output file
        :param pmodata: the data to write from
        :param bioid: the bioinformatics ID for the analysis to write out
        :param output_fnp: the file name path of the output file
        :param additional_specimen_infos_fields: any additional fields to write from the specimen_infos object
        :param additional_experiment_infos_fields: any additional fields to write from the experiment_samples object
        :param additional_microhap_fields: any additional fields to write from the microhap object
        :param additional_representative_infos_fields: any additional fields to write from the representative_microhaplotype_sequences object
        :param default_base_col_names: The default column name for the sample, locus and allele
        :return: none, data will be written to file
        """

        # check input
        checker = PMOChecker()
        checker.check_for_required_base_fields(pmodata)
        checker.check_bioinformatics_ids_consistency(pmodata)
        checker.check_for_bioinformatics_id(pmodata, bioid)

        # Check to see if at least 1 sample has supplied meta field
        # samples without this meta field will have NA
        if additional_specimen_infos_fields is not None:
            # Find meta fields that have at least some data
            meta_fields_with_data = {
                metafield
                for metafield in additional_specimen_infos_fields
                for specimen_data in pmodata["specimen_infos"].values()
                if metafield in specimen_data
            }

            # Determine meta fields with no samples having data
            meta_fields_with_no_samples = set(additional_specimen_infos_fields) - meta_fields_with_data

            if meta_fields_with_no_samples:
                raise Exception(f"No specimen_infos have data for fields: {', '.join(meta_fields_with_no_samples)}")
        # Check to see if at least 1 sample has supplied meta field
        # samples without this meta field will have NA
        if additional_experiment_infos_fields is not None:
            # Find meta fields that have at least some data
            meta_fields_with_data = {
                metafield
                for metafield in additional_experiment_infos_fields
                for experiment_data in pmodata["experiment_infos"].values()
                if metafield in experiment_data
            }
            # Determine meta fields with no samples having data
            meta_fields_with_no_samples = set(additional_experiment_infos_fields) - meta_fields_with_data

            if meta_fields_with_no_samples:
                raise Exception(f"No experiment_infos have data for fields: {', '.join(meta_fields_with_no_samples)}")

        # Check to see if at least 1 haplotype has this field
        # samples without this meta field will have NA
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
            additional_microhap_fields_with_no_samples = set(
                additional_microhap_fields) - additional_microhap_fields_with_data

            if additional_microhap_fields_with_no_samples:
                raise Exception(
                    f"No microhaplotypes_detected have data for fields: {', '.join(additional_microhap_fields_with_no_samples)}")
        # Check to see if at least 1 haplotype has this field
        # samples without this meta field will have NA
        if additional_representative_infos_fields is not None:
            # Find meta fields that have at least some data
            additional_microhap_fields_with_data = {
                additional_microhap_field
                for additional_microhap_field in additional_representative_infos_fields
                for target_data in pmodata["representative_microhaplotype_sequences"][bioid]["targets"].values()
                for microhap_data in target_data["seqs"].values()
                if additional_microhap_field in microhap_data
            }
            # Determine meta fields with no samples having data
            additional_microhap_fields_with_no_samples = set(
                additional_representative_infos_fields) - additional_microhap_fields_with_data

            if additional_microhap_fields_with_no_samples:
                raise Exception(
                    f"No representative_microhaplotype_sequences have data for fields: {', '.join(additional_microhap_fields_with_no_samples)}")

        if len(default_base_col_names) != 3:
            raise Exception("Must have 3 default columns for allele counts, not {}".format(len(default_base_col_names)))

        with open(output_fnp, mode='w') as f:
            f.write(output_delimiter.join(default_base_col_names))
            if additional_experiment_infos_fields is not None:
                for experiment_infos_field in additional_experiment_infos_fields:
                    f.write(f"{output_delimiter}{experiment_infos_field}")
            if additional_specimen_infos_fields is not None:
                for metafield in additional_specimen_infos_fields:
                    f.write(f"{output_delimiter}{metafield}")
            if additional_microhap_fields is not None:
                for additional_microhap_field in additional_microhap_fields:
                    f.write(f"{output_delimiter}{additional_microhap_field}")
            if additional_representative_infos_fields is not None:
                for representative_infos_field in additional_representative_infos_fields:
                    f.write(f"{output_delimiter}{representative_infos_field}")
            f.write("\n")

            specimen_infos = pmodata["specimen_infos"]
            experiment_infos = pmodata["experiment_infos"]
            detected_samples = pmodata["microhaplotypes_detected"][bioid]["experiment_samples"]
            rep_haps = pmodata["representative_microhaplotype_sequences"][pmodata["microhaplotypes_detected"][bioid].get("representative_microhaplotype_id")]["targets"]

            for sample, sample_data in detected_samples.items():
                specimen_id = experiment_infos[sample]["specimen_id"]
                specimen_meta = specimen_infos[specimen_id]

                for target, target_data in sample_data["target_results"].items():
                    for microhapid in target_data["microhaplotypes"]:
                        allele_id = microhapid["microhaplotype_id"]

                        f.write(f"{sample}{output_delimiter}{target}{output_delimiter}{allele_id}")
                        if additional_experiment_infos_fields is not None:
                            for experiment_infos_field in additional_experiment_infos_fields:
                                f.write(f"{output_delimiter}{experiment_infos[sample].get(experiment_infos_field, 'NA')}")
                        if additional_specimen_infos_fields is not None:
                            for metafield in additional_specimen_infos_fields:
                                f.write(f"{output_delimiter}{specimen_meta.get(metafield, 'NA')}")
                        if additional_microhap_fields is not None:
                            for additional_microhap_field in additional_microhap_fields:
                                f.write(f"{output_delimiter}{microhapid.get(additional_microhap_field, 'NA')}")
                        if additional_representative_infos_fields is not None:
                            for representative_infos_field in additional_representative_infos_fields:
                                f.write(f"{output_delimiter}{rep_haps[target]["seqs"][allele_id].get(representative_infos_field, 'NA')}")
                        f.write("\n")

    @staticmethod
    def write_allele_freq_output(output_fnp : str, allele_freqs : dict, output_delimiter: str = "\t", col_names : list[str] = ["locus", "allele", "freq"]):
        """
        Write out the allele frequencies to a file

        :param output_fnp: the output file
        :param allele_freqs: the frequencies to write out
        :param output_delimiter: the output delimiter
        :param col_names: the column names for the output
        :return: none, data will be written to file
        """
        with open(output_fnp, mode='w') as f:
            # writing target as locus as that is the default expectation for dcifer
            f.write(output_delimiter.join(col_names) + "\n")
            for target, freq_data in allele_freqs.items():
                for microhapid, freq in freq_data.items():
                    f.write(f"{target}{output_delimiter}{microhapid}{output_delimiter}{freq}\n")
