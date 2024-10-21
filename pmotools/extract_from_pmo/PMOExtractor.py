#!/usr/bin/env python3
import os

import pandas as pd
import json
from collections import defaultdict
from pmotools.utils.PMOChecker import PMOChecker


class PMOExtractor:

    @staticmethod
    def count_targets_per_sample(pmo, minimum_total_target_read_sum: float = 0.0):
        targets_per_sample_Count = defaultdict(lambda: defaultdict(int))
        for tar_amp_bioinformatics_info_id, results in pmo["microhaplotypes_detected"].items():

            for sample_id, sample in results["experiment_samples"].items():
                for target_id, targets in sample["target_results"].items():
                    read_count_sum = 0.0
                    for microhap in targets["microhaplotypes"]:
                        read_count_sum += microhap["read_count"]
                    if read_count_sum >= minimum_total_target_read_sum:
                        targets_per_sample_Count[tar_amp_bioinformatics_info_id][sample_id] += 1
        counts_df = pd.DataFrame(columns=["tar_amp_bioinformatics_info_id", "experiment_sample_id", "target_number"])
        for tar_amp_bioinformatics_info_id, sample_counts in targets_per_sample_Count.items():
            for sample_id, target_count in sample_counts.items():
                counts_df = pd.concat([counts_df, pd.DataFrame(
                    {"tar_amp_bioinformatics_info_id": tar_amp_bioinformatics_info_id,
                     "experiment_sample_id": sample_id,
                     "target_number": target_count}, index=[0])])
        return counts_df

    @staticmethod
    def count_samples_per_target(pmo, minimum_total_target_read_sum: float = 0.0):
        samples_per_target_count = defaultdict(lambda: defaultdict(int))
        for tar_amp_bioinformatics_info_id, results in pmo["microhaplotypes_detected"].items():

            for sample_id, sample in results["experiment_samples"].items():
                for target_id, targets in sample["target_results"].items():
                    read_count_sum = 0.0
                    for microhap in targets["microhaplotypes"]:
                        read_count_sum += microhap["read_count"]
                    if read_count_sum >= minimum_total_target_read_sum:
                        samples_per_target_count[tar_amp_bioinformatics_info_id][target_id] += 1
        counts_df = pd.DataFrame(columns=["tar_amp_bioinformatics_info_id", "target_id", "sample_number"])
        for tar_amp_bioinformatics_info_id, target_counts in samples_per_target_count.items():
            for target_id, target_count in target_counts.items():
                counts_df = pd.concat([counts_df, pd.DataFrame(
                    {"tar_amp_bioinformatics_info_id": tar_amp_bioinformatics_info_id,
                     "target_id": target_id,
                     "sample_number": target_count}, index=[0])])
        return counts_df

    @staticmethod
    def count_specimen_meta_subfields(pmo, meta_fields : list[str])-> pd.DataFrame:
        """
        Count the values of the meta fields, if a specimen doesn't have the field it will be entered as NA

        :param pmo: the pmo to count from
        :param meta_fields: the fields to get the counts for, if more than 1 field supplied will group the counts (e.g. will count field2 for all field1 values to get counts per sub-field gropsings)
        :return: counts for all sub-fields groups with column names being the meta fields plus the following: specimensCount, specimensFreq, totalSpecimenCount
        """
        field_counts = defaultdict(int)
        for specimen in pmo["specimen_infos"].values():
            value = ""
            for field in meta_fields:
                if "" != value:
                    value += "GOING_TO_SPLIT_ON_THIS_LATER"
                if field in specimen:
                    value += str(specimen[field])
                else:
                    value += "NA"
            field_counts[value] += 1

        counts_df = pd.DataFrame(columns=meta_fields + ["specimensCount", "specimensFreq", "totalSpecimenCount"])
        for field_name, field_count in field_counts.items():
            field_toks = field_name.split("GOING_TO_SPLIT_ON_THIS_LATER")
            new_row = {}
            for idx, meta_field in enumerate(meta_fields):
                new_row[meta_field] = field_toks[idx]
            new_row.update({"specimensCount": field_count,
                            "specimensFreq": field_count / len(pmo["specimen_infos"].values()),
                            "totalSpecimenCount": len(pmo["specimen_infos"].values())})
            counts_df.loc[len(counts_df)] = new_row
        counts_df.sort_values(by=meta_fields, inplace=True)
        return counts_df
    @staticmethod
    def count_specimen_meta_fields(pmo) -> pd.DataFrame:
        """
        Get a pandas dataframe of counts of the meta fields within the specimen_info section

        :param pmo: the pmo to count from
        :return: a pandas dataframe of counts with the following columns: field, presentInSpecimensCount, totalSpecimenCount
        """
        field_counts = defaultdict(int)
        for specimen in pmo["specimen_infos"].values():
            for meta_field in specimen:
                field_counts[meta_field] += 1
        counts_df = pd.DataFrame(columns=["field", "presentInSpecimensCount", "totalSpecimenCount"])
        for field_name, field_count in field_counts.items():
            counts_df.loc[len(counts_df)] = {"field": field_name, "presentInSpecimensCount": field_count,
                                             "totalSpecimenCount": len(pmo["specimen_infos"].values())}
        return counts_df

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

    @staticmethod
    def extract_from_pmo_select_specimens(pmo, specimen_ids: list[str]):
        # create a new pmo out
        # pmo_name, panel_info, sequencing_infos, taramp_bioinformatics_infos will stay the same
        # specimen_infos, experiment_infos, microhaplotypes_detected, representative_microhaplotype_sequences will be
        # created based on the supplied specimens

        # check to make sure the supplied specimens actually exist within the data
        warnings = []
        for specimen_id in specimen_ids:
            if specimen_id not in pmo["specimen_infos"]:
                warnings.append(f"{specimen_id} not in specimen_infos")
        if len(warnings) > 0:
            raise Exception("\n".join(warnings))

        pmo_out = {"pmo_name": pmo["pmo_name"],
                   "panel_info": pmo["panel_info"],
                   "sequencing_infos": pmo["sequencing_infos"],
                   "taramp_bioinformatics_infos": pmo["taramp_bioinformatics_infos"],
                   "specimen_infos": {},
                   "experiment_infos": {},
                   "microhaplotypes_detected": {},
                   "representative_microhaplotype_sequences": {}}

        # specimen_infos
        for specimen in pmo["specimen_infos"].values():
            if specimen["specimen_id"] in specimen_ids:
                pmo_out["specimen_infos"].update({specimen["specimen_id"]: specimen})
        # experiment_infos

        all_experiment_sample_ids = []
        for experiment in pmo["experiment_infos"].values():
            if experiment["specimen_id"] in specimen_ids:
                experiment_sample_id = experiment["experiment_sample_id"]
                all_experiment_sample_ids.append(experiment_sample_id)
                pmo_out["experiment_infos"].update({experiment_sample_id: experiment})

        # target_demultiplexed_experiment_samples
        if "target_demultiplexed_experiment_samples" in pmo:
            pmo_out["target_demultiplexed_experiment_samples"] = {}
            for bioid, demux_samples in pmo["target_demultiplexed_experiment_samples"].items():
                new_demux_samples = {"tar_amp_bioinformatics_info_id": bioid, "demultiplexed_experiment_samples" : {}}
                for exp_sample_id, sample in demux_samples.items():
                    if exp_sample_id in all_experiment_sample_ids:
                        new_demux_samples.update({exp_sample_id: sample})
                pmo_out["target_demultiplexed_experiment_samples"].update({bioid: new_demux_samples})

        # microhaplotypes_detected

        microhapids_for_tar = defaultdict(lambda: defaultdict(set))

        for id, microhaplotypes_detected in pmo["microhaplotypes_detected"].items():
            extracted_microhaps_for_id = {
                "tar_amp_bioinformatics_info_id": id,
                "representative_microhaplotype_id": microhaplotypes_detected["representative_microhaplotype_id"],
                "experiment_samples": {}}
            for experiment_sample_id, experiment in microhaplotypes_detected["experiment_samples"].items():
                if experiment_sample_id in all_experiment_sample_ids:
                    extracted_microhaps_for_id["experiment_samples"].update({experiment_sample_id: experiment})
                    for target_id, target in experiment["target_results"].items():
                        for microhap in target["microhaplotypes"]:
                            microhapids_for_tar[microhaplotypes_detected["representative_microhaplotype_id"]][target_id].add(microhap["microhaplotype_id"])
            pmo_out["microhaplotypes_detected"].update({id: extracted_microhaps_for_id})
        # representative_microhaplotype_sequences

        for rep_id, rep_haps in pmo["representative_microhaplotype_sequences"].items():
            rep_haps_for_id = {
                "representative_microhaplotype_id": rep_id,
                "targets": {}
            }
            for target_id, target in rep_haps["targets"].items():
                added_micro_haps = 0
                target_haps = {"target_id": target_id, "seqs": {}}
                for seq in target["seqs"].values():
                    if seq["microhaplotype_id"] in microhapids_for_tar[rep_id][target_id]:
                        target_haps["seqs"].update({seq["microhaplotype_id"]: seq})
                        added_micro_haps += 1
                if added_micro_haps > 0:
                    rep_haps_for_id["targets"].update({target_id: target_haps})
        return pmo_out

    @staticmethod
    def extract_from_pmo_select_targets(pmo, target_ids: list[str]):
        """
        Extract out data from the PMO for only select target IDs
        :param pmo: the pmo to extract data from
        :param target_ids: the target_ids to extract
        :return: a new pmo with the data for only the targets supplied
        """
        # create a new pmo out
        # pmo_name, panel_info, sequencing_infos, taramp_bioinformatics_infos will stay the same
        # specimen_infos, experiment_infos, microhaplotypes_detected, representative_microhaplotype_sequences will be
        # created based on the supplied specimens

        # check to make sure the supplied specimens actually exist within the data
        global target_id
        warnings = []
        for target_id in target_ids:
            found_in_rep_haps = False
            for rep_haps in pmo["representative_microhaplotype_sequences"].values():
                if target_id in rep_haps["targets"]:
                    found_in_rep_haps = True
                    break
            if not found_in_rep_haps:
                warnings.append(f"{target_id} not in representative_microhaplotype_sequences")
        if len(warnings) > 0:
            raise Exception("\n".join(warnings))

        pmo_out = {"pmo_name": pmo["pmo_name"],
                   "panel_info": pmo["panel_info"],
                   "sequencing_infos": pmo["sequencing_infos"],
                   "taramp_bioinformatics_infos": pmo["taramp_bioinformatics_infos"],
                   "specimen_infos": pmo["specimen_infos"],
                   "experiment_infos": pmo["experiment_infos"],
                   "microhaplotypes_detected": {},
                   "representative_microhaplotype_sequences": {}}

        # target_demultiplexed_experiment_samples
        if "target_demultiplexed_experiment_samples" in pmo:
            pmo_out["target_demultiplexed_experiment_samples"] = {}
            for bioid, demux_samples in pmo["target_demultiplexed_experiment_samples"].items():
                new_demux_samples = {"tar_amp_bioinformatics_info_id": bioid, "demultiplexed_experiment_samples" : {}}
                for sample_id, sample in demux_samples["demultiplexed_experiment_samples"].items():
                    targets_for_samples = {"experiment_sample_id": sample_id, "demultiplexed_targets": {}}
                    for target_id, target in sample.items():
                        if target_id in target_ids:
                            targets_for_samples["demultiplexed_targets"].update({target_id: target})
                    new_demux_samples.update({sample_id: targets_for_samples})
                pmo_out["target_demultiplexed_experiment_samples"].update({bioid: new_demux_samples})

        # microhaplotypes_detected
        for id, microhaplotypes_detected in pmo["microhaplotypes_detected"].items():
            extracted_microhaps_for_id = {
                "tar_amp_bioinformatics_info_id": id,
                "representative_microhaplotype_id": microhaplotypes_detected["representative_microhaplotype_id"],
                "experiment_samples": {}}
            for experiment_sample_id, experiment in microhaplotypes_detected["experiment_samples"].items():
                targets_for_samples = {"experiment_sample_id": experiment_sample_id, "target_results": {}}
                for target_id, target in experiment["target_results"].items():
                    if target_id in target_ids:
                        targets_for_samples["target_results"].update({target_id: target})
                extracted_microhaps_for_id["experiment_samples"].update({experiment_sample_id: targets_for_samples})
            pmo_out["microhaplotypes_detected"].update({id: extracted_microhaps_for_id})
        # representative_microhaplotype_sequences

        for rep_id, rep_haps in pmo["representative_microhaplotype_sequences"].items():
            rep_haps_for_id = {
                "representative_microhaplotype_id": rep_id,
                "targets": {}
            }
            for target_id, target in rep_haps["targets"].items():
                if target_id in target_ids:
                    rep_haps_for_id["targets"].update({target_id: target})
        return pmo_out

    @staticmethod
    def extract_from_pmo_samples_with_meta_groupings(pmo, meta_fields_values: str):
        """
        Extract out of a PMO the data associated with specimens that belong to specific meta data groupings
        :param pmo: the PMO to extract from
        :param meta_fields_values: Meta Fields to include, should either be a table with columns field, values (comma separated values) (and optionally group) or supplied command line as field1=value1,value2,value3:field2=value1,value2;field1=value5,value6, where each group is separated by a semicolon
        :return:
        """
        selected_meta_groups = {}
        # parse meta values
        if os.path.exists(meta_fields_values):
            selected_meta_groups = defaultdict(dict)
            meta_tab = pd.read_csv(meta_fields_values, sep='\t')
            if "field" not in meta_tab or "values" not in meta_tab:
                raise Exception(
                    meta_fields_values + " doesn't have columns field and values, has " + ",".join(meta_tab.columns))
            if "group" in meta_tab:
                for index, row in meta_tab.iterrows():
                    values_toks = row["values"].split(",")
                    selected_meta_groups[row["group"]][row["field"]] = values_toks
            else:
                group_criteria = {}
                for index, row in meta_tab.iterrows():
                    values_toks = row["values"].split(",")
                    group_criteria[row["field"]] = values_toks
                selected_meta_groups[0] = group_criteria
        else:
            group_toks = meta_fields_values.split(";")
            for idx, group_tok in enumerate(group_toks):
                group_criteria = {}
                field_with_values_toks = group_tok.split(":")
                for field_with_values_tok in field_with_values_toks:
                    field_values_toks = field_with_values_tok.split("=")
                    if len(field_values_toks) != 2:
                        raise Exception("error processing " + group_tok, " should be field and values separated by =")
                    values_toks = field_values_toks[1].split(",")
                    group_criteria[field_values_toks[0]] = values_toks
                selected_meta_groups[idx] = group_criteria

        # get count of fields
        fields_counts = PMOExtractor.count_specimen_meta_fields(pmo)

        # check to see if the fields supplied actually exit
        warnings = []
        fields_found = fields_counts["field"].tolist()
        for group in selected_meta_groups.values():
            for field in group.keys():
                if field not in fields_found:
                    warnings.append("missing the field: " + field + " in pmo")
        if len(warnings) > 0:
            raise Exception("\n".join(warnings))

        group_counts = defaultdict(int)
        all_specimen_ids = []
        for specimen in pmo["specimen_infos"].values():
            for group_name, meta in selected_meta_groups.items():
                passes_criteria = True
                for field, values in meta.items():
                    if not (field in specimen and str(specimen[field]) in values):
                        passes_criteria = False
                        break
                if passes_criteria:
                    group_counts[group_name] += 1
                    specimen_id = specimen["specimen_id"]
                    all_specimen_ids.append(specimen_id)

        pmo_out = PMOExtractor.extract_from_pmo_select_specimens(pmo, all_specimen_ids)

        return pmo_out, group_counts
