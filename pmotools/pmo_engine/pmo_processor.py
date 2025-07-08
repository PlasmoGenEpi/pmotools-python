#!/usr/bin/env python3
import contextlib
import gzip
import json
import os
import sys
from typing import NamedTuple

import pandas
import pandas as pd
from collections import defaultdict
from pmotools.pmo_engine.pmo_checker import PMOChecker

bed_loc_tuple = NamedTuple("bed_loc", [("chrom", str), ("start", int), ("end", int), ("name", str), ("score", float),
                                 ("strand", str), ("ref_seq", str), ("extra_info", str)])

class PMOProcessor:
    """
    A class to extract info out of a loaded PMO object
    """

    @staticmethod
    def get_index_key_of_specimen_names(pmodata):
        """
        Get key of specimen_name to index in pmodata["specimen_info"]
        :param pmodata: the PMO to get indexes from
        :return: a dictionary of indexes keyed by specimen_name
        """
        ret = {}
        for idx, specimen in enumerate(pmodata["specimen_info"]):
            ret[specimen["specimen_name"]] = idx
        return ret

    @staticmethod
    def get_index_key_of_experiment_sample_names(pmodata):
        """
        Get key of experiment_sample_name to index in pmodata["experiment_info"]
        :param pmodata: the PMO to get indexes from
        :return: a dictionary of indexes keyed by experiment_sample_name
        """
        ret = {}
        for idx, experiment in enumerate(pmodata["experiment_info"]):
            ret[experiment["experiment_sample_name"]] = idx
        return ret

    @staticmethod
    def get_index_key_of_target_names(pmodata):
        """
        Get key of target_name to index in pmodata["target_info"]
        :param pmodata: the PMO to get indexes from
        :return: a dictionary of indexes keyed by target_name
        """
        ret = {}
        for idx, target in enumerate(pmodata["target_info"]):
            ret[target["target_name"]] = idx
        return ret

    @staticmethod
    def get_index_key_of_panel_names(pmodata):
        """
        Get key of panel_name to index in pmodata["panel_info"]
        :param pmodata: the PMO to get indexes from
        :return: a dictionary of indexes keyed by panel_name
        """
        ret = {}
        for idx, panel in enumerate(pmodata["panel_info"]):
            ret[panel["panel_name"]] = idx
        return ret

    @staticmethod
    def get_index_key_of_target_in_microhaplotypes_info(pmodata):
        """
        Get key of target_name to index for the representative microhaplotypes for the target_name in pmodata["microhaplotypes_info"]
        :param pmodata: the PMO to get indexes from
        :return: a dictionary of indexes keyed by target_name
        """
        ret = {}
        for idx, microhaplotypes_info_for_target in enumerate(pmodata["microhaplotypes_info"]["targets"]):
            ret[pmodata["target_info"][microhaplotypes_info_for_target["target_id"]]["target_name"]] = idx
        return ret

    @staticmethod
    def get_index_of_specimen_names(pmodata, specimen_names: list[str]):
        """
        Get index of specimen_name in pmodata["specimen_info"]
        :param pmodata: the PMO to get indexes from
        :param specimen_names: a list of specimen_names
        :return: the index of specimen_names in pmodata["specimen_info"] returned in the same order as specimen_names
        """
        specimen_key = PMOProcessor.get_index_key_of_specimen_names(pmodata)
        return [specimen_key[name] for name in specimen_names]

    @staticmethod
    def get_index_of_experiment_sample_names(pmodata, experiment_sample_names: list[str]):
        """
        Get index of experiment_sample_name in pmodata["experiment_info"]
        :param pmodata: the PMO to get indexes from
        :param experiment_sample_names: a list of experiment_sample_names
        :return: the index of experiment_sample_names in pmodata["experiment_info"] returned in the same order as experiment_sample_names
        """
        experiment_sample_key = PMOProcessor.get_index_key_of_experiment_sample_names(pmodata)
        return [experiment_sample_key[name] for name in experiment_sample_names]

    @staticmethod
    def get_index_of_target_names(pmodata, target_names: list[str]):
        """
        Get index of target_name in pmodata["target_info"]
        :param pmodata: the PMO to get indexes from
        :param target_names: a list of target_names
        :return: the index of target_names in pmodata["target_info"] returned in the same order as target_names
        """
        target_key = PMOProcessor.get_index_key_of_target_names(pmodata)
        return [target_key[name] for name in target_names]

    @staticmethod
    def get_index_of_panel_names(pmodata, panel_names: list[str]):
        """
        Get index of panel_name in pmodata["panel_info"]
        :param pmodata: the PMO to get indexes from
        :param panel_names: a list of panel_names
        :return: the index of panel_names in pmodata["panel_info"] returned in the same order as panel_names
        """
        panel_key = PMOProcessor.get_index_key_of_panel_names(pmodata)
        return [panel_key[name] for name in panel_names]

    @staticmethod
    def get_index_of_target_in_microhaplotypes_info(pmodata, target_names: list[str]):
        """
        Get index of target_name in pmodata["microhaplotypes_info"]["targets"]
        :param pmodata: the PMO to get indexes from
        :param target_names: a list of target_names
        :return: the index of target_names in pmodata["microhaplotypes_info"]["targets"] returned in the same order as target_names
        """
        microhap_target_key = PMOProcessor.get_index_key_of_target_in_microhaplotypes_info(pmodata)
        return [microhap_target_key[name] for name in target_names]

    @staticmethod
    def get_experiment_ids_for_specimen_ids(pmodata, specimen_ids: set[int]):
        """
        get a dictionary that lists the experiment_ids for a specimen_id
        :param pmodata: the PMO to get indexes from
        :param specimen_ids: a set of specimen_ids
        :return: a dictionary that lists the experiment_ids for a specimen_id
        """
        ret = defaultdict(set)
        # check to make sure the supplied specimens actually exist within the data
        warnings = []
        for specimen_id in specimen_ids:
            if specimen_id > len(pmodata["specimen_info"]):
                warnings.append(f"{specimen_id} id is beyond the length of specimen_info: " + str(len(pmodata["specimen_info"])) )
        if len(warnings) > 0:
            raise Exception("\n".join(warnings))
        for experiment_sample_id, experiment_sample in enumerate(pmodata["experiment_info"]):
            if experiment_sample["specimen_id"] in specimen_ids:
                ret[experiment_sample["specimen_id"]].add(experiment_sample_id)
        return ret


    @staticmethod
    def count_targets_per_sample(pmodata, min_reads: float = 0.0) -> pd.DataFrame:
        """
        Count the number of targets per experimental sample, with optional collapsing across bioinformatics runs.

        :param pmodata: the loaded PMO
        :param min_reads: a minimum number of reads for a target in order for it to be counted
        :return: a pandas DataFrame, columns = [bioinformatics_run_id, experiment_sample_name, target_number]
        """
        records = []
        experiment_info = pmodata["experiment_info"]

        for result in pmodata["microhaplotypes_detected"]:
            run_id = result["bioinformatics_run_id"]
            for sample in result["experiment_samples"]:
                sample_id = sample["experiment_sample_id"]
                sample_name = experiment_info[sample_id]["experiment_sample_name"]

                target_count = sum(
                    sum(hap["reads"] for hap in target["haps"]) >= min_reads
                    for target in sample["target_results"]
                )

                record = {
                    "bioinformatics_run_id": run_id,
                    "experiment_sample_name": sample_name,
                    "target_number": target_count
                }

                records.append(record)
        return pd.DataFrame.from_records(records)

    @staticmethod
    def count_samples_per_target(pmodata, min_reads: float = 0.0,
                                 collapse_across_runs: bool = False) -> pd.DataFrame:
        """
        Count the number of experimental samples per target, optionally collapsing across bioinformatics runs.

        :param pmodata: the loaded PMO
        :param min_reads: the minimum number of reads for a target in order for it to be counted
        :param collapse_across_runs: if True, sums across bioinformatics_run_id per target
        :return: a pandas dataframe
                 - if collapse_across_runs=False: columns = [bioinformatics_run_id, target_name, sample_count]
                 - if collapse_across_runs=True:  columns = [target_name, sample_count]
        """
        records = []
        microhap_targets = pmodata["microhaplotypes_info"]["targets"]
        target_info = pmodata["target_info"]

        for result in pmodata["microhaplotypes_detected"]:
            run_id = result["bioinformatics_run_id"]
            target_sample_counts = defaultdict(int)

            for sample in result["experiment_samples"]:
                for target_result in sample["target_results"]:
                    if sum(hap["reads"] for hap in target_result["haps"]) >= min_reads:
                        mhaps_target_id = target_result["mhaps_target_id"]
                        target_id = microhap_targets[mhaps_target_id]["target_id"]
                        target_name = target_info[target_id]["target_name"]
                        target_sample_counts[target_name] += 1

            for target_name, count in target_sample_counts.items():
                records.append({
                    "bioinformatics_run_id": run_id,
                    "target_name": target_name,
                    "sample_count": count
                })

        ret = pd.DataFrame.from_records(records)

        if collapse_across_runs:
            ret = ret.groupby("target_name", as_index=False)["sample_count"].sum()
            ret = ret[["target_name", "sample_count"]]
            return ret.sort_values(by="target_name").reset_index(drop=True)
        return ret.sort_values(by=["bioinformatics_run_id", "target_name"]).reset_index(drop=True)

    @staticmethod
    def list_experiment_sample_ids_per_specimen_id(pmodata, select_specimen_ids: list[str] = None) -> pandas.DataFrame:
        """
        List the experiment_sample_id per specimen_id
        :param pmodata: the PMO
        :param select_specimen_ids: a list of specimen_ids to select, if None, all specimen_ids are used
        :return: a pandas dataframe with 3 columns, specimen_id, experiment_sample_id, and experiment_sample_id_count(the number of experiment_sample_ids per specimen_id)
        """
        exp_samples_per_spec = defaultdict(list[str])
        for exp_sample in pmodata["experiment_info"]:
            if select_specimen_ids is None or exp_sample["specimen_id"] in select_specimen_ids:
                exp_samples_per_spec[pmodata["specimen_info"][exp_sample["specimen_id"]]["specimen_name"]].append(exp_sample["experiment_sample_name"])

        specimens_not_list = []
        for specimen in pmodata["specimen_info"]:
            if specimen["specimen_name"] not in exp_samples_per_spec:
                specimens_not_list.append(specimen["specimen_name"])

        # Prepare the data for DataFrame creation
        data = []
        for specimen_name, experiment_sample_names in exp_samples_per_spec.items():
            for experiment_sample_name in experiment_sample_names:
                data.append({
                    "specimen_name": specimen_name,
                    "experiment_sample_name": experiment_sample_name,
                    "experiment_sample_count": len(experiment_sample_names)
                })

        # Create the DataFrame
        df = pd.DataFrame(data, columns=["specimen_name", "experiment_sample_name", "experiment_sample_count"])
        return df

    @staticmethod
    def count_specimen_per_meta_fields(pmodata) -> pd.DataFrame:
        """
        Get a pandas dataframe of counts of the meta fields within the specimen_info section
        :param pmodata: the pmo to count from
        :return: a pandas dataframe of counts with the following columns: field, present_in_specimens_count, total_specimen_count
        """
        field_counts = defaultdict(int)
        for specimen in pmodata["specimen_info"]:
            for meta_field in specimen:
                field_counts[meta_field] += 1
        counts_df = pd.DataFrame(columns=["field", "present_in_specimens_count", "total_specimen_count"])
        for field_name, field_count in field_counts.items():
            counts_df.loc[len(counts_df)] = {"field": field_name, "present_in_specimens_count": field_count,
                                             "total_specimen_count": len(pmodata["specimen_info"])}
        return counts_df

    @staticmethod
    def count_specimen_by_field_value(pmodata, meta_fields: list[str]) -> pd.DataFrame:
        """
        Count the values of the meta fields. If a specimen doesn't have a field, it is marked as 'NA'.
        Groups are combinations of all given meta fields.

        :param pmodata: the pmo to count from
        :param meta_fields: the fields to get counts for
        :return: counts for all sub-field groups, with metadata
        """
        total_specimens = len(pmodata["specimen_info"])
        field_counts = defaultdict(int)

        for specimen in pmodata["specimen_info"]:
            key = tuple(str(specimen.get(field, "NA")) for field in meta_fields)
            field_counts[key] += 1

        records = []
        for key, count in field_counts.items():
            record = dict(zip(meta_fields, key))
            record.update({
                "specimens_count": count,
                "specimens_freq": count / total_specimens,
                "total_specimen_count": total_specimens
            })
            records.append(record)

        return pd.DataFrame.from_records(records).sort_values(by=meta_fields).reset_index(drop=True)

    @staticmethod
    def extract_allele_counts_freq_from_pmo(pmodata,
                                            select_bioids: list[int] = None,
                                            select_exp_samples: list[str] = None,
                                            select_targets: list[str] = None,
                                            collapse_across_runs: bool = False) -> pd.DataFrame:
        """
        Extract allele counts from PMO data into a single DataFrame.

        :param pmodata: the pmo data structure
        :param select_bioids: optional list of bioinformatics_run_ids to include
        :param select_exp_samples: optional list of experiment_sample_names to include
        :param select_targets: optional list of target_names to include
        :param collapse_across_runs: whether to collapse count/freqs across bioinformatics_run_id runs
        :return: DataFrame with columns: bioinformatics_run_id, target, mhap_id, count, freq, target_total
        """
        with open(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                               "etc/portable_microhaplotype_object.schema.json")) as f:
            checker = PMOChecker(json.load(f))
            checker.check_for_required_base_fields(pmodata)

        allele_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        target_totals = defaultdict(lambda: defaultdict(int))

        for data_for_run in pmodata["microhaplotypes_detected"]:
            bioid = data_for_run["bioinformatics_run_id"]
            if select_bioids is not None and bioid not in select_bioids:
                continue
            for sample_data in data_for_run["experiment_samples"]:
                sample_name = pmodata["experiment_info"][sample_data["experiment_sample_id"]]["experiment_sample_name"]
                if select_exp_samples is not None and sample_name not in select_exp_samples:
                    continue
                for target_data in sample_data["target_results"]:
                    target_id = pmodata["microhaplotypes_info"]["targets"][target_data["mhaps_target_id"]]["target_id"]
                    target = pmodata["target_info"][target_id]["target_name"]
                    if select_targets is not None and target not in select_targets:
                        continue
                    for microhapid in target_data["haps"]:
                        mhap_id = microhapid["mhap_id"]
                        allele_counts[bioid][target][mhap_id] += 1
                        target_totals[bioid][target] += 1

        # Flatten into list of rows
        rows = []
        for bioid, targets in allele_counts.items():
            for target, mhap_counts in targets.items():
                total = target_totals[bioid][target]
                for mhap_id, count in mhap_counts.items():
                    freq = count / total if total > 0 else 0.0
                    rows.append({
                        "bioinformatics_run_id": bioid,
                        "target": target,
                        "mhap_id": mhap_id,
                        "count": count,
                        "freq": freq,
                        "target_total": total
                    })
        ret = pd.DataFrame(rows)

        if collapse_across_runs:
            # Aggregate counts across runs
            collapsed = (
                ret.groupby(["target", "mhap_id"], as_index=False)["count"]
                .sum()
            )
            # Recalculate target_total as sum of counts per target
            total_counts = (
                collapsed.groupby("target", as_index=False)["count"]
                .sum()
                .rename(columns={"count": "target_total"})
            )
            collapsed = collapsed.merge(total_counts, on="target", how="left")
            collapsed["freq"] = collapsed["count"] / collapsed["target_total"]

            # Sort output
            return collapsed.sort_values(["target", "mhap_id"]).reset_index(drop=True)[
                ["target", "mhap_id", "count", "freq", "target_total"]
            ]

        return ret.sort_values(["bioinformatics_run_id", "target", "mhap_id"]).reset_index(drop=True)


    @staticmethod
    def extract_alleles_per_sample_table(pmodata,
                                       additional_specimen_info_fields: list[str] = None,
                                       additional_experiment_info_fields: list[str] = None,
                                       additional_microhap_fields: list[str] = None,
                                       additional_representative_info_fields: list[str] = None,
                                       default_base_col_names: list[str] = ["sampleID", "locus", "allele"]) -> pd.DataFrame:
        """
        Create a pd.Dataframe of sample, target and allele. Can optionally add on any other additional fields

        :param output_delimiter: the delimiter used to write the output file
        :param pmodata: the data to write from
        :param additional_specimen_info_fields: any additional fields to write from the specimen_info object
        :param additional_experiment_info_fields: any additional fields to write from the experiment_samples object
        :param additional_microhap_fields: any additional fields to write from the microhap object
        :param additional_representative_info_fields: any additional fields to write from the representative_microhaplotype_sequences object
        :param default_base_col_names: The default column name for the sample, locus and allele
        :return: pandas dataframe
        """

        # check input
        with open(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                               "etc/portable_microhaplotype_object.schema.json")) as f:
            checker = PMOChecker(json.load(f))
            checker.check_for_required_base_fields(pmodata)

        # Check to see if at least 1 sample has supplied meta field
        # samples without this meta field will have NA
        if additional_specimen_info_fields is not None:
            # Find meta fields that have at least some data
            meta_fields_with_data = {
                metafield
                for metafield in additional_specimen_info_fields
                for specimen_data in pmodata["specimen_info"]
                if metafield in specimen_data
            }

            # Determine meta fields with no samples having data
            meta_fields_with_no_samples = set(additional_specimen_info_fields) - meta_fields_with_data

            if meta_fields_with_no_samples:
                raise Exception(f"No specimen_info have data for fields: {', '.join(meta_fields_with_no_samples)}")
        # Check to see if at least 1 sample has supplied meta field
        # samples without this meta field will have NA
        if additional_experiment_info_fields is not None:
            # Find meta fields that have at least some data
            meta_fields_with_data = {
                metafield
                for metafield in additional_experiment_info_fields
                for experiment_data in pmodata["experiment_info"]
                if metafield in experiment_data
            }
            # Determine meta fields with no samples having data
            meta_fields_with_no_samples = set(additional_experiment_info_fields) - meta_fields_with_data

            if meta_fields_with_no_samples:
                raise Exception(f"No experiment_info have data for fields: {', '.join(meta_fields_with_no_samples)}")

        # Check to see if at least 1 haplotype has this field
        # samples without this meta field will have NA
        if additional_microhap_fields is not None:
            # Find meta fields that have at least some data
            additional_microhap_fields_with_data = {
                additional_microhap_field
                for additional_microhap_field in additional_microhap_fields
                for microhaplotypes_detected in pmodata["microhaplotypes_detected"]
                for experiment_samples_data in microhaplotypes_detected["experiment_samples"]
                for target_data in experiment_samples_data["target_results"]
                for microhap_data in target_data["haps"]
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
        if additional_representative_info_fields is not None:
            # Find meta fields that have at least some data
            additional_microhap_fields_with_data = {
                additional_microhap_field
                for additional_microhap_field in additional_representative_info_fields
                for target_data in pmodata["microhaplotypes_info"]["targets"]
                for microhap_data in target_data["microhaplotypes"]
                if additional_microhap_field in microhap_data
            }
            # Determine meta fields with no samples having data
            additional_microhap_fields_with_no_samples = set(
                additional_representative_info_fields) - additional_microhap_fields_with_data

            if additional_microhap_fields_with_no_samples:
                raise Exception(
                    f"No representative_microhaplotype_sequences have data for fields: {', '.join(additional_microhap_fields_with_no_samples)}")

        if len(default_base_col_names) != 3:
            raise Exception("Must have 3 default columns for allele counts, not {}".format(len(default_base_col_names)))

        rows = []
        specimen_info = pmodata["specimen_info"]
        target_info = pmodata["target_info"]
        experiment_info = pmodata["experiment_info"]
        detected_microhaps = pmodata["microhaplotypes_detected"]
        rep_haps = pmodata["microhaplotypes_info"]["targets"]
        for bio_run_for_detected_microhaps in detected_microhaps:
            bioinformatics_run_id = bio_run_for_detected_microhaps["bioinformatics_run_id"]
            for sample_data in bio_run_for_detected_microhaps["experiment_samples"]:
                experiment_sample_id = sample_data["experiment_sample_id"]
                specimen_id = experiment_info[experiment_sample_id]["specimen_id"]
                experimental_meta = experiment_info[experiment_sample_id]
                specimen_meta = specimen_info[specimen_id]
                for target_data in sample_data["target_results"]:
                    target_name = target_info[rep_haps[target_data["mhaps_target_id"]]["target_id"]]["target_name"]
                    for microhap_data in target_data["haps"]:
                        allele_id = microhap_data["mhap_id"]
                        #print(rep_haps[target_data["mhaps_target_id"]])
                        rep_hap_meta = rep_haps[target_data["mhaps_target_id"]]["microhaplotypes"][allele_id]
                        row = {
                            "bioinformatics_run_id"  : bioinformatics_run_id,
                            default_base_col_names[0]: specimen_meta["specimen_name"],
                            default_base_col_names[1]: target_name,
                            default_base_col_names[2]: allele_id
                        }
                        if additional_experiment_info_fields is not None:
                            for field in additional_experiment_info_fields:
                                row[field] =experimental_meta.get(field, "NA")
                        if additional_specimen_info_fields is not None:
                            for field in additional_specimen_info_fields:
                                row[field] = specimen_meta.get(field, "NA")
                        if additional_microhap_fields is not None:
                            for field in additional_microhap_fields:
                                row[field] = microhap_data.get(field, "NA")
                        if additional_representative_info_fields is not None:
                            for field in additional_representative_info_fields:
                                row[field] = rep_hap_meta.get(field, "NA")
                        rows.append(row)
        # Build and return DataFrame
        return pd.DataFrame(rows)


    @staticmethod
    def filter_pmo_by_experiment_sample_ids(pmodata, experiment_sample_ids: set[int]):
        """
        Extract out of a load PMO the data associated with select experiment_sample_ids
        :param pmodata:the loaded PMO
        :param experiment_sample_ids: the experiment_sample_ids to extract the info for
        :return: a new PMO with only the data associated with the supplied experiment_sample_ids
        """

        # create a new pmo out
        # pmo_name, panel_info, sequencing_info, taramp_bioinformatics_info will stay the same
        # specimen_info, experiment_info, microhaplotypes_detected, representative_microhaplotype_sequences will be
        # created based on the supplied experiment ids

        # check to make sure the supplied specimens actually exist within the data
        warnings = []
        for experiment_sample_id in experiment_sample_ids:
            if experiment_sample_id > len(pmodata["experiment_info"]):
                warnings.append(f"{experiment_sample_id} id is beyond the length of experiment_info: " + str(len(pmodata["experiment_info"])) )
        if len(warnings) > 0:
            raise Exception("\n".join(warnings))

        pmo_out = {"pmo_header": pmodata["pmo_header"],
                   "panel_info": pmodata["panel_info"],
                   "sequencing_info": pmodata["sequencing_info"],
                   "target_info": pmodata["target_info"],
                   "targeted_genomes": pmodata["targeted_genomes"],
                   "microhaplotypes_info": pmodata["microhaplotypes_info"],
                   "bioinformatics_methods_info": pmodata["bioinformatics_methods_info"],
                   "bioinformatics_run_info": pmodata["bioinformatics_run_info"],
                   "specimen_info" : [],
                   "experiment_info" : [],
                   "microhaplotypes_detected" : []
                   }
        if "read_counts_by_stage" in pmodata:
            pmo_out["read_counts_by_stage"] = []
        # need to update read_counts_by_stage, experiment_info, specimen_info, microhaplotypes_detected

        # specimen_info
        # first get the specimen_ids needed and then build
        specimen_ids = set()
        specimen_id_index_key = {}
        for experiment_sample_id in experiment_sample_ids:
            specimen_ids.add(pmodata["experiment_info"][experiment_sample_id]["specimen_id"])
        for specimen_id in specimen_ids:
            specimen_id_index_key[specimen_id] = len(pmo_out["specimen_info"])
            pmo_out["specimen_info"].append(pmodata["specimen_info"][specimen_id])

        # experiment_info
        experiment_id_index_key = {}
        for experiment_sample_id in experiment_sample_ids:
            experiment_id_index_key[experiment_sample_id] = len(pmo_out["experiment_info"])
            pmo_out["experiment_info"].append(pmodata["experiment_info"][experiment_sample_id])
            # update specimen_id
            pmo_out["experiment_info"][len(pmo_out["experiment_info"])-1]["specimen_id"] = specimen_id_index_key[pmodata["experiment_info"][experiment_sample_id]["specimen_id"]]

        # microhaplotypes_detected
        for microhaplotypes_detected in pmodata["microhaplotypes_detected"]:
            new_microhaplotypes_detected = {"bioinformatics_run_id": microhaplotypes_detected["bioinformatics_run_id"],
                                            "experiment_samples": []}
            for sample in microhaplotypes_detected["experiment_samples"]:
                if sample["experiment_sample_id"] in experiment_sample_ids:
                    new_microhaplotypes_detected["experiment_samples"].append(sample)
                    # update experiment_sample_id
                    new_microhaplotypes_detected["experiment_samples"][len(new_microhaplotypes_detected["experiment_samples"]) - 1]["experiment_sample_id"] = experiment_id_index_key[sample["experiment_sample_id"]]
            pmo_out["microhaplotypes_detected"].append(new_microhaplotypes_detected)
        # read_counts_by_stage
        if "read_counts_by_stage" in pmodata:
            for read_count in pmodata["read_counts_by_stage"]:
                new_read_count = {"bioinformatics_run_id": read_count["bioinformatics_run_id"],
                                  "read_counts_by_experimental_sample_by_stage": []}
                for sample in read_count["read_counts_by_experimental_sample_by_stage"]:
                    if sample["experiment_sample_id"] in experiment_sample_ids:
                        new_read_count["read_counts_by_experimental_sample_by_stage"].append(sample)
                        # update experiment_sample_id
                        new_read_count["read_counts_by_experimental_sample_by_stage"][len(new_read_count["read_counts_by_experimental_sample_by_stage"]) - 1]["experiment_sample_id"] = experiment_id_index_key[sample["experiment_sample_id"]]
        return pmo_out
    @staticmethod
    def filter_pmo_by_experiment_sample_names(pmodata, experiment_sample_names: set[str]):
        """
        Filters pmodata by experiment sample names
        :param pmodata: the pmodata object
        :param experiment_sample_names: set of experiment sample names, will be converted into indexes to extract out
        :return: filtered pmodata object containing only the indexes
        """
        experiment_sample_names_list = sorted(list(experiment_sample_names))
        experiment_sample_ids_list = PMOProcessor.get_index_of_experiment_sample_names(pmodata,experiment_sample_names_list)
        return PMOProcessor.filter_pmo_by_experiment_sample_ids(pmodata, set(experiment_sample_ids_list))

    @staticmethod
    def filter_pmo_by_specimen_ids(pmodata, specimen_ids: set[int]):
        """
        Extract out of a load PMO the data associated with select specimen_ids
        :param pmodata:the loaded PMO
        :param specimen_ids: the specimen_ids to extract the info for
        :return: a new PMO with only the data associated with the supplied specimen_ids
        """
        # check to make sure the supplied specimens actually exist within the data
        warnings = []
        for specimen_id in specimen_ids:
            if specimen_id > len(pmodata["specimen_info"]):
                warnings.append(f"{specimen_id} id is beyond the length of specimen_info: " + str(len(pmodata["specimen_info"])) )
        if len(warnings) > 0:
            raise Exception("\n".join(warnings))
        experiment_sample_ids_for_specimen_ids = PMOProcessor.get_experiment_ids_for_specimen_ids(pmodata, specimen_ids)
        all_experiment_sample_ids = {exp_samp for spec in experiment_sample_ids_for_specimen_ids.values() for exp_samp in spec}
        return PMOProcessor.filter_pmo_by_experiment_sample_ids(pmodata, all_experiment_sample_ids)
    @staticmethod
    def filter_pmo_by_specimen_names(pmodata, specimen_names: set[str]):
        """
        Extract out of a load PMO the data associated with select specimen_ids
        :param pmodata:the loaded PMO
        :param specimen_names: the specimen_names to extract the info for
        :return: a new PMO with only the data associated with the supplied specimen_names
        """
        specimen_names_list = sorted(list(specimen_names))
        specimen_ids_list = PMOProcessor.get_index_of_specimen_names(pmodata,specimen_names_list)
        return PMOProcessor.filter_pmo_by_specimen_ids(pmodata, set(specimen_ids_list))

    @staticmethod
    def filter_pmo_by_target_ids(pmodata, target_ids: set[int]):
        """
        Extract out data from the PMO for only select target IDs
        :param pmodata: the pmo to extract data from
        :param target_ids: the target_ids to extract
        :return: a new pmo with the data for only the targets supplied
        """
        # create a new pmo out

        # check to make sure the supplied specimens actually exist within the data
        warnings = []
        for target_id in target_ids:
            if target_id >= len(pmodata["target_info"]):
                warnings.append(f"{target_id} out of range of target_info, length is {len(pmodata['target_info'])}")
        target_ids_in_microhaplotypes_info = []
        for target in pmodata["microhaplotypes_info"]["targets"]:
            target_ids_in_microhaplotypes_info.append(target["target_id"])
        for target_id in target_ids:
            if target_id not in target_ids_in_microhaplotypes_info:
                warnings.append(f"{target_id} not in pmodata[\"microhaplotypes_info\"]")

        if len(warnings) > 0:
            raise Exception("\n".join(warnings))

        pmo_out = {"pmo_header": pmodata["pmo_header"],
                   "sequencing_info": pmodata["sequencing_info"],
                   "specimen_info": pmodata["specimen_info"],
                   "experiment_info": pmodata["experiment_info"],
                   "bioinformatics_methods_info": pmodata["bioinformatics_methods_info"],
                   "bioinformatics_run_info": pmodata["bioinformatics_run_info"],
                   "targeted_genomes": pmodata["targeted_genomes"], "target_info": []}
        # will need to update target_info, panel_info, microhaplotypes_info, microhaplotypes_detected, read_counts_by_stage based
        # on target_ids selecting for
        # first update microhaplotypes_info, microhaplotypes_detected, read_counts_by_stage
        # then update target_info, panel_info
        # then update the target_ids

        # target_info
        pmo_out["target_info"] = []
        target_info_index_key = {}
        for target_info_id, target_info in enumerate(pmodata["target_info"]):
            if target_info_id in target_ids:
                target_info_index_key[target_info_id] = len(pmo_out["target_info"])
                pmo_out["target_info"].append(target_info)

        # panel_info
        pmo_out["panel_info"] = []
        for panel_info in pmodata["panel_info"]:
            new_panel_info = {"panel_name": panel_info["panel_name"],"reactions": []}
            for reaction in panel_info["reactions"]:
                new_reaction = {"reaction_name": reaction["reaction_name"], "panel_targets" : []}
                for panel_target_id in reaction["panel_targets"]:
                    if panel_target_id in target_ids:
                        # add new updated target_id index
                        new_reaction["panel_targets"].append(target_info_index_key[panel_target_id])
                if len(new_reaction["panel_targets"]) > 0:
                    new_panel_info["reactions"].append(new_reaction)
        # microhaplotypes_info
        pmo_out["microhaplotypes_info"] = {"targets" :[]}
        # key=old_mhaps_target_id, value = new_mhaps_target_id
        mhaps_target_id_new_key = {}
        for microhap_info_index, microhap_info in enumerate(pmodata["microhaplotypes_info"]["targets"]):
            if microhap_info["target_id"] in target_ids:
                mhaps_target_id_new_key[microhap_info_index] = len(pmo_out["microhaplotypes_info"]["targets"])
                # update new target_id index
                microhap_info["target_id"] = target_info_index_key[microhap_info["target_id"]]
                pmo_out["microhaplotypes_info"]["targets"].append(microhap_info)
        # microhaplotypes_info
        pmo_out["microhaplotypes_detected"] = []
        for microhaplotypes_detected in pmodata["microhaplotypes_detected"]:
            new_microhaplotypes_detected = {"bioinformatics_run_id": microhaplotypes_detected["bioinformatics_run_id"],
                                            "experiment_samples" : []}
            for sample in microhaplotypes_detected["experiment_samples"]:
                new_sample = {"experiment_sample_id": sample["experiment_sample_id"],
                              "target_results" : []}
                for target in sample["target_results"]:
                    if target["mhaps_target_id"] in mhaps_target_id_new_key:
                        # update with new mhaps_target_id id
                        target["mhaps_target_id"] = mhaps_target_id_new_key[target["mhaps_target_id"]]
                        new_sample["target_results"].append(target)
                new_microhaplotypes_detected["experiment_samples"].append(new_sample)
            pmo_out["microhaplotypes_detected"].append(new_microhaplotypes_detected)

        # read_counts_by_stage
        if "read_counts_by_stage" in pmodata:
            pmo_out["read_counts_by_stage"] = []
            for read_counts_by_bioid in pmodata["read_counts_by_stage"]:
                new_read_counts_by_bioid = {"bioinformatics_run_id": read_counts_by_bioid["bioinformatics_run_id"], "read_counts_by_experimental_sample_by_stage" : []}
                for sample in read_counts_by_bioid["read_counts_by_experimental_sample_by_stage"]:
                    new_samples = {"experiment_sample_id": sample["experiment_sample_id"],
                                   "total_raw_count" : sample["total_raw_count"]}
                    if "read_counts_for_targets" in sample:
                        new_samples["read_counts_for_targets"] = []
                        for target in sample["read_counts_for_targets"]:
                            if target["target_id"] in target_ids:
                                # update with new target_id index
                                target["target_id"] = target_info_index_key[target["target_id"]]
                                new_samples["read_counts_for_targets"].append(target)
                    new_read_counts_by_bioid["read_counts_by_experimental_sample_by_stage"].append(new_samples)
                pmo_out["read_counts_by_stage"].append(new_read_counts_by_bioid)
        return pmo_out
    @staticmethod
    def filter_pmo_by_target_names(pmodata, target_names: set[str]):
        """
        Extract out data from the PMO for only select target IDs
        :param pmodata: the pmo to extract data from
        :param target_names: the target_names to extract
        :return: a new pmo with the data for only the targets supplied
        """
        target_names_list = sorted(list(target_names))
        target_ids_list = PMOProcessor.get_index_of_target_names(pmodata,target_names_list)
        return PMOProcessor.filter_pmo_by_target_ids(pmodata, set(target_ids_list))

    @staticmethod
    def extract_from_pmo_samples_with_meta_groupings(pmodata, meta_fields_values: str):
        """
        Extract out of a PMO the data associated with specimens that belong to specific meta data groupings
        :param pmodata: the PMO to extract from
        :param meta_fields_values: Meta Fields to include, should either be a table with columns field, values (comma separated values) (and optionally group) or supplied command line as field1=value1,value2,value3:field2=value1,value2;field1=value5,value6, where each group is separated by a semicolon
        :return: a pmodata with the input meta
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
        fields_counts = PMOProcessor.count_specimen_per_meta_fields(pmodata)

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
        all_specimen_names = []
        for specimen in pmodata["specimen_info"]:
            for group_name, meta in selected_meta_groups.items():
                passes_criteria = True
                for field, values in meta.items():
                    if not (field in specimen and str(specimen[field]) in values):
                        passes_criteria = False
                        break
                if passes_criteria:
                    group_counts[group_name] += 1
                    specimen_name = specimen["specimen_name"]
                    all_specimen_names.append(specimen_name)
        # Convert selected_meta_groups to a DataFrame
        group_counts_df = pd.DataFrame.from_dict(selected_meta_groups, orient='index')

        # Add the values from count_dict as a new column in df
        group_counts_df['count'] = group_counts_df.index.map(group_counts)

        # Display the resulting DataFrame
        # Collapse lists into comma-separated strings
        group_counts_df = group_counts_df.map(lambda x: ','.join(x) if isinstance(x, list) else x)
        group_counts_df.index.name = "group"

        all_specimen_ids = set(PMOProcessor.get_index_of_specimen_names(pmodata, all_specimen_names))

        pmo_out = PMOProcessor.filter_pmo_by_specimen_ids(pmodata, all_specimen_ids)

        return pmo_out, group_counts_df

    @staticmethod
    def extract_from_pmo_with_read_filter(pmodata, read_filter: float):
        """
        Extract out data from the PMO with inconclusive read filter
        :param pmodata: the pmo to extract data from
        :param read_filter: the read filter to use, inconclusive filter
        :return: a new pmodata with the data only with detected microhaplotypes above this read filter
        """
        # create a new pmo out
        # majority will be the same, just filtering detected microhaplotypes based on read counts
        # @todo consider updating microhaplotypes_info if certain microhaplotypes are no longer detected in any sample with the given filter
        pmo_out = {"pmo_header": pmodata["pmo_header"],
                   "panel_info": pmodata["panel_info"],
                   "sequencing_info": pmodata["sequencing_info"],
                   "target_info": pmodata["target_info"],
                   "specimen_info": pmodata["specimen_info"],
                   "experiment_info": pmodata["experiment_info"],

                   "targeted_genomes": pmodata["targeted_genomes"],
                   "microhaplotypes_info": pmodata["microhaplotypes_info"],
                   "bioinformatics_methods_info": pmodata["bioinformatics_methods_info"],
                   "bioinformatics_run_info": pmodata["bioinformatics_run_info"],
                   "microhaplotypes_detected": []
                   }
        # if has optional read_counts_by_stage then add as well
        # if does contain, @todo consider updating with new counts now that a filter has been applied

        if "read_counts_by_stage" in pmodata:
            pmo_out["read_counts_by_stage"] = pmodata["read_counts_by_stage"]

        # microhaplotypes_detected
        for microhaplotypes_detected in pmodata["microhaplotypes_detected"]:
            extracted_microhaps_for_id = {
                "bioinformatics_run_id": microhaplotypes_detected["bioinformatics_run_id"],
                "experiment_samples": []}
            for experiment in microhaplotypes_detected["experiment_samples"]:
                targets_for_samples = {"experiment_sample_id": experiment["experiment_sample_id"],
                                       "target_results": []}
                for target in experiment["target_results"]:
                    microhaps_for_target = []
                    for microhap in target["haps"]:
                        if microhap["reads"] >= read_filter:
                            microhaps_for_target.append(microhap)
                    if len(microhaps_for_target) > 0:
                        targets_for_samples["target_results"].append(
                            {"mhaps_target_id" : target["mhaps_target_id"],
                            "haps": microhaps_for_target})
                if len(targets_for_samples["target_results"]) > 0:
                    extracted_microhaps_for_id["experiment_samples"].append(targets_for_samples)
            pmo_out["microhaplotypes_detected"].append(extracted_microhaps_for_id)
        return pmo_out

    @staticmethod
    def write_bed_locs(bed_locs : list[bed_loc_tuple], fnp, add_header: bool = False):
        """
        Write out a list of bed_loc_tuple to a file, will auto overwrite it
        :param bed_locs: a list of bed_loc_tuple
        :param fnp: output file path, will be overwritten if it exists
        :param add_header: add header of #chrom,start end,name,score,strand,ref_seq,extra_info, starts with comment so tools will treat it as a comment line
        """
        with open(fnp, "w") as f:
            if add_header:
                f.write("\t".join(["#chrom", "start", "end", "name", "score", "strand", "ref_seq", "extra_info"]))
            for bed_loc in bed_locs:
                f.write("\t".join(
                    [bed_loc.chrom, str(bed_loc.start), str(bed_loc.end), bed_loc.name,
                     str(bed_loc.score), bed_loc.strand,
                     str(bed_loc.ref_seq), bed_loc.extra_info]))
                f.write("\n")

    @staticmethod
    def extract_targets_insert_bed_loc(pmodata, select_target_ids: list[int] = None, sort_output: bool = True):
        """
        Extract out of a PMO the insert location for targets, will add ref seq if loaded into PMO
        :param pmodata: the PMO to extract from
        :param select_target_ids: a list of target ids to select, if None will select all targets
        :param sort_output: whether to sort output by genomic location
        :return: a list of target inserts, with named tuples with fields: chrom, start, end, name, score, strand, extra_info, ref_seq
        """
        # bed_loc = NamedTuple("bed_loc", [("chrom", str), ("start", int), ("end", int), ("name", str), ("score", float), ("strand", str), ("extra_info", str), ("ref_seq", str)])
        bed_loc_out = []
        if select_target_ids is None:
            select_target_ids = list(range(len(pmodata["target_info"])))
        for target_id in select_target_ids:
            tar = pmodata["target_info"][target_id]
            if "insert_location" not in tar:
                raise Exception("no insert_location in pmodata for target id " + str(target_id) + " target_name " + str(tar["target_name"]) + ", cannot extract insert_location")
            genome_info = pmodata["targeted_genomes"][tar["insert_location"]["genome_id"]]
            genome_name_version = genome_info["name"]  + "_" + genome_info["genome_version"]
            extra_info = str("[") + str("genome_name_version=") + genome_name_version + ";]"
            strand = "+" if "strand" not in tar["insert_location"] else tar["insert_location"]["strand"]
            ref_seq =  "" if "ref_seq" not in tar["insert_location"] else tar["insert_location"]["ref_seq"]
            bed_loc_out.append(bed_loc_tuple(tar["insert_location"]["chrom"], tar["insert_location"]["start"], tar["insert_location"]["end"], tar["target_name"], tar["insert_location"]["end"] - tar["insert_location"]["start"], strand, ref_seq, extra_info))
        if sort_output:
            return sorted(bed_loc_out, key=lambda bed: (bed.chrom, bed.start, bed.end))
        return bed_loc_out

    @staticmethod
    def extract_panels_insert_bed_loc(pmodata, select_panel_ids: list[int] = None, sort_output: bool = True):
        """
        Extract out of a PMO the insert location for panels, will add ref seq if loaded into PMO
        :param pmodata: the PMO to extract from
        :param select_panel_ids: a list of panels ids to select, if None will select all panels
        :param sort_output: whether to sort output by genomic location
        :return: a list of target inserts, with named tuples with fields: chrom, start, end, name, score, strand, extra_info, ref_seq
        """
        bed_loc_out = {}
        if select_panel_ids is None:
            select_panel_ids = list(range(len(pmodata["panel_info"])))
        for panel_id in select_panel_ids:
            bed_loc_out_per_panel = []
            for reaction_id in range(len(pmodata["panel_info"][panel_id]["reactions"])):
                for target_id in pmodata["panel_info"][panel_id]["reactions"][reaction_id]["panel_targets"]:
                    tar = pmodata["target_info"][target_id]
                    if "insert_location" not in tar:
                        raise Exception("no insert_location in pmodata for target id " + str(target_id) + " target_name " + str(tar["target_name"]) + ", cannot extract insert_location")
                    genome_info = pmodata["targeted_genomes"][tar["insert_location"]["genome_id"]]
                    genome_name_version = genome_info["name"]  + "_" + genome_info["genome_version"]
                    extra_info = (str("[") +
                                  "genome_name_version=" + genome_name_version + ";" +
                                  "panel=" + pmodata["panel_info"][panel_id]["panel_name"] + ";" +
                                  "reaction=" + pmodata["panel_info"][panel_id]["reactions"][reaction_id]["reaction_name"] + ";" +
                                  "]")
                    strand = "+" if "strand" not in tar["insert_location"] else tar["insert_location"]["strand"]
                    ref_seq =  "" if "ref_seq" not in tar["insert_location"] else tar["insert_location"]["ref_seq"]
                    bed_loc_out_per_panel.append(bed_loc_tuple(tar["insert_location"]["chrom"], tar["insert_location"]["start"], tar["insert_location"]["end"], tar["target_name"], tar["insert_location"]["end"] - tar["insert_location"]["start"], strand, ref_seq, extra_info))
                if sort_output:
                    return sorted(bed_loc_out_per_panel, key=lambda bed: (bed.chrom, bed.start, bed.end))
            bed_loc_out[panel_id] = bed_loc_out_per_panel
        return bed_loc_out

