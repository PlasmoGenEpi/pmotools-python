#!/usr/bin/env python3
import os
import copy
import pandas as pd

from collections import defaultdict


class PMOProcessor:
    """
    A class to extract info out of a loaded PMO object
    """

    @staticmethod
    def get_index_key_of_bioinformatics_run_names(pmodata):
        """
        Get key of bioinformatics_run_name to index in pmodata["bioinformatics_run_info"]
        :param pmodata: the PMO to get indexes from
        :return: a dictionary of indexes keyed by bioinformatics_run_name
        """
        ret = {}
        for idx, bioinformatics_run in enumerate(pmodata["bioinformatics_run_info"]):
            ret[bioinformatics_run["bioinformatics_run_name"]] = idx
        return ret

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
    def get_index_key_of_library_sample_names(pmodata):
        """
        Get key of library_sample_name to index in pmodata["library_sample_info"]
        :param pmodata: the PMO to get indexes from
        :return: a dictionary of indexes keyed by library_sample_name
        """
        ret = {}
        for idx, library in enumerate(pmodata["library_sample_info"]):
            ret[library["library_sample_name"]] = idx
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
    def get_sorted_bioinformatics_run_names(pmodata) -> list[str]:
        """
        Get a name sorted list of bioinformatics_run_names in pmodata["bioinformatics_run_info"]
        :param pmodata: the PMO to get bioinformatics_run_names from
        :return: a list of all bioinformatics_run_names
        """
        return sorted(
            PMOProcessor.get_index_key_of_bioinformatics_run_names(pmodata).keys()
        )

    @staticmethod
    def get_sorted_specimen_names(pmodata) -> list[str]:
        """
        Get a name sorted list of specimen_names in pmodata["specimen_info"]
        :param pmodata: the PMO to get specimen_names from
        :return: a list of all specimen_names
        """
        return sorted(PMOProcessor.get_index_key_of_specimen_names(pmodata).keys())

    @staticmethod
    def get_sorted_library_sample_names(pmodata) -> list[str]:
        """
        Get a name sorted list of library_sample_names in pmodata["library_sample_info"]
        :param pmodata: the PMO to get library_sample_names from
        :return: a list of all library_sample_names
        """
        return sorted(
            PMOProcessor.get_index_key_of_library_sample_names(pmodata).keys()
        )

    @staticmethod
    def get_sorted_target_names(pmodata) -> list[str]:
        """
        Get a name sorted list of target_names in pmodata["target_info"]
        :param pmodata: the PMO to get target_names from
        :return: a list of all target_names
        """
        return sorted(PMOProcessor.get_index_key_of_target_names(pmodata).keys())

    @staticmethod
    def get_sorted_panel_names(pmodata) -> list[str]:
        """
        Get a name sorted list of panel_names in pmodata["panel_info"]
        :param pmodata: the PMO to get panel_names from
        :return: a list of all panel_names
        """
        return sorted(PMOProcessor.get_index_key_of_panel_names(pmodata).keys())

    @staticmethod
    def get_bioinformatics_run_names(pmodata) -> list[str]:
        """
        Get a list of bioinformatics_run_names in pmodata["bioinformatics_run_info"] in order they appear
        :param pmodata: the PMO to get bioinformatics_run_names from
        :return: a list of all bioinformatics_run_names
        """
        ret = []
        for bioinformatics_run in pmodata["bioinformatics_run_info"]:
            ret.append(bioinformatics_run["bioinformatics_run_name"])
        return ret

    @staticmethod
    def get_specimen_names(pmodata) -> list[str]:
        """
        Get a list of specimen_names in pmodata["specimen_info"] in the order they appear
        :param pmodata: the PMO to get specimen_names from
        :return: a list of all specimen_names
        """
        ret = []
        for specimen in pmodata["specimen_info"]:
            ret.append(specimen["specimen_name"])
        return ret

    @staticmethod
    def get_library_sample_names(pmodata) -> list[str]:
        """
        Get a list of library_sample_names in pmodata["library_sample_info"] in the order they appear
        :param pmodata: the PMO to get library_sample_names from
        :return: a list of all library_sample_names
        """
        ret = []
        for library_sample in pmodata["library_sample_info"]:
            ret.append(library_sample["library_sample_name"])
        return ret

    @staticmethod
    def get_target_names(pmodata) -> list[str]:
        """
        Get a list of target_names in pmodata["target_info"] in the order they appear
        :param pmodata: the PMO to get target_names from
        :return: a list of all target_names
        """
        ret = []
        for target in pmodata["target_info"]:
            ret.append(target["target_name"])
        return ret

    @staticmethod
    def get_panel_names(pmodata) -> list[str]:
        """
        Get a list of panel_names in pmodata["panel_info"] in the order they appear
        :param pmodata: the PMO to get panel_names from
        :return: a list of all panel_names
        """
        ret = []
        for panel in pmodata["panel_info"]:
            ret.append(panel["panel_name"])
        return ret

    @staticmethod
    def get_index_key_of_target_in_representative_microhaplotypes(pmodata):
        """
        Get key of target_name to index for the representative microhaplotypes for the target_name in pmodata["representative_microhaplotypes"]
        :param pmodata: the PMO to get indexes from
        :return: a dictionary of indexes keyed by target_name
        """
        ret = {}
        for idx, representative_microhaplotypes_for_target in enumerate(
            pmodata["representative_microhaplotypes"]["targets"]
        ):
            ret[
                pmodata["target_info"][
                    representative_microhaplotypes_for_target["target_id"]
                ]["target_name"]
            ] = idx
        return ret

    @staticmethod
    def get_index_of_bioinformatics_run_names(
        pmodata, bioinformatics_run_names: list[str]
    ):
        """
        Get index of bioinformatics_run_name in pmodata["bioinformatics_run_info"]
        :param pmodata: the PMO to get indexes from
        :param bioinformatics_run_names: a list of bioinformatics_run_names
        :return: the index of bioinformatics_run_names in pmodata["bioinformatics_run_name"] returned in the same order as bioinformatics_run_names
        """
        bioinformatics_run_key = PMOProcessor.get_index_key_of_bioinformatics_run_names(
            pmodata
        )
        return [bioinformatics_run_key[name] for name in bioinformatics_run_names]

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
    def get_index_of_library_sample_names(pmodata, library_sample_names: list[str]):
        """
        Get index of library_sample_name in pmodata["library_sample_info"]
        :param pmodata: the PMO to get indexes from
        :param library_sample_names: a list of library_sample_names
        :return: the index of library_sample_names in pmodata["library_sample_info"] returned in the same order as library_sample_names
        """
        library_sample_key = PMOProcessor.get_index_key_of_library_sample_names(pmodata)
        return [library_sample_key[name] for name in library_sample_names]

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
    def get_index_of_target_in_representative_microhaplotypes(
        pmodata, target_names: list[str]
    ):
        """
        Get index of target_name in pmodata["representative_microhaplotypes"]["targets"]
        :param pmodata: the PMO to get indexes from
        :param target_names: a list of target_names
        :return: the index of target_names in pmodata["representative_microhaplotypes"]["targets"] returned in the same order as target_names
        """
        microhap_target_key = (
            PMOProcessor.get_index_key_of_target_in_representative_microhaplotypes(
                pmodata
            )
        )
        return [microhap_target_key[name] for name in target_names]

    @staticmethod
    def get_library_ids_for_specimen_ids(pmodata, specimen_ids: set[int]):
        """
        get a dictionary that lists the library_ids for a specimen_id
        :param pmodata: the PMO to get indexes from
        :param specimen_ids: a set of specimen_ids
        :return: a dictionary that lists the library_ids for a specimen_id
        """
        ret = defaultdict(set)
        # check to make sure the supplied specimens actually exist within the data
        warnings = []
        for specimen_id in specimen_ids:
            if specimen_id > len(pmodata["specimen_info"]):
                warnings.append(
                    f"{specimen_id} id is beyond the length of specimen_info: "
                    + str(len(pmodata["specimen_info"]))
                )
        if len(warnings) > 0:
            raise Exception("\n".join(warnings))
        for library_sample_id, library_sample in enumerate(
            pmodata["library_sample_info"]
        ):
            if library_sample["specimen_id"] in specimen_ids:
                ret[library_sample["specimen_id"]].add(library_sample_id)
        return ret

    @staticmethod
    def count_targets_per_library_sample(
        pmodata, min_reads: float = 0.0
    ) -> pd.DataFrame:
        """
        Count the number of targets per library sample, with optional collapsing across bioinformatics runs.

        :param pmodata: the loaded PMO
        :param min_reads: a minimum number of reads for a target in order for it to be counted
        :return: a pandas DataFrame, columns = [bioinformatics_run_id, library_sample_name, target_number]
        """
        records = []
        library_sample_info = pmodata["library_sample_info"]

        for result in pmodata["detected_microhaplotypes"]:
            run_id = result["bioinformatics_run_id"]
            for sample in result["library_samples"]:
                sample_id = sample["library_sample_id"]
                sample_name = library_sample_info[sample_id]["library_sample_name"]

                target_count = sum(
                    sum(hap["reads"] for hap in target["mhaps"]) >= min_reads
                    for target in sample["target_results"]
                )

                record = {
                    "bioinformatics_run_id": run_id,
                    "library_sample_name": sample_name,
                    "target_number": target_count,
                }

                records.append(record)
        return pd.DataFrame.from_records(records)

    @staticmethod
    def count_library_samples_per_target(
        pmodata, min_reads: float = 0.0, collapse_across_runs: bool = False
    ) -> pd.DataFrame:
        """
        Count the number of library samples per target, optionally collapsing across bioinformatics runs.

        :param pmodata: the loaded PMO
        :param min_reads: the minimum number of reads for a target in order for it to be counted
        :param collapse_across_runs: if True, sums across bioinformatics_run_id per target
        :return: a pandas dataframe
                 - if collapse_across_runs=False: columns = [bioinformatics_run_id, target_name, sample_count]
                 - if collapse_across_runs=True:  columns = [target_name, sample_count]
        """
        records = []
        microhap_targets = pmodata["representative_microhaplotypes"]["targets"]
        target_info = pmodata["target_info"]

        for result in pmodata["detected_microhaplotypes"]:
            run_id = result["bioinformatics_run_id"]
            target_sample_counts = defaultdict(int)

            for sample in result["library_samples"]:
                for target_result in sample["target_results"]:
                    if sum(hap["reads"] for hap in target_result["mhaps"]) >= min_reads:
                        mhaps_target_id = target_result["mhaps_target_id"]
                        target_id = microhap_targets[mhaps_target_id]["target_id"]
                        target_name = target_info[target_id]["target_name"]
                        target_sample_counts[target_name] += 1

            for target_name, count in target_sample_counts.items():
                records.append(
                    {
                        "bioinformatics_run_id": run_id,
                        "target_name": target_name,
                        "sample_count": count,
                    }
                )

        ret = pd.DataFrame.from_records(records)

        if collapse_across_runs:
            ret = ret.groupby("target_name", as_index=False)["sample_count"].sum()
            ret = ret[["target_name", "sample_count"]]
            return ret.sort_values(by="target_name").reset_index(drop=True)
        return ret.sort_values(by=["bioinformatics_run_id", "target_name"]).reset_index(
            drop=True
        )

    @staticmethod
    def count_targets_per_panel(pmodata) -> pd.DataFrame:
        """
        Count the targets per panel.

        :param pmodata: the pmo to count from
        :return: counts for each panel
        """
        # how many targets in each panel
        panels = []
        target_count = []
        for panel in pmodata["panel_info"]:
            panel_targets = []
            panels.append(panel["panel_name"])
            for reaction in panel["reactions"]:
                panel_targets.extend(reaction["panel_targets"])
            target_count.append(len(set(panel_targets)))
        return pd.DataFrame(
            data={"panel_name": panels, "panel_target_count": target_count}
        )

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
        counts_df = pd.DataFrame(
            columns=["field", "present_in_specimens_count", "total_specimen_count"]
        )
        for field_name, field_count in field_counts.items():
            counts_df.loc[len(counts_df)] = {
                "field": field_name,
                "present_in_specimens_count": field_count,
                "total_specimen_count": len(pmodata["specimen_info"]),
            }
        return counts_df

    @staticmethod
    def count_specimen_by_field_value(pmodata, meta_fields: list[str]) -> pd.DataFrame:
        """
        Count the values of the meta fields. If a specimen doesn't have a field, it is marked as 'NA'.
        Groups are combinations of all given meta fields.

        :param pmodata: the pmo to count from
        :param meta_fields: a list of meta fields to count
        :type meta_fields: list[str]
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
            record.update(
                {
                    "specimens_count": count,
                    "specimens_freq": count / total_specimens,
                    "total_specimen_count": total_specimens,
                }
            )
            records.append(record)

        return (
            pd.DataFrame.from_records(records)
            .sort_values(by=meta_fields)
            .reset_index(drop=True)
        )

    @staticmethod
    def extract_allele_counts_freq_from_pmo(
        pmodata,
        bioinformatics_run_ids: list[int] = None,
        library_sample_names: list[str] = None,
        target_names: list[str] = None,
        collapse_across_runs: bool = False,
    ) -> pd.DataFrame:
        """
        Extract allele counts from PMO data into a single DataFrame.

        :param pmodata: the pmo data structure
        :param bioinformatics_run_ids: optional list of bioinformatics_run_ids to include
        :param library_sample_names: optional list of library_sample_names to include
        :param target_names: optional list of target_names to include
        :param collapse_across_runs: whether to collapse count/freqs across bioinformatics_run_id runs
        :return: DataFrame with columns: bioinformatics_run_id, target, mhap_id, count, freq, target_total
        """

        allele_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        target_totals = defaultdict(lambda: defaultdict(int))

        for data_for_run in pmodata["detected_microhaplotypes"]:
            bioid = data_for_run["bioinformatics_run_id"]
            if (
                bioinformatics_run_ids is not None
                and bioid not in bioinformatics_run_ids
            ):
                continue
            for sample_data in data_for_run["library_samples"]:
                sample_name = pmodata["library_sample_info"][
                    sample_data["library_sample_id"]
                ]["library_sample_name"]
                if (
                    library_sample_names is not None
                    and sample_name not in library_sample_names
                ):
                    continue
                for target_data in sample_data["target_results"]:
                    target_id = pmodata["representative_microhaplotypes"]["targets"][
                        target_data["mhaps_target_id"]
                    ]["target_id"]
                    target = pmodata["target_info"][target_id]["target_name"]
                    if target_names is not None and target not in target_names:
                        continue
                    for microhapid in target_data["mhaps"]:
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
                    rows.append(
                        {
                            "bioinformatics_run_id": bioid,
                            "target_name": target,
                            "mhap_id": mhap_id,
                            "count": count,
                            "freq": freq,
                            "total_haps_per_target": total,
                        }
                    )
        ret = pd.DataFrame(rows)

        if collapse_across_runs:
            # Aggregate counts across runs
            collapsed = ret.groupby(["target_name", "mhap_id"], as_index=False)[
                "count"
            ].sum()
            # Recalculate target_total as sum of counts per target
            total_counts = (
                collapsed.groupby("target_name", as_index=False)["count"]
                .sum()
                .rename(columns={"count": "target_total"})
            )
            collapsed = collapsed.merge(total_counts, on="target_name", how="left")
            collapsed["freq"] = collapsed["count"] / collapsed["target_total"]

            # Sort output
            return collapsed.sort_values(["target_name", "mhap_id"]).reset_index(
                drop=True
            )[["target_name", "mhap_id", "count", "freq", "target_total"]]

        return ret.sort_values(
            ["bioinformatics_run_id", "target_name", "mhap_id"]
        ).reset_index(drop=True)

    @staticmethod
    def filter_pmo_by_library_sample_ids(pmodata, library_sample_ids: set[int]):
        """
        Extract out of a load PMO the data associated with select library_sample_ids
        :param pmodata:the loaded PMO
        :param library_sample_ids: the library_sample_ids to extract the info for
        :return: a new PMO with only the data associated with the supplied library_sample_ids
        """

        # create a new pmo out
        # pmo_name, panel_info, sequencing_info, taramp_bioinformatics_info will stay the same
        # specimen_info, library_sample_info, detected_microhaplotypes, representative_microhaplotype_sequences will be
        # created based on the supplied library ids

        # check to make sure the supplied specimens actually exist within the data
        warnings = []
        for library_sample_id in library_sample_ids:
            if library_sample_id > len(pmodata["library_sample_info"]):
                warnings.append(
                    f"{library_sample_id} id is beyond the length of library_sample_info: "
                    + str(len(pmodata["library_sample_info"]))
                )
        if len(warnings) > 0:
            raise Exception("\n".join(warnings))

        pmo_out = {
            "pmo_header": copy.deepcopy(pmodata["pmo_header"]),
            "panel_info": copy.deepcopy(pmodata["panel_info"]),
            "sequencing_info": copy.deepcopy(pmodata["sequencing_info"]),
            "target_info": copy.deepcopy(pmodata["target_info"]),
            "targeted_genomes": copy.deepcopy(pmodata["targeted_genomes"]),
            "representative_microhaplotypes": copy.deepcopy(
                pmodata["representative_microhaplotypes"]
            ),
            "bioinformatics_methods_info": copy.deepcopy(
                pmodata["bioinformatics_methods_info"]
            ),
            "bioinformatics_run_info": copy.deepcopy(
                pmodata["bioinformatics_run_info"]
            ),
            "specimen_info": [],
            "library_sample_info": [],
            "project_info": copy.deepcopy(pmodata["project_info"]),
            "detected_microhaplotypes": [],
        }
        if "read_counts_by_stage" in pmodata:
            pmo_out["read_counts_by_stage"] = []
        # need to update read_counts_by_stage, library_sample_info, specimen_info, detected_microhaplotypes

        # specimen_info
        # first get the specimen_ids needed and then build
        specimen_ids = set()
        specimen_id_index_key = {}
        for library_sample_id in library_sample_ids:
            specimen_ids.add(
                pmodata["library_sample_info"][library_sample_id]["specimen_id"]
            )
        for specimen_id in specimen_ids:
            specimen_id_index_key[specimen_id] = len(pmo_out["specimen_info"])
            pmo_out["specimen_info"].append(
                copy.deepcopy(pmodata["specimen_info"][specimen_id])
            )

        # library_sample_info
        library_id_index_key = {}
        for library_sample_id in library_sample_ids:
            library_id_index_key[library_sample_id] = len(
                pmo_out["library_sample_info"]
            )
            pmo_out["library_sample_info"].append(
                copy.deepcopy(pmodata["library_sample_info"][library_sample_id])
            )
            # update specimen_id
            pmo_out["library_sample_info"][len(pmo_out["library_sample_info"]) - 1][
                "specimen_id"
            ] = specimen_id_index_key[
                pmodata["library_sample_info"][library_sample_id]["specimen_id"]
            ]

        # detected_microhaplotypes
        for detected_microhaplotypes in pmodata["detected_microhaplotypes"]:
            new_detected_microhaplotypes = {
                "bioinformatics_run_id": detected_microhaplotypes[
                    "bioinformatics_run_id"
                ],
                "library_samples": [],
            }
            for sample in detected_microhaplotypes["library_samples"]:
                if sample["library_sample_id"] in library_sample_ids:
                    new_detected_microhaplotypes["library_samples"].append(
                        copy.deepcopy(sample)
                    )
                    # update library_sample_id
                    new_detected_microhaplotypes["library_samples"][
                        len(new_detected_microhaplotypes["library_samples"]) - 1
                    ]["library_sample_id"] = library_id_index_key[
                        sample["library_sample_id"]
                    ]
            pmo_out["detected_microhaplotypes"].append(new_detected_microhaplotypes)
        # read_counts_by_stage
        if "read_counts_by_stage" in pmodata:
            for read_count in pmodata["read_counts_by_stage"]:
                new_read_count = {
                    "bioinformatics_run_id": read_count["bioinformatics_run_id"],
                    "read_counts_by_library_sample_by_stage": [],
                }
                for sample in read_count["read_counts_by_library_sample_by_stage"]:
                    if sample["library_sample_id"] in library_sample_ids:
                        new_read_count["read_counts_by_library_sample_by_stage"].append(
                            copy.deepcopy(sample)
                        )
                        # update library_sample_id
                        new_read_count["read_counts_by_library_sample_by_stage"][
                            len(
                                new_read_count["read_counts_by_library_sample_by_stage"]
                            )
                            - 1
                        ]["library_sample_id"] = library_id_index_key[
                            sample["library_sample_id"]
                        ]
        return pmo_out

    @staticmethod
    def filter_pmo_by_library_sample_names(pmodata, library_sample_names: set[str]):
        """
        Filters pmodata by library sample names
        :param pmodata: the pmodata object
        :param library_sample_names: set of library sample names, will be converted into indexes to extract out
        :return: filtered pmodata object containing only the indexes
        """
        library_sample_names_list = sorted(list(library_sample_names))
        library_sample_ids_list = PMOProcessor.get_index_of_library_sample_names(
            pmodata, library_sample_names_list
        )
        return PMOProcessor.filter_pmo_by_library_sample_ids(
            pmodata, set(library_sample_ids_list)
        )

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
                warnings.append(
                    f"{specimen_id} id is beyond the length of specimen_info: "
                    + str(len(pmodata["specimen_info"]))
                )
        if len(warnings) > 0:
            raise Exception("\n".join(warnings))
        library_sample_ids_for_specimen_ids = (
            PMOProcessor.get_library_ids_for_specimen_ids(pmodata, specimen_ids)
        )
        all_library_sample_ids = {
            exp_samp
            for spec in library_sample_ids_for_specimen_ids.values()
            for exp_samp in spec
        }
        return PMOProcessor.filter_pmo_by_library_sample_ids(
            pmodata, all_library_sample_ids
        )

    @staticmethod
    def filter_pmo_by_specimen_names(pmodata, specimen_names: set[str]):
        """
        Extract out of a load PMO the data associated with select specimen_ids
        :param pmodata:the loaded PMO
        :param specimen_names: the specimen_names to extract the info for
        :return: a new PMO with only the data associated with the supplied specimen_names
        """
        specimen_names_list = sorted(list(specimen_names))
        specimen_ids_list = PMOProcessor.get_index_of_specimen_names(
            pmodata, specimen_names_list
        )
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
                warnings.append(
                    f"{target_id} out of range of target_info, length is {len(pmodata['target_info'])}"
                )
        target_ids_in_representative_microhaplotypes = []
        for target in pmodata["representative_microhaplotypes"]["targets"]:
            target_ids_in_representative_microhaplotypes.append(target["target_id"])
        for target_id in target_ids:
            if target_id not in target_ids_in_representative_microhaplotypes:
                warnings.append(
                    f'{target_id} not in pmodata["representative_microhaplotypes"]'
                )

        if len(warnings) > 0:
            raise Exception("\n".join(warnings))

        pmo_out = {
            "pmo_header": copy.deepcopy(pmodata["pmo_header"]),
            "sequencing_info": copy.deepcopy(pmodata["sequencing_info"]),
            "specimen_info": copy.deepcopy(pmodata["specimen_info"]),
            "project_info": copy.deepcopy(pmodata["project_info"]),
            "library_sample_info": copy.deepcopy(pmodata["library_sample_info"]),
            "bioinformatics_methods_info": copy.deepcopy(
                pmodata["bioinformatics_methods_info"]
            ),
            "bioinformatics_run_info": copy.deepcopy(
                pmodata["bioinformatics_run_info"]
            ),
            "targeted_genomes": copy.deepcopy(pmodata["targeted_genomes"]),
            "target_info": [],
        }
        # function will update target_info, panel_info, representative_microhaplotypes, detected_microhaplotypes, read_counts_by_stage based
        # on target_ids selecting for first update representative_microhaplotypes, detected_microhaplotypes, read_counts_by_stage
        # then update target_info, panel_info
        # then update the target_ids

        # target_info
        pmo_out["target_info"] = []
        target_info_index_key = {}
        for target_info_id, target_info in enumerate(pmodata["target_info"]):
            if target_info_id in target_ids:
                target_info_index_key[target_info_id] = len(pmo_out["target_info"])
                pmo_out["target_info"].append(copy.deepcopy(target_info))

        # panel_info
        pmo_out["panel_info"] = []
        for panel_info in pmodata["panel_info"]:
            new_panel_info = {"panel_name": panel_info["panel_name"], "reactions": []}
            for reaction in panel_info["reactions"]:
                new_reaction = {
                    "reaction_name": reaction["reaction_name"],
                    "panel_targets": [],
                }
                for panel_target_id in reaction["panel_targets"]:
                    if panel_target_id in target_ids:
                        # add new updated target_id index
                        new_reaction["panel_targets"].append(
                            target_info_index_key[panel_target_id]
                        )
                if len(new_reaction["panel_targets"]) > 0:
                    new_panel_info["reactions"].append(new_reaction)
        # representative_microhaplotypes
        pmo_out["representative_microhaplotypes"] = {"targets": []}
        # key=old_mhaps_target_id, value = new_mhaps_target_id
        mhaps_target_id_new_key = {}
        for microhap_info_index, microhap_info in enumerate(
            pmodata["representative_microhaplotypes"]["targets"]
        ):
            if microhap_info["target_id"] in target_ids:
                mhaps_target_id_new_key[microhap_info_index] = len(
                    pmo_out["representative_microhaplotypes"]["targets"]
                )
                # update new target_id index
                microhap_info["target_id"] = target_info_index_key[
                    microhap_info["target_id"]
                ]
                pmo_out["representative_microhaplotypes"]["targets"].append(
                    copy.deepcopy(microhap_info)
                )
        # representative_microhaplotypes
        pmo_out["detected_microhaplotypes"] = []
        for detected_microhaplotypes in pmodata["detected_microhaplotypes"]:
            new_detected_microhaplotypes = {
                "bioinformatics_run_id": detected_microhaplotypes[
                    "bioinformatics_run_id"
                ],
                "library_samples": [],
            }
            for sample in detected_microhaplotypes["library_samples"]:
                new_sample = {
                    "library_sample_id": sample["library_sample_id"],
                    "target_results": [],
                }
                for target in sample["target_results"]:
                    if target["mhaps_target_id"] in mhaps_target_id_new_key:
                        # update with new mhaps_target_id id
                        target["mhaps_target_id"] = mhaps_target_id_new_key[
                            target["mhaps_target_id"]
                        ]
                        new_sample["target_results"].append(copy.deepcopy(target))
                new_detected_microhaplotypes["library_samples"].append(new_sample)
            pmo_out["detected_microhaplotypes"].append(new_detected_microhaplotypes)

        # read_counts_by_stage
        if "read_counts_by_stage" in pmodata:
            pmo_out["read_counts_by_stage"] = []
            for read_counts_by_bioid in pmodata["read_counts_by_stage"]:
                new_read_counts_by_bioid = {
                    "bioinformatics_run_id": read_counts_by_bioid[
                        "bioinformatics_run_id"
                    ],
                    "read_counts_by_library_sample_by_stage": [],
                }
                for sample in read_counts_by_bioid[
                    "read_counts_by_library_sample_by_stage"
                ]:
                    new_samples = {
                        "library_sample_id": sample["library_sample_id"],
                        "total_raw_count": sample["total_raw_count"],
                    }
                    if "read_counts_for_targets" in sample:
                        new_samples["read_counts_for_targets"] = []
                        for target in sample["read_counts_for_targets"]:
                            if target["target_id"] in target_ids:
                                # update with new target_id index
                                target["target_id"] = target_info_index_key[
                                    target["target_id"]
                                ]
                                new_samples["read_counts_for_targets"].append(
                                    copy.deepcopy(target)
                                )
                    new_read_counts_by_bioid[
                        "read_counts_by_library_sample_by_stage"
                    ].append(new_samples)
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
        target_ids_list = PMOProcessor.get_index_of_target_names(
            pmodata, target_names_list
        )
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
            meta_tab = pd.read_csv(meta_fields_values, sep="\t")
            if "field" not in meta_tab or "values" not in meta_tab:
                raise Exception(
                    meta_fields_values
                    + " doesn't have columns field and values, has "
                    + ",".join(meta_tab.columns)
                )
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
                        raise Exception(
                            "error processing " + group_tok,
                            " should be field and values separated by =",
                        )
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
        group_counts_df = pd.DataFrame.from_dict(selected_meta_groups, orient="index")

        # Add the values from count_dict as a new column in df
        group_counts_df["count"] = group_counts_df.index.map(group_counts)

        # Display the resulting DataFrame
        # Collapse lists into comma-separated strings
        group_counts_df = group_counts_df.map(
            lambda x: ",".join(x) if isinstance(x, list) else x
        )
        group_counts_df.index.name = "group"

        all_specimen_ids = set(
            PMOProcessor.get_index_of_specimen_names(pmodata, all_specimen_names)
        )

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
        # @todo consider updating representative_microhaplotypes if certain microhaplotypes are no longer detected in any sample with the given filter
        pmo_out = {
            "pmo_header": copy.deepcopy(pmodata["pmo_header"]),
            "panel_info": copy.deepcopy(pmodata["panel_info"]),
            "sequencing_info": copy.deepcopy(pmodata["sequencing_info"]),
            "target_info": copy.deepcopy(pmodata["target_info"]),
            "specimen_info": copy.deepcopy(pmodata["specimen_info"]),
            "library_sample_info": copy.deepcopy(pmodata["library_sample_info"]),
            "project_info": copy.deepcopy(pmodata["project_info"]),
            "targeted_genomes": copy.deepcopy(pmodata["targeted_genomes"]),
            "representative_microhaplotypes": copy.deepcopy(
                pmodata["representative_microhaplotypes"]
            ),
            "bioinformatics_methods_info": copy.deepcopy(
                pmodata["bioinformatics_methods_info"]
            ),
            "bioinformatics_run_info": copy.deepcopy(
                pmodata["bioinformatics_run_info"]
            ),
            "detected_microhaplotypes": [],
        }
        # if has optional read_counts_by_stage then add as well
        # if does contain, @todo consider updating with new counts now that a filter has been applied

        if "read_counts_by_stage" in pmodata:
            pmo_out["read_counts_by_stage"] = copy.deepcopy(
                pmodata["read_counts_by_stage"]
            )

        # detected_microhaplotypes
        for detected_microhaplotypes in pmodata["detected_microhaplotypes"]:
            extracted_microhaps_for_id = {
                "bioinformatics_run_id": detected_microhaplotypes[
                    "bioinformatics_run_id"
                ],
                "library_samples": [],
            }
            for library in detected_microhaplotypes["library_samples"]:
                targets_for_samples = {
                    "library_sample_id": library["library_sample_id"],
                    "target_results": [],
                }
                for target in library["target_results"]:
                    microhaps_for_target = []
                    for microhap in target["mhaps"]:
                        if microhap["reads"] >= read_filter:
                            microhaps_for_target.append(copy.deepcopy(microhap))
                    if len(microhaps_for_target) > 0:
                        targets_for_samples["target_results"].append(
                            {
                                "mhaps_target_id": target["mhaps_target_id"],
                                "mhaps": microhaps_for_target,
                            }
                        )
                if len(targets_for_samples["target_results"]) > 0:
                    extracted_microhaps_for_id["library_samples"].append(
                        targets_for_samples
                    )
            pmo_out["detected_microhaplotypes"].append(extracted_microhaps_for_id)
        return pmo_out
