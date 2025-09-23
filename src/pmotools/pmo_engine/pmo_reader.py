#!/usr/bin/env python3
import copy
import datetime
import json
import gzip
import os
import sys
from collections import defaultdict


class PMOReader:
    """
    A class for reading in PMO from files
    """

    @staticmethod
    def read_in_pmo(fnp: str | os.PathLike[str]):
        """
        Read in a PMO file, can either be compressed(.gz) or uncompressed
        :param fnp: the file name path of the PMO file to read in
        :return: a PMO like object
        """
        if "STDIN" == fnp:
            pmo_data = json.load(sys.stdin)
        else:
            if fnp.endswith(".gz"):
                with gzip.open(fnp) as f:
                    pmo_data = json.load(f)
            else:
                with open(fnp) as f:
                    pmo_data = json.load(f)
        return pmo_data

    @staticmethod
    def read_in_pmos(fnps: list[str] | list[os.PathLike[str]]):
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
    def combine_multiple_pmos(pmos: list[dict]):
        """
        Combine multiple PMOs into one pmo
        :param pmos: a list of PMO objects
        :return: a combined PMO
        """
        if len(pmos) <= 1:
            raise Exception(
                "Only supplied "
                + str(len(pmos))
                + " but multiple PMO objects were expected"
            )
        # create new pmo out
        pmo_out = {}
        # create new pmo_header
        # currently losing all info about previous header info,
        # consider coming up with something in standard that might preserve this info if needed
        pmo_out["pmo_header"] = {
            "pmo_version": "v1.0.0",
            "creation_date": datetime.datetime.now().strftime("%Y-%m-%d"),
            "generation_method": {
                "program_name": "pmotools-python.PMOReader.combine_multiple_pmos",
                "program_version": "v1.0.0",
            },
        }

        # combine targeted_genomes fields
        pmo_out["targeted_genomes"] = copy.deepcopy(pmos[0]["targeted_genomes"])
        # key: genome name + _ + genome_version, val: index
        targeted_genomes_out_index_key = {}
        for genome_info_index, genome in enumerate(pmos[0]["targeted_genomes"]):
            targeted_genomes_out_index_key[
                genome["name"] + "_" + genome["genome_version"]
            ] = genome_info_index
        # key1 pmo_index, key2 old_index, val new_index
        targeted_genomes_old_index_key = defaultdict(dict)
        for pmo_index, pmo in enumerate(pmos[1:], start=1):
            for genome_index, genome in enumerate(pmo["targeted_genomes"]):
                genome_id = genome["name"] + "_" + genome["genome_version"]
                if genome_id in targeted_genomes_out_index_key:
                    targeted_genomes_old_index_key[pmo_index][
                        genome_index
                    ] = targeted_genomes_out_index_key[genome_id]
                else:
                    new_index = len(pmo_out["targeted_genomes"])
                    pmo_out["targeted_genomes"].append(copy.deepcopy(genome))
                    targeted_genomes_out_index_key[genome_id] = new_index
                    targeted_genomes_old_index_key[pmo_index][genome_index] = new_index

        # combine target_info fields
        pmo_out["target_info"] = copy.deepcopy(pmos[0]["target_info"])
        # key: target_name, val: index
        target_info_out_index_key = {}
        for target_info_index, target_info in enumerate(pmos[0]["target_info"]):
            target_info_out_index_key[target_info["target_name"]] = target_info_index
        # key1 pmo_index, key2 old_index, val new_index
        target_info_old_index_key = defaultdict(dict)
        for pmo_index, pmo in enumerate(pmos[1:], start=1):
            for target_index, target in enumerate(pmo["target_info"]):
                if target["target_name"] in target_info_out_index_key:
                    target_info_old_index_key[pmo_index][
                        target_index
                    ] = target_info_out_index_key[target["target_name"]]
                else:
                    new_index = len(pmo_out["target_info"])
                    target_copy = copy.deepcopy(target)
                    # update genome_id if adding new target
                    if len(pmo_out["targeted_genomes"]) > 1:
                        if "insert_location" in target_copy:
                            # update genome_id
                            target_copy["insert_location"][
                                "genome_id"
                            ] = targeted_genomes_old_index_key[pmo_index][
                                target_copy["insert_location"]["genome_id"]
                            ]
                        if "location" in target_copy["forward_primer"]:
                            # update genome_id
                            target_copy["forward_primer"]["location"][
                                "genome_id"
                            ] = targeted_genomes_old_index_key[pmo_index][
                                target_copy["forward_primer"]["location"]["genome_id"]
                            ]
                        if "location" in target_copy["reverse_primer"]:
                            # update genome_id
                            target_copy["reverse_primer"]["location"][
                                "genome_id"
                            ] = targeted_genomes_old_index_key[pmo_index][
                                target_copy["reverse_primer"]["location"]["genome_id"]
                            ]
                    pmo_out["target_info"].append(target_copy)
                    target_info_out_index_key[target_copy["target_name"]] = new_index
                    target_info_old_index_key[pmo_index][target_index] = new_index

        # combine panel_info
        # todo, more extensive testing than just panel name, make sure reactions and targets are the same
        pmo_out["panel_info"] = copy.deepcopy(pmos[0]["panel_info"])
        # key: panel_name, val: index
        panel_info_out_index_key = {}
        for panel_info_index, panel_info in enumerate(pmos[0]["panel_info"]):
            panel_info_out_index_key[panel_info["panel_name"]] = panel_info_index
        # key1 pmo_index, key2 old_index, val new_index
        panel_info_old_index_key = defaultdict(dict)
        for pmo_index, pmo in enumerate(pmos[1:], start=1):
            for panel_index, panel in enumerate(pmo["panel_info"]):
                panel_copy = copy.deepcopy(panel)
                if panel_copy["panel_name"] in panel_info_out_index_key:
                    panel_info_old_index_key[pmo_index][
                        panel_index
                    ] = panel_info_out_index_key[panel_copy["panel_name"]]
                else:
                    new_index = len(pmo_out["panel_info"])
                    # update target indexes to make sure reaction points to the right target indexes
                    for reaction in panel_copy["reactions"]:
                        for target_id_reaction in range(
                            len(reaction["target_id_reactions"])
                        ):
                            reaction["target_id_reactions"][
                                target_id_reaction
                            ] = target_info_old_index_key[pmo_index][
                                reaction["target_id_reactions"][target_id_reaction]
                            ]
                    pmo_out["panel_info"].append(panel_copy)
                    panel_info_out_index_key[panel_copy["panel_name"]] = new_index
                    panel_info_old_index_key[pmo_index][panel_index] = new_index

        # combine sequencing_info
        # really shouldn't be possible to have the same sequencing_info in different pmos so
        # just concatenate sequencing infos. Only way this could have happened is if files were split into different
        # pmos and then rejoined but even if we concatenate sequencing_info of the same, they will still properly
        # have the right info per library
        pmo_out["sequencing_info"] = copy.deepcopy(pmos[0]["sequencing_info"])
        # key1 pmo_index, key2 old_index, val new_index
        sequencing_info_old_index_key = defaultdict(dict)
        for pmo_index, pmo in enumerate(pmos[1:], start=1):
            for sequencing_info_index, sequencing_info in enumerate(
                pmo["sequencing_info"]
            ):
                new_index = len(pmo_out["sequencing_info"])
                pmo_out["sequencing_info"].append(copy.deepcopy(sequencing_info))
                sequencing_info_old_index_key[pmo_index][
                    sequencing_info_index
                ] = new_index

        # combine project_info
        # could be possible to be combining PMOs across one project so check if project already exists
        pmo_out["project_info"] = copy.deepcopy(pmos[0]["project_info"])
        project_info_old_index_key = defaultdict(dict)
        for pmo_index, pmo in enumerate(pmos[1:], start=1):
            for project_info_index, project_info in enumerate(pmo["project_info"]):
                # check to see if project already exists
                found_project_info = False
                for current_project_id, current_project_info in enumerate(
                    pmo_out["project_info"]
                ):
                    if (
                        current_project_info["project_name"]
                        == project_info["project_name"]
                    ):
                        if (
                            current_project_info["project_description"]
                            != project_info["project_description"]
                        ):
                            raise Exception(
                                "Project description mismatch for project_name: "
                                + project_info["project_name"]
                            )
                        else:
                            project_info_old_index_key[pmo_index][
                                project_info_index
                            ] = current_project_id
                            found_project_info = True
                if not found_project_info:
                    new_index = len(pmo_out["project_info"])
                    pmo_out["project_info"].append(copy.deepcopy(project_info))
                    project_info_old_index_key[pmo_index][
                        project_info_index
                    ] = new_index

        # combine specimen_info and library_sample_info
        # update project_id
        pmo_out["specimen_info"] = copy.deepcopy(pmos[0]["specimen_info"])
        specimen_names = []
        duplicate_specimen_names = []
        for specimen in pmo_out["specimen_info"]:
            if specimen["specimen_name"] in specimen_names:
                duplicate_specimen_names.append(specimen["specimen_name"])
            specimen_names.append(specimen["specimen_name"])

        # key1 pmo_index, key2 old_index, val new_index
        specimen_info_old_index_key = defaultdict(dict)
        for pmo_index, pmo in enumerate(pmos[1:], start=1):
            for specimen_info_index, specimen_info in enumerate(pmo["specimen_info"]):
                # checkin for duplicates
                if specimen_info["specimen_name"] in specimen_names:
                    duplicate_specimen_names.append(specimen_info["specimen_name"])
                specimen_names.append(specimen_info["specimen_name"])
                new_index = len(pmo_out["specimen_info"])
                # update project_id
                specimen_info_copy = copy.deepcopy(specimen_info)
                specimen_info_copy["project_id"] = project_info_old_index_key[
                    pmo_index
                ][specimen_info_copy["project_id"]]
                pmo_out["specimen_info"].append(specimen_info_copy)
                specimen_info_old_index_key[pmo_index][specimen_info_index] = new_index

        ## library_sample_info
        pmo_out["library_sample_info"] = copy.deepcopy(pmos[0]["library_sample_info"])
        # key1 pmo_index, key2 old_index, val new_index
        library_sample_info_old_index_key = defaultdict(dict)
        duplicate_library_sample_names = []
        library_sample_names = []
        for pmo_index, pmo in enumerate(pmos[1:], start=1):
            for library_sample_info_index, library_sample_info in enumerate(
                pmo["library_sample_info"]
            ):
                # checkin for duplicates
                if library_sample_info["library_sample_name"] in library_sample_names:
                    duplicate_library_sample_names.append(
                        library_sample_info["library_sample_name"]
                    )
                library_sample_names.append(library_sample_info["library_sample_name"])
                # update indexes
                library_sample_info_copy = copy.deepcopy(library_sample_info)
                library_sample_info_copy["specimen_id"] = specimen_info_old_index_key[
                    pmo_index
                ][library_sample_info_copy["specimen_id"]]
                library_sample_info_copy["panel_id"] = panel_info_old_index_key[
                    pmo_index
                ][library_sample_info_copy["panel_id"]]
                library_sample_info_copy[
                    "sequencing_info_id"
                ] = sequencing_info_old_index_key[pmo_index][
                    library_sample_info_copy["sequencing_info_id"]
                ]
                # append to the out library_sample_info_copy after getting new index
                new_index = len(pmo_out["library_sample_info"])
                pmo_out["library_sample_info"].append(library_sample_info_copy)
                library_sample_info_old_index_key[pmo_index][
                    library_sample_info_index
                ] = new_index

        warnings = []
        if len(duplicate_specimen_names) > 0:
            warnings.append(
                "Duplicate specimen names were supplied for the following specimens: "
                + ",".join(duplicate_specimen_names)
            )
        if len(duplicate_library_sample_names) > 0:
            warnings.append(
                "Duplicate library sample names were supplied for the following librarys: "
                + ",".join(duplicate_library_sample_names)
            )
        if len(warnings) > 0:
            raise Exception("\n".join(warnings))

        # update bioinformatics_methods_info
        # the different bioinformatics_methods_info might be the same but there's no easy way to perfectly match up right now
        pmo_out["bioinformatics_methods_info"] = copy.deepcopy(
            pmos[0]["bioinformatics_methods_info"]
        )
        # key1 pmo_index, key2 old_index, val new_index
        bioinformatics_methods_info_old_index_key = defaultdict(dict)
        for pmo_index, pmo in enumerate(pmos[1:], start=1):
            for (
                bioinformatics_methods_info_index,
                bioinformatics_methods_info,
            ) in enumerate(pmo["bioinformatics_methods_info"]):
                new_index = len(pmo_out["bioinformatics_methods_info"])
                pmo_out["bioinformatics_methods_info"].append(
                    copy.deepcopy(bioinformatics_methods_info)
                )
                bioinformatics_methods_info_old_index_key[pmo_index][
                    bioinformatics_methods_info_index
                ] = new_index

        # update bioinformatics_run_info
        pmo_out["bioinformatics_run_info"] = copy.deepcopy(
            pmos[0]["bioinformatics_run_info"]
        )
        # key1 pmo_index, key2 old_index, val new_index
        bioinformatics_run_info_old_index_key = defaultdict(dict)
        for pmo_index, pmo in enumerate(pmos[1:], start=1):
            for bioinformatics_run_info_index, bioinformatics_run_info in enumerate(
                pmo["bioinformatics_run_info"]
            ):
                bioinformatics_run_info_copy = copy.deepcopy(bioinformatics_run_info)
                bioinformatics_run_info_copy[
                    "bioinformatics_methods_id"
                ] = bioinformatics_methods_info_old_index_key[pmo_index][
                    bioinformatics_run_info_index
                ]
                new_index = len(pmo_out["bioinformatics_run_info"])
                pmo_out["bioinformatics_run_info"].append(bioinformatics_run_info_copy)
                bioinformatics_run_info_old_index_key[pmo_index][
                    bioinformatics_run_info_index
                ] = new_index

        # update representative_microhaplotypes
        pmo_out["representative_microhaplotypes"] = pmos[0][
            "representative_microhaplotypes"
        ]
        # key: target_name (not index), val: index in representative_microhaplotypes
        representative_microhaplotypes_out_index_key = {}
        for (
            representative_microhaplotypes_index,
            representative_microhaplotypes,
        ) in enumerate(pmo_out["representative_microhaplotypes"]["targets"]):
            representative_microhaplotypes_out_index_key[
                pmo_out["target_info"][representative_microhaplotypes["target_id"]][
                    "target_name"
                ]
            ] = representative_microhaplotypes_index
        # key1: pmo_index, key2: old_mhaps_target_id, val: new_mhaps_target_id
        representative_microhaplotypes_old_index_key = defaultdict(dict)
        # key1: pmo_index, key2: old_mhaps_target_id, key3: old_mhap_id, val: new_mhap_id
        representative_microhaplotypes_hmap_for_target_index_old_index_key = (
            defaultdict(lambda: defaultdict(dict))
        )
        for pmo_index, pmo in enumerate(pmos[1:], start=1):
            for (
                representative_microhaplotypes_index,
                representative_microhaplotypes,
            ) in enumerate(pmo["representative_microhaplotypes"]["targets"]):
                if (
                    pmo["target_info"][representative_microhaplotypes["target_id"]][
                        "target_name"
                    ]
                    in representative_microhaplotypes_out_index_key
                ):
                    representative_microhaplotypes_old_index_key[pmo_index][
                        representative_microhaplotypes_index
                    ] = representative_microhaplotypes_out_index_key[
                        pmo["target_info"][representative_microhaplotypes["target_id"]][
                            "target_name"
                        ]
                    ]
                    # now update per microhaplotype
                    for adding_microhap_index, adding_microhap in enumerate(
                        representative_microhaplotypes["microhaplotypes"]
                    ):
                        found = False
                        # print(pmo_out["representative_microhaplotypes"]["targets"][representative_microhaplotypes_out_index_key[pmo["target_info"][representative_microhaplotypes["target_id"]]["target_name"]]]["microhaplotypes"])
                        # print("\n\n\n")
                        for (
                            already_have_microhap_index,
                            already_have_microhap,
                        ) in enumerate(
                            pmo_out["representative_microhaplotypes"]["targets"][
                                representative_microhaplotypes_out_index_key[
                                    pmo["target_info"][
                                        representative_microhaplotypes["target_id"]
                                    ]["target_name"]
                                ]
                            ]["microhaplotypes"]
                        ):
                            # print(already_have_microhap)
                            # print("\n")
                            if adding_microhap["seq"] == already_have_microhap["seq"]:
                                representative_microhaplotypes_hmap_for_target_index_old_index_key[
                                    pmo_index
                                ][representative_microhaplotypes_index][
                                    adding_microhap_index
                                ] = already_have_microhap_index
                                found = True
                                break
                        if not found:
                            new_index = len(
                                pmo_out["representative_microhaplotypes"]["targets"][
                                    representative_microhaplotypes_out_index_key[
                                        pmo["target_info"][
                                            representative_microhaplotypes["target_id"]
                                        ]["target_name"]
                                    ]
                                ]["microhaplotypes"]
                            )
                            pmo_out["representative_microhaplotypes"]["targets"][
                                representative_microhaplotypes_out_index_key[
                                    pmo["target_info"][
                                        representative_microhaplotypes["target_id"]
                                    ]["target_name"]
                                ]
                            ]["microhaplotypes"].append(copy.deepcopy(adding_microhap))
                            representative_microhaplotypes_hmap_for_target_index_old_index_key[
                                pmo_index
                            ][representative_microhaplotypes_index][
                                adding_microhap_index
                            ] = new_index
                else:
                    # if not currently in representative_microhaplotypes, update keys and look-ups
                    new_mhaps_target_index = len(
                        pmo_out["representative_microhaplotypes"]["targets"]
                    )
                    pmo_out["representative_microhaplotypes"]["targets"].append(
                        copy.deepcopy(representative_microhaplotypes)
                    )
                    representative_microhaplotypes_old_index_key[pmo_index][
                        representative_microhaplotypes_index
                    ] = new_mhaps_target_index
                    representative_microhaplotypes_out_index_key[
                        pmo["target_info"][representative_microhaplotypes["target_id"]][
                            "target_name"
                        ]
                    ] = new_mhaps_target_index
                    for adding_microhap_index, adding_microhap in enumerate(
                        representative_microhaplotypes["microhaplotypes"]
                    ):
                        representative_microhaplotypes_hmap_for_target_index_old_index_key[
                            pmo_index
                        ][representative_microhaplotypes_index][
                            adding_microhap_index
                        ] = adding_microhap_index
        # print(representative_microhaplotypes_hmap_for_target_index_old_index_key)
        # update detected_microhaplotypes
        pmo_out["detected_microhaplotypes"] = copy.deepcopy(
            pmos[0]["detected_microhaplotypes"]
        )
        for pmo_index, pmo in enumerate(pmos[1:], start=1):
            # update indexes
            for detected_microhaplotypes in pmo["detected_microhaplotypes"]:
                detected_microhaplotypes_copy = copy.deepcopy(detected_microhaplotypes)
                for library_sample in detected_microhaplotypes_copy["library_samples"]:
                    for target_result in library_sample["target_results"]:
                        for hap in target_result["mhaps"]:
                            hap[
                                "mhap_id"
                            ] = representative_microhaplotypes_hmap_for_target_index_old_index_key[
                                pmo_index
                            ][target_result["mhaps_target_id"]][hap["mhap_id"]]
                        target_result[
                            "mhaps_target_id"
                        ] = representative_microhaplotypes_old_index_key[pmo_index][
                            target_result["mhaps_target_id"]
                        ]
                    library_sample[
                        "library_sample_id"
                    ] = library_sample_info_old_index_key[pmo_index][
                        library_sample["library_sample_id"]
                    ]
                detected_microhaplotypes_copy[
                    "bioinformatics_run_id"
                ] = bioinformatics_run_info_old_index_key[pmo_index][
                    detected_microhaplotypes_copy["bioinformatics_run_id"]
                ]
                # append after the indexes have been updated
                pmo_out["detected_microhaplotypes"].append(
                    detected_microhaplotypes_copy
                )

        pmo_indexes_with_read_counts_by_stage = []
        for pmo_index, pmo in enumerate(pmos):
            if "read_counts_by_stage" in pmo:
                pmo_indexes_with_read_counts_by_stage.append(pmo_index)
        if 0 not in pmo_indexes_with_read_counts_by_stage:
            pmo_out["read_counts_by_stage"] = []
        for pmo_index in pmo_indexes_with_read_counts_by_stage:
            # if read_counts_by_stage is in pmos[0] then no indexes need to be updated
            if 0 == pmo_index:
                pmo_out["read_counts_by_stage"] = copy.deepcopy(
                    pmos[pmo_index]["read_counts_by_stage"]
                )
            else:
                # update index and then append to out
                for read_counts_by_stage in pmos[pmo_index]["read_counts_by_stage"]:
                    read_counts_by_stage_copy = copy.deepcopy(read_counts_by_stage)
                    for (
                        read_counts_by_library_sample_by_stage
                    ) in read_counts_by_stage_copy[
                        "read_counts_by_library_sample_by_stage"
                    ]:
                        if (
                            "read_counts_for_targets"
                            in read_counts_by_library_sample_by_stage
                        ):
                            for (
                                read_counts_for_target
                            ) in read_counts_by_library_sample_by_stage[
                                "read_counts_for_targets"
                            ]:
                                read_counts_for_target[
                                    "target_id"
                                ] = target_info_old_index_key[pmo_index][
                                    read_counts_for_target["target_id"]
                                ]
                        read_counts_by_library_sample_by_stage[
                            "library_sample_id"
                        ] = library_sample_info_old_index_key[pmo_index][
                            read_counts_by_library_sample_by_stage["library_sample_id"]
                        ]
                    read_counts_by_stage_copy[
                        "bioinformatics_run_id"
                    ] = bioinformatics_run_info_old_index_key[pmo_index][
                        read_counts_by_stage_copy["bioinformatics_run_id"]
                    ]
                    pmo_out["read_counts_by_stage"].append(read_counts_by_stage_copy)
        return pmo_out
