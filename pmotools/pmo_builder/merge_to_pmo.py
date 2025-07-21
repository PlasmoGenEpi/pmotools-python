#!/usr/bin/env python3
import pandas as pd
import numpy as np
from datetime import date
import json
import warnings

from ..pmo_builder.json_convert_utils import check_additional_columns_exist
from ..pmo_engine.pmo_processor import PMOProcessor


def panel_info_table_to_json(
    specimen_info: list,
    library_info: list,
    sequencing_info: list,
    panel_info: list,
    mhap_info: dict,
    bioinfo_method_info: list,
    bioinfo_run_info: list,
    project_info: list,
):
    # TODO: make sure it isn't editing the original
    # TODO: write function to generate look ups
    # TODO: add checks if the entry doesn't exist in the lookup table
    # TODO: decide whether to take json or dict

    # SPECIMEN INFO
    # replace name with project ID
    project_lookup = {proj["project_name"]
        : idx for idx, proj in enumerate(project_info)}
    for specimen_dict in specimen_info:
        name = specimen_dict["project_name"]
        specimen_dict["project_id"] = project_lookup.get(name)
        del specimen_dict["project_name"]

    # EXPERIMENT INFO
    # replace with sequencing_info_id, specimen_id, panel_id
    seq_lookup = {seq_run['sequencing_info_name']
        : idx for idx, seq_run in enumerate(sequencing_info)}
    spec_lookup = {specimen['specimen_name']
        : idx for idx, specimen in enumerate(specimen_info)}
    panel_lookup = {panel['panel_name']: idx for idx,
                    panel in enumerate(panel_info["panel_info"])}
    for library_dict in library_info:
        library_dict['sequencing_info_id'] = seq_lookup.get(
            library_dict["sequencing_info_name"])
        library_dict['specimen_id'] = spec_lookup.get(
            library_dict["specimen_name"])
        library_dict['panel_id'] = panel_lookup.get(library_dict["panel_name"])
        del library_dict["sequencing_info_name"]
        del library_dict["specimen_name"]
        del library_dict["panel_name"]

    # REP MHAPS
    # target id
    target_lookup = {target["target_name"]: idx for idx,
                     target in enumerate(panel_info["target_info"])}
    rep_dicts = mhap_info["representative_microhaplotypes"]['targets']
    for target_dict in rep_dicts:
        target_dict["target_id"] = target_lookup.get(
            str(target_dict["target_name"]))
        del target_dict["target_name"]

    # DETECTED MHAPS
    # set bioinformatics_run_id and library_sample_id
    bioinfo_run_lookup = {run['bioinformatics_run_name']: idx for idx, run in enumerate(bioinfo_run_info)}
    library_sample_lookup = {
        lib_sample['library_sample_name']: idx for idx, lib_sample in enumerate(library_info)}
    detected_dicts = mhap_info["detected_microhaplotypes"]
    for detected_dict in detected_dicts:
        detected_dict["bioinformatics_run_id"] = bioinfo_run_lookup.get(
            detected_dict["bioinformatics_run_name"])
        del detected_dict["bioinformatics_run_name"]
        for sample_dict in detected_dict["library_samples"]:
            sample_dict["library_sample_id"] = library_sample_lookup.get(
                sample_dict["library_sample_name"])
            del sample_dict["library_sample_name"]

    pmo_header = generate_pmo_header()
    pmo = {
        "library_sample_info": library_info,
        "specimen_info": specimen_info,
        "sequencing_info": sequencing_info,
        "bioinformatics_methods_info": bioinfo_method_info,
        "bioinformatics_run_info": bioinfo_run_info,
        "project_info": project_info,
        "pmo_header": pmo_header
    } | panel_info | mhap_info

    return pmo


def generate_pmo_header():
    today = date.today().isoformat()
    # TODO: update to grab pmo version
    pmo_header = {'pmo_version': '1.0.0', 'creation_date': today, 'generation_method': {
        'program_name': 'pmotools-python', 'program_version': 'pmo_version'}}
    return pmo_header
