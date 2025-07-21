#!/usr/bin/env python3
from datetime import date
import json
import warnings


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
    # Make copies to avoid editing input
    specimen_info = [dict(d) for d in specimen_info]
    library_info = [dict(d) for d in library_info]
    sequencing_info = [dict(d) for d in sequencing_info]
    bioinfo_method_info = [dict(d) for d in bioinfo_method_info]
    bioinfo_run_info = [dict(d) for d in bioinfo_run_info]
    project_info = [dict(d) for d in project_info]
    panel_info = json.loads(json.dumps(panel_info))
    mhap_info = json.loads(json.dumps(mhap_info))
    # TODO: add checks if the entry doesn't exist in the lookup table
    # TODO: decide whether to take json or dict

    # SPECIMEN INFO
    # replace name with project ID
    replace_key_with_id(specimen_info, project_info,
                        "project_name", "project_id")

    # EXPERIMENT INFO
    # replace with sequencing_info_id, specimen_id, panel_id
    replace_key_with_id(library_info, sequencing_info,
                        'sequencing_info_name', "sequencing_info_id")
    replace_key_with_id(library_info, specimen_info,
                        'specimen_name', "specimen_id")
    replace_key_with_id(
        library_info, panel_info["panel_info"], 'panel_name', "panel_id")

    # REP MHAPS
    # replace target_name with ID
    replace_key_with_id(
        mhap_info["representative_microhaplotypes"]["targets"],
        panel_info["target_info"],
        "target_name",
        "target_id"
    )

    # DETECTED MHAPS
    # Replace library_sample_name and bioinformatics_run_name
    replace_key_with_id(
        mhap_info["detected_microhaplotypes"],
        bioinfo_run_info,
        "bioinformatics_run_name",
        "bioinformatics_run_id"
    )
    lib_sample_lookup = make_lookup(library_info, 'library_sample_name')
    for detected in mhap_info["detected_microhaplotypes"]:
        replace_key_with_id(detected["library_samples"],
                            library_info, 'library_sample_name', 'library_sample_id', lookup=lib_sample_lookup)

    # Build PMO
    pmo_header = generate_pmo_header()
    pmo = {
        "pmo_header": pmo_header,
        "library_sample_info": library_info,
        "specimen_info": specimen_info,
        "sequencing_info": sequencing_info,
        "bioinformatics_methods_info": bioinfo_method_info,
        "bioinformatics_run_info": bioinfo_run_info,
        "project_info": project_info,
    } | panel_info | mhap_info
    pmo_json = json.dumps(pmo, indent=4)
    return pmo_json


def make_lookup(dict, key):
    lookup = {entry[key]: idx for idx,
              entry in enumerate(dict)}
    return lookup


def replace_key_with_id(
    target_list,
    reference_list,
    name_key,
    id_key,
    lookup=None
):
    """
    Replaces name_key in target_list with id_key, based on lookup from reference_list.
    """
    if not lookup:
        lookup = make_lookup(reference_list, name_key)
    for entry in target_list:
        entry[id_key] = lookup.get(str(entry.pop(name_key)))


def generate_pmo_header():
    today = date.today().isoformat()
    # TODO: update to grab pmo version
    pmo_header = {'pmo_version': '1.0.0', 'creation_date': today, 'generation_method': {
        'program_name': 'pmotools-python', 'program_version': 'pmo_version'}}
    return pmo_header
