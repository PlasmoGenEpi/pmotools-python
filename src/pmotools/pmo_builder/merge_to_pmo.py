#!/usr/bin/env python3
from datetime import date
import numpy as np
import json
from pmotools import __version__ as __pmotools_version__

def _convert_numpy_scalars(obj):
    """Recursively convert numpy scalar types to native Python types."""
    if isinstance(obj, dict):
        return {key: _convert_numpy_scalars(value) for key, value in obj.items()}
    if isinstance(obj, list):
        return [_convert_numpy_scalars(value) for value in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


def merge_to_pmo(
    specimen_info: list,
    library_sample_info: list,
    sequencing_info: list,
    panel_info: dict,
    mhap_info: dict,
    bioinfo_method_info: list,
    bioinfo_run_info: list,
    project_info: list,
    read_counts_by_stage_info: list | None = None,
):
    """
    Merge components into PMO, replacing names with indeces.

    :param specimen_info (list): a list of all the specimens within this project
    :param library_sample_info (list) : a list of library samples within this project
    :param sequencing_info (list) : a list of sequencing info for this project
    :param panel_info (list) : a dictionary containing the panel and target information for this project
    :param mhap_info (list) : a dictionary containing the microhaplotypes within this project, both detected and representative
    :param bioinfo_method_info (list) : the bioinformatics pipeline/methods used to generated the amplicon analysis for this project
    :param bioinfo_run_info (list) : the runtime info for the bioinformatics pipeline used to generated the amplicon analysis for this project
    :param project_info (list) : the information about the projects stored in this PMO
    :param read_counts_by_stage_info (Optional[list]) : the read counts by stage information for this project

    :return: a json formatted PMO string.
    """
    # Make copies to avoid editing input
    specimen_info = [dict(d) for d in specimen_info]
    library_sample_info = [dict(d) for d in library_sample_info]
    sequencing_info = [dict(d) for d in sequencing_info]
    bioinfo_method_info = [dict(d) for d in bioinfo_method_info]
    bioinfo_run_info = [dict(d) for d in bioinfo_run_info]
    project_info = [dict(d) for d in project_info]
    panel_info = _convert_numpy_scalars(panel_info)
    mhap_info = _convert_numpy_scalars(mhap_info)

    # Handle read_counts_by_stage_info if provided
    if read_counts_by_stage_info is not None:
        read_counts_by_stage_info = [
            _convert_numpy_scalars(dict(d)) for d in read_counts_by_stage_info
        ]

    specimen_info = _convert_numpy_scalars(specimen_info)
    library_sample_info = _convert_numpy_scalars(library_sample_info)
    sequencing_info = _convert_numpy_scalars(sequencing_info)
    bioinfo_method_info = _convert_numpy_scalars(bioinfo_method_info)
    bioinfo_run_info = _convert_numpy_scalars(bioinfo_run_info)
    project_info = _convert_numpy_scalars(project_info)

    _replace_names_with_IDs(
        specimen_info,
        project_info,
        library_sample_info,
        sequencing_info,
        panel_info,
        mhap_info,
        bioinfo_run_info,
        read_counts_by_stage_info,
    )

    # Build PMO
    pmo_header = _generate_pmo_header()
    pmo = (
        {
            "pmo_header": pmo_header,
            "library_sample_info": library_sample_info,
            "specimen_info": specimen_info,
            "sequencing_info": sequencing_info,
            "bioinformatics_methods_info": bioinfo_method_info,
            "bioinformatics_run_info": bioinfo_run_info,
            "project_info": project_info,
        }
        | panel_info
        | mhap_info
    )

    # Add read_counts_by_stage_info if provided
    if read_counts_by_stage_info is not None:
        pmo["read_counts_by_stage"] = read_counts_by_stage_info

    return _convert_numpy_scalars(pmo)


def _make_lookup(dict, key):
    lookup = {entry[key]: idx for idx, entry in enumerate(dict)}
    return lookup


def _replace_key_with_id(target_list, reference_list, name_key, id_key, lookup=None):
    """
    Replaces name_key in target_list with id_key, based on lookup from reference_list.
    """
    if not lookup:
        lookup = _make_lookup(reference_list, name_key)
    unique_names = set()
    for entry in target_list:
        name = str(entry.pop(name_key))
        unique_names.add(name)
        entry[id_key] = lookup.get(name)
    missing_items = list(unique_names - lookup.keys())
    return missing_items


def _generate_pmo_header():
    today = date.today().isoformat()
    # TODO: update to grab pmo version - will put this in a seperate PR
    pmo_header = {
        "pmo_version": __pmotools_version__,
        "creation_date": today,
        "generation_method": {
            "program_name": "pmotools-python",
            "program_version": __pmotools_version__,
        },
    }
    return pmo_header


def _report_missing_IDs(
    missing_projects,
    missing_sequencing,
    missing_specimen,
    missing_panels,
    missing_targets,
    missing_bioinfo_runs,
    missing_libs,
    missing_read_counts_bioinfo_runs,
    missing_read_counts_libs,
    missing_read_counts_targets,
):
    if any(
        [
            missing_projects,
            missing_sequencing,
            missing_specimen,
            missing_panels,
            missing_targets,
            missing_bioinfo_runs,
            missing_libs,
            missing_read_counts_bioinfo_runs,
            missing_read_counts_libs,
            missing_read_counts_targets,
        ]
    ):
        error_message = (
            "The following fields were found in one table and not another:\n"
        )
        if missing_projects:
            error_message += f"Project names in Specimen Info not in Project Info: {missing_projects}\n"
        if missing_sequencing:
            error_message += f"Sequencing names in Library Sample Info not in Sequencing Info: {missing_sequencing}\n"
        if missing_specimen:
            error_message += f"Specimen names in Library Sample Info not in Specimen Info: {missing_specimen}\n"
        if missing_panels:
            error_message += f"Panel names in Library Sample Info not in Panel Info: {missing_panels}\n"
        if missing_targets:
            error_message += f"Target names in Representative Microhaplotypes not in Target Info: {missing_targets}\n"
        if missing_bioinfo_runs:
            error_message += f"Bioinformatics run names in Detected Microhaplotypes not in Bioinformatic Run Info: {missing_bioinfo_runs}\n"
        if missing_libs:
            error_message += f"Library Sample names in Detected Microhaplotypes not in Library Sample Info: {missing_libs}\n"
        if missing_read_counts_bioinfo_runs:
            error_message += f"Bioinformatics run names in Read Counts by Stage not in Bioinformatic Run Info: {missing_read_counts_bioinfo_runs}\n"
        if missing_read_counts_libs:
            error_message += f"Library Sample names in Read Counts by Stage not in Library Sample Info: {missing_read_counts_libs}\n"
        if missing_read_counts_targets:
            error_message += f"Target names in Read Counts by Stage not in Target Info: {missing_read_counts_targets}\n"
        raise ValueError(error_message)


def _replace_names_with_IDs(
    specimen_info,
    project_info,
    library_sample_info,
    sequencing_info,
    panel_info,
    mhap_info,
    bioinfo_run_info,
    read_counts_by_stage_info,
):
    # SPECIMEN INFO
    # replace name with project ID
    missing_projects = _replace_key_with_id(
        specimen_info, project_info, "project_name", "project_id"
    )

    # LIBRARY SAMPLE INFO
    # replace with sequencing_info_id, specimen_id, panel_id
    missing_sequencing = _replace_key_with_id(
        library_sample_info,
        sequencing_info,
        "sequencing_info_name",
        "sequencing_info_id",
    )
    missing_specimen = _replace_key_with_id(
        library_sample_info, specimen_info, "specimen_name", "specimen_id"
    )
    missing_panels = _replace_key_with_id(
        library_sample_info, panel_info["panel_info"], "panel_name", "panel_id"
    )

    # REP MHAPS
    # replace target_name with ID
    missing_targets = _replace_key_with_id(
        mhap_info["representative_microhaplotypes"]["targets"],
        panel_info["target_info"],
        "target_name",
        "target_id",
    )

    # DETECTED MHAPS
    # Replace library_sample_name and bioinformatics_run_name
    missing_bioinfo_runs = _replace_key_with_id(
        mhap_info["detected_microhaplotypes"],
        bioinfo_run_info,
        "bioinformatics_run_name",
        "bioinformatics_run_id",
    )
    lib_sample_lookup = _make_lookup(library_sample_info, "library_sample_name")
    missing_libs = []
    for detected in mhap_info["detected_microhaplotypes"]:
        missing_libs += _replace_key_with_id(
            detected["library_samples"],
            library_sample_info,
            "library_sample_name",
            "library_sample_id",
            lookup=lib_sample_lookup,
        )

    # READ COUNTS BY STAGE
    # Replace bioinformatics_run_name and library_sample_name if provided
    missing_read_counts_bioinfo_runs = []
    missing_read_counts_libs = []
    missing_read_counts_targets = []
    target_lookup = _make_lookup(panel_info["target_info"], "target_name")
    if read_counts_by_stage_info is not None:
        # Replace bioinformatics_run_name with bioinformatics_run_id
        missing_read_counts_bioinfo_runs = _replace_key_with_id(
            read_counts_by_stage_info,
            bioinfo_run_info,
            "bioinformatics_run_name",
            "bioinformatics_run_id",
        )

        # Replace library_sample_name with library_sample_id in each run and map targets
        for read_counts_run in read_counts_by_stage_info:
            missing_read_counts_libs += _replace_key_with_id(
                read_counts_run["read_counts_by_library_sample_by_stage"],
                library_sample_info,
                "library_sample_name",
                "library_sample_id",
                lookup=lib_sample_lookup,
            )

            for library_entry in read_counts_run.get(
                "read_counts_by_library_sample_by_stage", []
            ):
                target_entries = library_entry.get("read_counts_for_targets") or []
                for target_entry in target_entries:
                    target_name = target_entry.pop("target_name", None)
                    if target_name is None:
                        continue
                    if target_name in target_lookup:
                        target_entry["target_id"] = target_lookup[target_name]
                    else:
                        missing_read_counts_targets.append(target_name)

    # If any names were missing from reference tables error
    _report_missing_IDs(
        missing_projects,
        missing_sequencing,
        missing_specimen,
        missing_panels,
        missing_targets,
        missing_bioinfo_runs,
        missing_libs,
        missing_read_counts_bioinfo_runs,
        missing_read_counts_libs,
        missing_read_counts_targets,
    )
