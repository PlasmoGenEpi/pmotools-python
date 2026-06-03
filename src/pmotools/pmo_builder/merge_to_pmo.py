#!/usr/bin/env python3
from datetime import date
import numpy as np
from pmotools import __version__ as __pmotools_version__
from pmotools.pmo_builder.mhap_table_to_pmo import (
    create_minimum_library_specimen_dict_from_mhap_table,
)


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
    mhap_info: dict,
    panel_target_info: dict,
    specimen_info: list | None = None,
    library_sample_info: list | None = None,
    sequencing_info: list | None = None,
    bioinfo_method_info: list | None = None,
    bioinfo_run_info: list | None = None,
    project_info: list | None = None,
    read_counts_by_stage_info: list | None = None,
):
    """
    Merge components into PMO, replacing names with indices.

    The required input are ``mhap_info`` (must have fields:detected_microhaplotypes and representative_microhaplotypes) and ``panel_target_info`` (must have fields: target_info and panel_info). If no ``library_sample_info`` or ``specimen_info`` are provided, they will be automatically generated from the detected_microhaplotypes. Is also possible to provide only ``specimen_info`` or ``library_sample_info`` but their names must match up with the detected_microhaplotypes names.

    Args:
        mhap_info (dict): microhaplotypes within this project, both detected
            and representative; must contain ``detected_microhaplotypes`` and
            ``representative_microhaplotypes``.
        panel_target_info (dict): panel and target information; must contain
            ``target_info`` and ``panel_info``.
        specimen_info (list, optional): all the specimens within this project.
        library_sample_info (list, optional): library samples within this project.
        sequencing_info (list, optional): sequencing info for this project.
        bioinfo_method_info (list, optional): bioinformatics pipeline/methods.
        bioinfo_run_info (list, optional): runtime info for the pipeline.
        project_info (list, optional): info about projects stored in this PMO.
        read_counts_by_stage_info (list, optional): read counts by stage.

    Returns:
        str: a JSON-formatted PMO string.
    """

    missing_fields = []
    if "panel_info" not in panel_target_info:
        missing_fields.append("panel_info")
    if "target_info" not in panel_target_info:
        missing_fields.append("target_info")
    if "representative_microhaplotypes" not in mhap_info:
        missing_fields.append("representative_microhaplotypes")
    if "detected_microhaplotypes" not in mhap_info:
        missing_fields.append("detected_microhaplotypes")
    if missing_fields:
        raise ValueError(
            f"Missing required fields for panel_target_info or mhap_info: {missing_fields}"
        )

    if bioinfo_run_info is not None and bioinfo_method_info is None:
        raise ValueError(
            "bioinfo_method_info must be provided if bioinfo_run_info is provided"
        )
    if specimen_info is not None and library_sample_info is None:
        spec_and_lib_info = create_minimum_library_specimen_dict_from_mhap_table(
            mhap_info["detected_microhaplotypes"],
            panel_target_info["panel_info"][0]["panel_name"],
        )
        library_sample_info = spec_and_lib_info["library_sample_info"]

        # Validate that library_sample_info and specimen_info have matching names
        library_sample_names = {
            item["library_sample_name"] for item in library_sample_info
        }
        specimen_names = {item["specimen_name"] for item in specimen_info}

        # Check for names in library_sample_info that are not in specimen_info
        missing_in_specimen = library_sample_names - specimen_names
        if missing_in_specimen:
            raise ValueError(
                f"library_sample_names found in the detected_microhaplotypes that don't have corresponding supplied specimen_names: {sorted(missing_in_specimen)}"
            )

        # Check for names in specimen_info that are not in library_sample_info
        missing_in_library = specimen_names - library_sample_names
        if missing_in_library:
            raise ValueError(
                f"specimen_name were supplied that don't have corresponding library_sample_names in detected_microhaplotypes: {sorted(missing_in_library)}"
            )

        # names match, so the supplied specimen_info names will match up with the names in library_sample_info
    if specimen_info is None and library_sample_info is None:
        if len(panel_target_info["panel_info"]) > 1:
            raise Exception(
                "If multiple panels are included in the panel information,specimen_info and library_sample_info must also be provided to indicate which panel was used for each specimen. Otherwise, provide only a single panel. Panels found: "
                + str(len(panel_target_info["panel_info"]))
            )
        spec_and_lib_info = create_minimum_library_specimen_dict_from_mhap_table(
            mhap_info["detected_microhaplotypes"],
            panel_target_info["panel_info"][0]["panel_name"],
        )
        specimen_info = spec_and_lib_info["specimen_info"]
        library_sample_info = spec_and_lib_info["library_sample_info"]
    else:
        # Make copies to avoid editing input
        library_sample_info = [dict(d) for d in library_sample_info]
        if specimen_info is not None:
            specimen_info = [dict(d) for d in specimen_info]
        else:
            # if giving only library sample info can default to the specimen being just the library_sample_names
            for library_sample in library_sample_info:
                library_sample["specimen_name"] = library_sample["library_sample_name"]
            specimen_info = [
                {"specimen_name": library_sample["library_sample_name"]}
                for library_sample in library_sample_info
            ]

    panel_target_info = _convert_numpy_scalars(panel_target_info)
    mhap_info = _convert_numpy_scalars(mhap_info)
    # optional
    if sequencing_info is not None:
        sequencing_info = [dict(d) for d in sequencing_info]
        sequencing_info = _convert_numpy_scalars(sequencing_info)
    if bioinfo_method_info is not None:
        bioinfo_method_info = [dict(d) for d in bioinfo_method_info]
        bioinfo_method_info = _convert_numpy_scalars(bioinfo_method_info)
    if bioinfo_run_info is not None:
        bioinfo_run_info = [dict(d) for d in bioinfo_run_info]
        bioinfo_run_info = _convert_numpy_scalars(bioinfo_run_info)
    if project_info is not None:
        project_info = [dict(d) for d in project_info]
        project_info = _convert_numpy_scalars(project_info)

    # Handle read_counts_by_stage_info if provided
    if read_counts_by_stage_info is not None:
        read_counts_by_stage_info = [
            _convert_numpy_scalars(dict(d)) for d in read_counts_by_stage_info
        ]

    specimen_info = _convert_numpy_scalars(specimen_info)
    library_sample_info = _convert_numpy_scalars(library_sample_info)

    _replace_names_with_IDs(
        specimen_info=specimen_info,
        project_info=project_info,
        library_sample_info=library_sample_info,
        sequencing_info=sequencing_info,
        panel_target_info=panel_target_info,
        mhap_info=mhap_info,
        bioinfo_run_info=bioinfo_run_info,
        read_counts_by_stage_info=read_counts_by_stage_info,
    )

    # Build PMO
    pmo_header = _generate_pmo_header()
    pmo = (
        {
            "pmo_header": pmo_header,
            "library_sample_info": library_sample_info,
            "specimen_info": specimen_info,
        }
        | panel_target_info
        | mhap_info
    )

    if sequencing_info:
        pmo["sequencing_info"] = sequencing_info
    if bioinfo_method_info:
        pmo["bioinformatics_methods_info"] = bioinfo_method_info
    if bioinfo_run_info:
        pmo["bioinformatics_run_info"] = bioinfo_run_info
    if project_info:
        pmo["project_info"] = project_info

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


def _generate_pmo_header(version=__pmotools_version__):
    today = date.today().isoformat()
    # TODO: update to grab pmo version - will put this in a seperate PR
    pmo_header = {
        "pmo_version": version,
        "creation_date": today,
        "generation_method": {
            "program_name": "pmotools-python",
            "program_version": version,
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
    specimen_info: list,
    panel_target_info: dict,
    mhap_info: dict,
    project_info: list | None = None,
    library_sample_info: list | None = None,
    sequencing_info: list | None = None,
    bioinfo_run_info: list | None = None,
    read_counts_by_stage_info: list | None = None,
):
    # SPECIMEN INFO
    # replace name with project ID
    if project_info is not None:
        missing_projects = _replace_key_with_id(
            specimen_info, project_info, "project_name", "project_id"
        )
    else:
        missing_projects = []

    # LIBRARY SAMPLE INFO
    # replace with sequencing_info_id, specimen_id, panel_id
    missing_specimen = _replace_key_with_id(
        library_sample_info, specimen_info, "specimen_name", "specimen_id"
    )
    missing_panels = _replace_key_with_id(
        library_sample_info,
        panel_target_info["panel_info"],
        "panel_name",
        "panel_id",
    )
    if sequencing_info is not None:
        missing_sequencing = _replace_key_with_id(
            library_sample_info,
            sequencing_info,
            "sequencing_info_name",
            "sequencing_info_id",
        )
    else:
        missing_sequencing = []

    # REP MHAPS
    # replace target_name with ID
    missing_targets = _replace_key_with_id(
        mhap_info["representative_microhaplotypes"]["targets"],
        panel_target_info["target_info"],
        "target_name",
        "target_id",
    )

    # DETECTED MHAPS
    # Replace library_sample_name and bioinformatics_run_name
    if bioinfo_run_info is not None:
        missing_bioinfo_runs = _replace_key_with_id(
            mhap_info["detected_microhaplotypes"],
            bioinfo_run_info,
            "bioinformatics_run_name",
            "bioinformatics_run_id",
        )
    else:
        missing_bioinfo_runs = []
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
    target_lookup = _make_lookup(panel_target_info["target_info"], "target_name")
    if read_counts_by_stage_info is not None:
        if bioinfo_run_info is not None:
            # Replace bioinformatics_run_name with bioinformatics_run_id
            missing_read_counts_bioinfo_runs = _replace_key_with_id(
                read_counts_by_stage_info,
                bioinfo_run_info,
                "bioinformatics_run_name",
                "bioinformatics_run_id",
            )
        else:
            missing_read_counts_bioinfo_runs = []

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
