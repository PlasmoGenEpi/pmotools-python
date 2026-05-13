#!/usr/bin/env python3
import json
from pathlib import Path

import pandas as pd

from pmotools.pmo_builder.metatable_to_pmo import (
    specimen_info_table_to_pmo,
)
from pmotools.pmo_builder.panel_information_to_pmo import panel_info_table_to_pmo
from pmotools.pmo_builder.mhap_table_to_pmo import (
    mhap_table_to_pmo,
    create_minimum_library_specimen_dict_from_mhap_table,
)
from pmotools.pmo_builder.pmo_updater import PMOUpdater
from pmotools.pmo_builder.merge_to_pmo import merge_to_pmo
from pmotools.pmo_engine.pmo_checker import PMOChecker


def test_minimal_pmo_creation():
    """Build PMO just from the smallest amount of data from  allele table and target/panel info."""

    # Panel and target information with optional fields
    target_df = pd.DataFrame(
        {
            "target_name": ["target1"],
            "fwd_primer": ["ATGCATGC"],
            "rev_primer": ["GCATGCAT"],
        }
    )
    panel_and_target_info = panel_info_table_to_pmo(
        target_table=target_df,
        panel_name="panel1",
        target_name_col="target_name",
        forward_primers_seq_col="fwd_primer",
        reverse_primers_seq_col="rev_primer",
    )

    # Microhaplotype information with optional details
    mhap_df = pd.DataFrame(
        {
            "library_sample_name": ["lib1"],
            "target_name": ["target1"],
            "seq": ["ATGCATGC"],
            "reads": [42],
        }
    )

    mhap_info = mhap_table_to_pmo(
        microhaplotype_table=mhap_df,
        library_sample_name_col="library_sample_name",
        target_name_col="target_name",
        seq_col="seq",
        reads_col="reads",
    )

    lib_and_spec_infos = create_minimum_library_specimen_dict_from_mhap_table(
        mhap_info["detected_microhaplotypes"],
        panel_name="panel1",
        library_sample_specimen_key={"lib1": "specimen1"},
    )
    # Merge into PMO structure
    no_spec_meta_pmo = merge_to_pmo(
        specimen_info=lib_and_spec_infos["specimen_info"],
        library_sample_info=lib_and_spec_infos["library_sample_info"],
        panel_and_target_info=panel_and_target_info,
        mhap_info=mhap_info,
    )

    # Specimen information with optional fields
    specimen_df = pd.DataFrame(
        {
            "specimen_name": ["specimen1"],
            "specimen_taxon_id": [[5900]],
            "host_taxon_id": [9606],
            "collection_date": ["2024-01-01"],
            "collection_country": ["Wonderland"],
            "project_name": ["Test Project"],
            "host_age": [35],
            "host_sex": ["female"],
            "lat_lon": ["37.77,-122.42"],
            "specimen_collect_device": ["venipuncture"],
            "specimen_comments": [["no issues"]],
            "specimen_store_loc": ["Freezer 1"],
            "drug_usage": [["DrugX"]],
            "env_broad_scale": ["Urban"],
            "env_local_scale": ["Clinic"],
            "env_medium": ["Blood"],
            "alternate_ids": [["ALT1", "ALT2"]],
            "custom_note": ["Important specimen"],
            "parasite_density": [1200],
            "parasite_density_method": ["microscopy"],
        }
    )

    specimen_info = specimen_info_table_to_pmo(
        specimen_df,
        specimen_name_col="specimen_name",
        specimen_taxon_id_col="specimen_taxon_id",
        host_taxon_id_col="host_taxon_id",
        collection_date_col="collection_date",
        collection_country_col="collection_country",
        project_name_col="project_name",
        alternate_identifiers_col="alternate_ids",
        drug_usage_col="drug_usage",
        env_broad_scale_col="env_broad_scale",
        env_local_scale_col="env_local_scale",
        env_medium_col="env_medium",
        host_age_col="host_age",
        host_sex_col="host_sex",
        specimen_collect_device_col="specimen_collect_device",
        specimen_comments_col="specimen_comments",
        specimen_store_loc_col="specimen_store_loc",
        lat_lon_col="lat_lon",
        parasite_density_col="parasite_density",
        parasite_density_method_col="parasite_density_method",
        additional_specimen_cols=["custom_note"],
    )

    # merging in specimen meta and merge into PMO structure
    with_spec_meta_pmo = merge_to_pmo(
        specimen_info=PMOUpdater.merge_dicts_by_key(
            lib_and_spec_infos["specimen_info"], specimen_info, "specimen_name"
        ),
        library_sample_info=lib_and_spec_infos["library_sample_info"],
        panel_and_target_info=panel_and_target_info,
        mhap_info=mhap_info,
    )

    # Load the schema and validate using PMOChecker
    # checking against both versions of PMO, the above builds a "full" PMO and want to check if new schema still validates
    # this old format
    for schema_version in ["1.1.0"]:
        schemas_dir = (
            Path(__file__).resolve().parents[2] / "src" / "pmotools" / "schemas"
        )
        schema_filename = (
            f"portable_microhaplotype_object_v{schema_version}.schema.json"
        )
        schema_path = schemas_dir / schema_filename
        if not schema_path.exists():
            available_schemas = sorted(
                schemas_dir.glob("portable_microhaplotype_object_*.schema.json")
            )
            if not available_schemas:
                raise FileNotFoundError(
                    f"No schema files found in {schemas_dir} matching "
                    "'portable_microhaplotype_object_*.schema.json'"
                )
            schema_path = available_schemas[-1]
        with schema_path.open(encoding="utf-8") as schema_file:
            schema = json.load(schema_file)

        checker = PMOChecker(schema)
        checker.check_for_required_base_fields(no_spec_meta_pmo)
        checker.validate_pmo_json(no_spec_meta_pmo)
        checker.check_for_required_base_fields(with_spec_meta_pmo)
        checker.validate_pmo_json(with_spec_meta_pmo)

    # Validate optional fields propagated through update function
    specimen_entry = with_spec_meta_pmo["specimen_info"][0]
    assert specimen_entry["collection_country"] == "Wonderland"
    assert specimen_entry["collection_date"] == "2024-01-01"
    assert specimen_entry["host_age"] == 35
    assert specimen_entry["specimen_store_loc"] == "Freezer 1"
    assert specimen_entry["custom_note"] == "Important specimen"
    assert specimen_entry["parasite_density_info"][0]["parasite_density"] == 1200
    assert (
        specimen_entry["parasite_density_info"][0]["parasite_density_method"]
        == "microscopy"
    )
