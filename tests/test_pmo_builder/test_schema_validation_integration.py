#!/usr/bin/env python3
import json
from pathlib import Path

import pandas as pd

from pmotools import __version__
from pmotools.pmo_builder.metatable_to_pmo import (
    library_sample_info_table_to_pmo,
    specimen_info_table_to_pmo,
)
from pmotools.pmo_builder.panel_information_to_pmo import panel_info_table_to_pmo
from pmotools.pmo_builder.mhap_table_to_pmo import mhap_table_to_pmo
from pmotools.pmo_builder.read_count_by_stage_table_to_pmo import (
    read_count_by_stage_table_to_pmo,
)
from pmotools.pmo_builder.merge_to_pmo import merge_to_pmo
from pmotools.pmo_engine.pmo_checker import PMOChecker


def test_toy_pmo_validates_against_schema():
    """Build a toy PMO with builder functions and validate against the schema."""
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

    # Sequencing information (manually constructed to satisfy schema requirements)
    sequencing_info = [
        {
            "sequencing_info_name": "seq_run1",
            "seq_platform": "ILLUMINA",
            "seq_instrument_model": "NextSeq 2000",
            "library_layout": "PAIRED",
            "library_strategy": "AMPLICON",
            "library_source": "GENOMIC",
            "library_selection": "PCR",
            "library_kit": "MiSeq Reagent Kit v3",
            "seq_center": "UCSF Core",
            "seq_date": "2024-01-15",
            "nucl_acid_amp": "https://example.org/amp",
            "nucl_acid_ext": "https://example.org/ext",
        }
    ]

    # Library sample information with optional fields
    library_df = pd.DataFrame(
        {
            "library_sample_name": ["lib1"],
            "sequencing_info_name": ["seq_run1"],
            "specimen_name": ["specimen1"],
            "panel_name": ["panel1"],
            "accession": ["ACC123"],
            "prep_plate_name": ["PlateA"],
            "prep_plate_row": ["B"],
            "prep_plate_col": [3],
            "library_note": ["High quality"],
        }
    )
    library_sample_info = library_sample_info_table_to_pmo(
        library_df,
        accession_col="accession",
        library_prep_plate_name_col="prep_plate_name",
        library_prep_plate_row_col="prep_plate_row",
        library_prep_plate_col_col="prep_plate_col",
        additional_library_sample_info_cols=["library_note"],
    )

    # Panel and target information with optional fields
    target_df = pd.DataFrame(
        {
            "target_name": ["target1"],
            "fwd_primer": ["ATGCATGC"],
            "rev_primer": ["GCATGCAT"],
            "reaction": ["rxn1"],
            "forward_start": [100],
            "forward_end": [107],
            "reverse_start": [150],
            "reverse_end": [157],
            "insert_start": [110],
            "insert_end": [149],
            "chrom": ["chr1"],
            "strand": ["+"],
            "ref_seq": ["ATGCGCTA"],
            "gene_name": ["geneA"],
            "target_attributes": [["marker"]],
            "amplicon_length": [250],
        }
    )
    genome_info = {
        "name": "TestGenome",
        "genome_version": "v1.0",
        "taxon_id": [9999],
        "url": "https://example.org/genome",
        "chromosomes": ["chr1"],
        "gff_url": "https://example.org/genome.gff",
    }
    panel_info = panel_info_table_to_pmo(
        target_table=target_df,
        panel_name="panel1",
        genome_info=genome_info,
        reaction_name_col="reaction",
        forward_primers_start_col="forward_start",
        forward_primers_end_col="forward_end",
        reverse_primers_start_col="reverse_start",
        reverse_primers_end_col="reverse_end",
        insert_start_col="insert_start",
        insert_end_col="insert_end",
        chrom_col="chrom",
        strand_col="strand",
        ref_seq_col="ref_seq",
        gene_name_col="gene_name",
        target_attributes_col="target_attributes",
        additional_target_info_cols=["amplicon_length"],
    )

    # Microhaplotype information with optional details
    mhap_df = pd.DataFrame(
        {
            "library_sample_name": ["lib1"],
            "target_name": ["target1"],
            "seq": ["ATGCATGC"],
            "reads": [42],
            "bioinformatics_run_name": ["run1"],
            "umis": [10],
            "microhap_name": ["mh1"],
            "pseudocigar": ["8M"],
            "quality": ["ABCD"],
            "mask_start": ["1"],
            "mask_segment": ["2"],
            "mask_replacement": ["2"],
            "chrom": ["chr1"],
            "start": [120],
            "end": [127],
            "ref_seq": ["ATGC"],
            "strand": ["+"],
            "custom_annotation": ["custom"],
            "custom_detected": ["det-note"],
        }
    )
    mhap_info = mhap_table_to_pmo(
        microhaplotype_table=mhap_df,
        bioinformatics_run_name="bioinformatics_run_name",
        umis_col="umis",
        microhaplotype_name_col="microhap_name",
        pseudocigar_col="pseudocigar",
        quality_col="quality",
        masking_seq_start_col="mask_start",
        masking_seq_segment_size_col="mask_segment",
        masking_replacement_size_col="mask_replacement",
        chrom_col="chrom",
        start_col="start",
        end_col="end",
        ref_seq_col="ref_seq",
        strand_col="strand",
        additional_representative_mhap_cols=["custom_annotation"],
        additional_mhap_detected_cols=["custom_detected"],
    )

    # Read counts by stage with additional columns
    total_raw_counts_df = pd.DataFrame(
        {
            "library_sample_name": ["lib1"],
            "total_raw_count": [1000],
            "bioinformatics_run_name": ["run1"],
            "library_comment": ["primary run"],
        }
    )
    reads_by_stage_df = pd.DataFrame(
        {
            "library_sample_name": ["lib1"],
            "target_name": ["target1"],
            "stage": ["RawReads"],
            "read_count": [800],
            "bioinformatics_run_name": ["run1"],
            "quality_score": [0.98],
            "coverage_depth": [150],
        }
    )
    read_counts_by_stage_info = read_count_by_stage_table_to_pmo(
        bioinformatics_run_name="bioinformatics_run_name",
        total_raw_count_table=total_raw_counts_df,
        reads_by_stage_table=reads_by_stage_df,
        additional_library_sample_cols=["library_comment"],
        additional_target_cols=["quality_score", "coverage_depth"],
    )

    # Bioinformatics method and run information
    bioinfo_methods_info = [
        {
            "methods": [
                {
                    "program": "Cutadapt",
                    "program_version": "v1.0.0",
                },
                {
                    "program": "DADA2",
                    "program_version": "v1.0.0",
                },
                {
                    "program": "CustomFilter",
                    "program_version": "v2.0.0",
                },
            ]
        }
    ]
    bioinfo_run_info = [
        {
            "bioinformatics_run_name": "run1",
            "bioinformatics_methods_id": 0,
            "run_date": "2024-02-01",
            "run_comments": "First batch",
        }
    ]

    # Project info with optional metadata
    project_info = [
        {
            "project_name": "Test Project",
            "project_description": "Toy project for schema validation.",
            "project_type": "Surveillance",
            "project_contributors": ["Alice", "Bob"],
            "BioProject_accession": "PRJNA12345",
        }
    ]

    # Merge into PMO structure
    pmo = merge_to_pmo(
        specimen_info=specimen_info,
        library_sample_info=library_sample_info,
        sequencing_info=sequencing_info,
        panel_info=panel_info,
        mhap_info=mhap_info,
        bioinfo_method_info=bioinfo_methods_info,
        bioinfo_run_info=bioinfo_run_info,
        project_info=project_info,
        read_counts_by_stage_info=read_counts_by_stage_info,
    )

    # Load the schema and validate using PMOChecker
    schemas_dir = Path(__file__).resolve().parents[2] / "src" / "pmotools" / "schemas"
    schema_filename = f"portable_microhaplotype_object_v{__version__}.schema.json"
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
    checker.check_for_required_base_fields(pmo)
    checker.validate_pmo_json(pmo)

    # Validate optional fields propagated through builders
    specimen_entry = pmo["specimen_info"][0]
    assert specimen_entry["host_age"] == 35
    assert specimen_entry["specimen_store_loc"] == "Freezer 1"
    assert specimen_entry["custom_note"] == "Important specimen"
    assert specimen_entry["parasite_density_info"][0]["parasite_density"] == 1200
    assert (
        specimen_entry["parasite_density_info"][0]["parasite_density_method"]
        == "microscopy"
    )

    sequencing_entry = pmo["sequencing_info"][0]
    assert sequencing_entry["seq_center"] == "UCSF Core"
    assert sequencing_entry["library_kit"] == "MiSeq Reagent Kit v3"

    library_entry = pmo["library_sample_info"][0]
    assert library_entry["accession"] == "ACC123"
    assert library_entry["library_note"] == "High quality"

    genome_entry = pmo["targeted_genomes"][0]
    assert genome_entry["chromosomes"] == ["chr1"]
    assert genome_entry["gff_url"] == "https://example.org/genome.gff"

    target_entry = pmo["target_info"][0]
    assert target_entry["gene_name"] == "geneA"
    assert target_entry["amplicon_length"] == 250
    assert target_entry["insert_location"]["chrom"] == "chr1"

    representative_mhap = pmo["representative_microhaplotypes"]["targets"][0][
        "microhaplotypes"
    ][0]
    assert representative_mhap["microhaplotype_name"] == "mh1"
    assert representative_mhap["masking"][0]["seq_segment_size"] == 2
    assert representative_mhap["custom_annotation"] == "custom"

    detected_mhap = pmo["detected_microhaplotypes"][0]["library_samples"][0][
        "target_results"
    ][0]["mhaps"][0]
    assert detected_mhap["umis"] == 10
    assert detected_mhap["custom_detected"] == "det-note"

    read_counts_run = pmo["read_counts_by_stage"][0]
    assert read_counts_run["bioinformatics_run_id"] == 0
    library_counts = read_counts_run["read_counts_by_library_sample_by_stage"][0]
    assert library_counts["library_sample_id"] == 0
    target_counts = library_counts["read_counts_for_targets"][0]
    assert target_counts["target_id"] == 0
    stage_entry = target_counts["stages"][0]
    assert stage_entry["quality_score"] == 0.98
    assert stage_entry["coverage_depth"] == 150

    assert pmo["project_info"][0]["project_type"] == "Surveillance"
