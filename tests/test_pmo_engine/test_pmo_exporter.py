#!/usr/bin/env python3
import copy
import gzip
import hashlib
import os
import tempfile
import unittest
import json
import pandas as pd
from pmotools.pmo_engine.pmo_exporter import PMOExporter


def md5sum_of_fnp(filename):
    with open(filename, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


class TestPMOExporter(unittest.TestCase):
    def setUp(self):
        self.working_dir = os.path.dirname(os.path.abspath(__file__))
        self.test_dir = tempfile.TemporaryDirectory()
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/combined_pmo_example.json"
            )
        ) as f:
            self.combined_pmo_data = json.load(f)
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/minimum_pmo_example.json"
            )
        ) as f:
            self.small_example_pmo_data = json.load(f)
        with gzip.open(
            os.path.join(
                os.path.dirname(self.working_dir),
                "data/minimum_fields_pmo_example1.json.gz",
            ),
            "rt",
        ) as f:
            self.minimum_fields_v1_1_0_pmo_data = json.load(f)
        # a real, schema-valid v1.1.0 PMO with only the required top-level sections
        with gzip.open(
            os.path.join(
                os.path.dirname(self.working_dir),
                "data/minimum_Furstenau2025_PMO.json.gz",
            ),
            "rt",
        ) as f:
            self.minimum_v1_1_0_pmo_data = json.load(f)

    def tearDown(self):
        self.test_dir.cleanup()

    def test_list_library_sample_names_per_specimen_name(self):
        id_counts = PMOExporter.list_library_sample_names_per_specimen_name(
            self.small_example_pmo_data
        )
        id_counts_check_data = {
            "specimen_name": ["8025874217", "8025874266"],
            "library_sample_name": ["8025874217_lib_name", "8025874266_lib_name"],
            "library_sample_count": [1, 1],
        }
        id_counts_check_df = pd.DataFrame(id_counts_check_data)
        pd.testing.assert_frame_equal(id_counts, id_counts_check_df)

        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json"
            )
        ) as f:
            pmo_data_2 = json.load(f)
        id_counts_2 = PMOExporter.list_library_sample_names_per_specimen_name(
            pmo_data_2
        )
        id_counts_check_data_2 = {
            "specimen_name": ["5tbx", "XUC009"],
            "library_sample_name": ["5tbx_lib_name", "XUC009_lib_name"],
            "library_sample_count": [1, 1],
        }
        id_counts_check_df_2 = pd.DataFrame(id_counts_check_data_2)
        pd.testing.assert_frame_equal(id_counts_check_df_2, id_counts_2)

    def test_extract_targets_insert_bed_loc(self):
        all_target_inserts = PMOExporter.extract_targets_insert_bed_loc(
            self.combined_pmo_data, sort_output=True
        )
        output1_fnp = os.path.join(self.test_dir.name, "all_target_inserts_test1.bed")
        PMOExporter.write_bed_locs(all_target_inserts, output1_fnp)
        self.assertEqual("b7e477fe327ad7ae85f78ddaa66c313c", md5sum_of_fnp(output1_fnp))

        all_target_inserts_8_6_10 = PMOExporter.extract_targets_insert_bed_loc(
            self.combined_pmo_data, [8, 6, 10, 20], sort_output=False
        )
        output2_fnp = os.path.join(
            self.test_dir.name, "all_target_inserts_test2_8_6_10.bed"
        )
        PMOExporter.write_bed_locs(all_target_inserts_8_6_10, output2_fnp)
        self.assertEqual("9be1da64a4794489e08eca11c240d879", md5sum_of_fnp(output2_fnp))

        all_target_inserts_8_6_10_sorted = PMOExporter.extract_targets_insert_bed_loc(
            self.combined_pmo_data, [8, 6, 10, 20], sort_output=True
        )
        output3_fnp = os.path.join(
            self.test_dir.name, "all_target_inserts_test2_8_6_10_sorted.bed"
        )
        PMOExporter.write_bed_locs(all_target_inserts_8_6_10_sorted, output3_fnp)
        self.assertEqual("1831cc6a6f9f2bc4036d1dfad90771a1", md5sum_of_fnp(output3_fnp))

    def test_extract_panels_insert_bed_loc(self):
        all_target_inserts = PMOExporter.extract_panels_insert_bed_loc(
            self.combined_pmo_data, sort_output=True
        )
        output_fnp = os.path.join(self.test_dir.name, "all_panel_inserts_test1.bed")
        PMOExporter.write_bed_locs(all_target_inserts, output_fnp)
        self.assertEqual("52b1f79a3a89f8265573fa54b5a7ce57", md5sum_of_fnp(output_fnp))

    def test_extract_panels_insert_bed_loc_covers_all_reactions(self):
        # regression: previously the function returned early inside the reaction
        # loop, so only the first reaction of the first panel was returned.
        import copy

        pmo = copy.deepcopy(self.combined_pmo_data)
        panel = pmo["panel_info"][0]
        targets = panel["reactions"][0]["panel_targets"]
        half = len(targets) // 2
        panel["reactions"] = [
            {"reaction_name": "pool1", "panel_targets": targets[:half]},
            {"reaction_name": "pool2", "panel_targets": targets[half:]},
        ]
        bed_locs = PMOExporter.extract_panels_insert_bed_loc(pmo, sort_output=False)
        # every target across BOTH reactions is present
        self.assertEqual(len(bed_locs), len(targets))
        self.assertEqual(
            {b.name for b in bed_locs},
            {pmo["target_info"][t]["target_name"] for t in targets},
        )

    def test_write_bed_locs_header_on_own_line(self):
        # regression: header was written without a trailing newline, gluing the
        # first data row onto it.
        bed_locs = PMOExporter.extract_targets_insert_bed_loc(
            self.combined_pmo_data, sort_output=True
        )
        out_fnp = os.path.join(self.test_dir.name, "with_header.bed")
        PMOExporter.write_bed_locs(bed_locs, out_fnp, add_header=True)
        with open(out_fnp) as f:
            lines = f.read().splitlines()
        self.assertEqual(
            lines[0],
            "\t".join(
                [
                    "#chrom",
                    "start",
                    "end",
                    "name",
                    "score",
                    "strand",
                    "ref_seq",
                    "extra_info",
                ]
            ),
        )
        self.assertEqual(len(lines), len(bed_locs) + 1)

    def test_extract_alleles_per_sample_table(self):
        allele_data = PMOExporter.extract_alleles_per_sample_table(
            self.combined_pmo_data,
            additional_microhap_fields=["mhap_id"],
        ).sort_values(
            by=[
                "bioinformatics_run_name",
                "library_sample_name",
                "target_name",
                "mhap_id",
            ]
        )
        output_fnp = os.path.join(
            self.test_dir.name, "extracted_alleles_per_sample_table_no_extra_args.csv"
        )
        allele_data.to_csv(output_fnp, index=False)
        self.assertEqual("0c9df242c0f990682b32e8211bfa198c", md5sum_of_fnp(output_fnp))

        allele_data_with_seq_reads = PMOExporter.extract_alleles_per_sample_table(
            self.combined_pmo_data, additional_microhap_fields=["reads", "mhap_id"]
        ).sort_values(
            by=[
                "bioinformatics_run_name",
                "library_sample_name",
                "target_name",
                "mhap_id",
            ]
        )
        output_fnp = os.path.join(
            self.test_dir.name,
            "extracted_alleles_per_sample_table_no_extra_args_with_seq_reads.csv",
        )
        allele_data_with_seq_reads.to_csv(output_fnp, index=False)
        self.assertEqual("742c2c40546bdb25d1e2d517174120bb", md5sum_of_fnp(output_fnp))

        allele_data_with_seq_reads_panel_id_collection_country = (
            PMOExporter.extract_alleles_per_sample_table(
                self.combined_pmo_data,
                additional_microhap_fields=["reads", "mhap_id"],
                additional_library_sample_info_fields=["panel_id"],
                additional_specimen_info_fields=["collection_country"],
            ).sort_values(
                by=[
                    "bioinformatics_run_name",
                    "library_sample_name",
                    "target_name",
                    "mhap_id",
                ]
            )
        )
        output_fnp = os.path.join(
            self.test_dir.name,
            "extracted_alleles_per_sample_table_no_extra_args_with_seq_reads_panel_id_collection_country.csv",
        )
        allele_data_with_seq_reads_panel_id_collection_country.to_csv(
            output_fnp, index=False
        )
        self.assertEqual("189f1c73418c3cb85fcec2a736ff23b9", md5sum_of_fnp(output_fnp))

    def test_extract_alleles_per_sample_table_minimum_fields_pmo_input(self):
        allele_data = PMOExporter.extract_alleles_per_sample_table(
            self.minimum_fields_v1_1_0_pmo_data,
            additional_microhap_fields=["mhap_id"],
        ).sort_values(
            by=[
                "bioinformatics_run_name",
                "library_sample_name",
                "target_name",
                "mhap_id",
            ]
        )
        output_fnp = os.path.join(
            self.test_dir.name,
            "extracted_alleles_per_sample_table_no_extra_args_on_minimum_fields_pmo.csv",
        )
        allele_data.to_csv(output_fnp, index=False)
        self.assertEqual("cd016ab8d619328b32f11506b908e05f", md5sum_of_fnp(output_fnp))

        allele_data_with_seq_reads = PMOExporter.extract_alleles_per_sample_table(
            self.minimum_fields_v1_1_0_pmo_data,
            additional_microhap_fields=["reads", "mhap_id"],
        ).sort_values(
            by=[
                "bioinformatics_run_name",
                "library_sample_name",
                "target_name",
                "mhap_id",
            ]
        )
        output_fnp = os.path.join(
            self.test_dir.name,
            "extracted_alleles_per_sample_table_no_extra_args_with_seq_reads_on_minimum_fields_pmo.csv",
        )
        allele_data_with_seq_reads.to_csv(output_fnp, index=False)
        self.assertEqual("2810b465c005d4c1acbed12856ddfd26", md5sum_of_fnp(output_fnp))

    def test_export_specimen_meta_table(self):
        spec_table = PMOExporter.export_specimen_meta_table(self.small_example_pmo_data)
        spec_table.to_csv(os.path.join(self.test_dir.name, "specimen_meta_table.csv"))
        self.assertEqual(
            "8f94b8b774696e26c4ff6c8086e616a4",
            md5sum_of_fnp(os.path.join(self.test_dir.name, "specimen_meta_table.csv")),
        )

    def test_export_target_info_meta_table(self):
        target_info_table = PMOExporter.export_target_info_meta_table(
            self.small_example_pmo_data
        )
        target_info_table.to_csv(
            os.path.join(self.test_dir.name, "target_info_table.csv")
        )
        self.assertEqual(
            "2397407dcff8be3fdf54d27ba9a9cbff",
            md5sum_of_fnp(os.path.join(self.test_dir.name, "target_info_table.csv")),
        )

    def test_export_panel_info_meta_table(self):
        panel_info_table = PMOExporter.export_panel_info_meta_table(
            self.small_example_pmo_data
        )
        panel_info_table.to_csv(
            os.path.join(self.test_dir.name, "panel_info_table.csv")
        )
        self.assertEqual(
            "e5127ecaf7fe7950395d6f3d45f1c82a",
            md5sum_of_fnp(os.path.join(self.test_dir.name, "panel_info_table.csv")),
        )

    def test_export_library_sample_meta_table(self):
        library_sample_table = PMOExporter.export_library_sample_meta_table(
            self.small_example_pmo_data
        )
        library_sample_table.to_csv(
            os.path.join(self.test_dir.name, "library_sample_table.csv")
        )
        self.assertEqual(
            "7c433a74d215708e9339b5f6dece0bf3",
            md5sum_of_fnp(os.path.join(self.test_dir.name, "library_sample_table.csv")),
        )

    def test_export_meta_tables_without_optional_cross_ref_sections(self):
        # project_info and sequencing_info are optional as of v1.1.0; the specimen and
        # library_sample meta-table exporters must not raise a raw KeyError when a
        # referencing id is present but the referenced optional section is absent
        self.assertNotIn("project_info", self.minimum_v1_1_0_pmo_data)
        self.assertNotIn("sequencing_info", self.minimum_v1_1_0_pmo_data)
        # a clean minimal PMO (no dangling ids) should export fine
        PMOExporter.export_specimen_meta_table(self.minimum_v1_1_0_pmo_data)
        PMOExporter.export_library_sample_meta_table(self.minimum_v1_1_0_pmo_data)

        # inject the referencing ids without the optional sections (referential dangling)
        # the raw id is kept rather than resolving an unavailable name
        dangling = copy.deepcopy(self.minimum_v1_1_0_pmo_data)
        dangling["specimen_info"][0]["project_id"] = 0
        dangling["library_sample_info"][0]["sequencing_info_id"] = 0
        spec_table = PMOExporter.export_specimen_meta_table(dangling)
        self.assertIn("project_id", spec_table.columns)
        self.assertNotIn("project_name", spec_table.columns)
        library_table = PMOExporter.export_library_sample_meta_table(dangling)
        self.assertIn("sequencing_info_id", library_table.columns)
        self.assertNotIn("sequencing_info_name", library_table.columns)

    def test_export_sequencing_info_meta_table(self):
        sequencing_info_table = PMOExporter.export_sequencing_info_meta_table(
            self.small_example_pmo_data
        )
        sequencing_info_table.to_csv(
            os.path.join(self.test_dir.name, "sequencing_info_table.csv")
        )
        self.assertEqual(
            "1cc6fb83227752454cfc3ba63eac503b",
            md5sum_of_fnp(
                os.path.join(self.test_dir.name, "sequencing_info_table.csv")
            ),
        )

    def test_export_project_info_meta_table(self):
        project_info_table = PMOExporter.export_project_info_meta_table(
            self.small_example_pmo_data
        )
        project_info_table.to_csv(
            os.path.join(self.test_dir.name, "project_info_table.csv")
        )
        self.assertEqual(
            "e533098411cbd96de2733668e8475ab8",
            md5sum_of_fnp(os.path.join(self.test_dir.name, "project_info_table.csv")),
        )

    def test_export_specimen_travel_meta_table(self):
        test_pmo_with_travel_info = {
            "specimen_info": [
                {
                    "specimen_name": "spec1",
                    "travel_out_six_month": [
                        {
                            "travel_country": "Kenya",
                            "travel_start_date": "2024-01",
                            "travel_end_date": "2024-02",
                        },
                        {
                            "travel_country": "Kenya",
                            "travel_start_date": "2024-04",
                            "travel_end_date": "2024-06",
                        },
                    ],
                },
                {
                    "specimen_name": "spec2",
                    "travel_out_six_month": [
                        {
                            "travel_country": "Tanzania",
                            "travel_start_date": "2024-02-15",
                            "travel_end_date": "2024-02-27",
                        }
                    ],
                },
            ]
        }
        specimen_trable_info_table = PMOExporter.export_specimen_travel_meta_table(
            test_pmo_with_travel_info
        )
        specimen_trable_info_table.to_csv(
            os.path.join(self.test_dir.name, "specimen_trable_info_table.csv")
        )
        self.assertEqual(
            "0305350d655184aa385d3d1ddc9b3600",
            md5sum_of_fnp(
                os.path.join(self.test_dir.name, "specimen_trable_info_table.csv")
            ),
        )

    def test_basic_structure_from_minimum_example(self):
        """DataFrame has expected columns and one row for the single genome."""
        df = PMOExporter.export_targeted_genomes_meta_table(self.small_example_pmo_data)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1
        assert "genome_id" in df.columns
        assert "name" in df.columns

    def test_genome_id_assigned(self):
        """genome_id starts at 0 and increments per genome."""
        df = PMOExporter.export_targeted_genomes_meta_table(self.small_example_pmo_data)
        assert df["genome_id"].iloc[0] == 0

    def test_raises_on_missing_targeted_genomes(self):
        """Raises ValueError when targeted_genomes key is absent."""
        with self.assertRaises(ValueError) as context:
            PMOExporter.export_targeted_genomes_meta_table({"pmo_header": {}})
        self.assertIn(
            "no targeted_genomes found",
            str(context.exception),
        )

    def test_multiple_genomes_from_combined_example(self):
        """One row is produced per genome entry."""
        df = PMOExporter.export_targeted_genomes_meta_table(self.combined_pmo_data)
        assert len(df) == len(self.combined_pmo_data["targeted_genomes"])
        assert list(df["genome_id"]) == list(range(len(df)))

    def test_empty_targeted_genomes_list(self):
        """Empty targeted_genomes list produces an empty DataFrame."""
        pmodata = {"targeted_genomes": []}
        df = PMOExporter.export_targeted_genomes_meta_table(pmodata)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0

    def test_basic_structure_from_minimum_example_bioinformatics_run_info(self):
        """DataFrame has expected columns and correct row count."""
        df = PMOExporter.export_bioinformatics_run_info_meta_table(
            self.small_example_pmo_data
        )
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1

    def test_raises_on_missing_bioinformatics_run_info(self):
        """Raises ValueError when bioinformatics_run_info key is absent."""
        with self.assertRaises(ValueError) as context:
            PMOExporter.export_bioinformatics_run_info_meta_table({"pmo_header": {}})
        self.assertIn(
            "no bioinformatics_run_info found",
            str(context.exception),
        )

    def test_empty_bioinformatics_run_info_list(self):
        """Empty bioinformatics_run_info list produces an empty DataFrame."""
        pmodata = {"bioinformatics_run_info": []}
        df = PMOExporter.export_bioinformatics_run_info_meta_table(pmodata)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0

    def test_multiple_runs_from_combined_example_bioinformatics_run_info(self):
        """One row is produced per bioinformatics_run_info entry."""
        df = PMOExporter.export_bioinformatics_run_info_meta_table(
            self.combined_pmo_data
        )
        assert len(df) == len(self.combined_pmo_data["bioinformatics_run_info"])
        assert list(df["run_id"]) == list(range(len(df)))

    def test_custom_separator_bioinformatics_run_info(self):
        """Custom separator is used for any list fields."""
        pmodata = {
            "bioinformatics_run_info": [
                {
                    "bioinformatics_run_name": "test-run",
                    "run_date": "2024-01-01",
                    "some_list_field": ["a", "b", "c"],
                }
            ]
        }
        df_pipe = PMOExporter.export_bioinformatics_run_info_meta_table(
            pmodata, separator="|"
        )
        assert df_pipe["some_list_field"].iloc[0] == "a|b|c"
        df_comma = PMOExporter.export_bioinformatics_run_info_meta_table(
            pmodata, separator=","
        )
        assert df_comma["some_list_field"].iloc[0] == "a,b,c"

    def test_basic_structure_from_minimum_example_bioinformatics_methods_info(self):
        """DataFrame is a pandas DataFrame with expected columns."""
        df = PMOExporter.export_bioinformatics_methods_info_meta_table(
            self.small_example_pmo_data
        )
        assert isinstance(df, pd.DataFrame)
        assert "bioinformatics_methods_id" in df.columns
        assert "method_id" in df.columns

    def test_row_count_matches_total_methods_bioinformatics_methods_info(self):
        """Total rows equals the sum of methods across all bioinformatics_methods_info entries."""
        df = PMOExporter.export_bioinformatics_methods_info_meta_table(
            self.small_example_pmo_data
        )
        expected_row_count = sum(
            len(entry["methods"])
            for entry in self.small_example_pmo_data["bioinformatics_methods_info"]
        )
        assert len(df) == expected_row_count

    def test_row_count_matches_total_methods_combined_example_bioinformatics_methods_info(
        self,
    ):
        """Row count matches sum of methods in combined example."""
        df = PMOExporter.export_bioinformatics_methods_info_meta_table(
            self.combined_pmo_data
        )
        expected_row_count = sum(
            len(entry["methods"])
            for entry in self.combined_pmo_data["bioinformatics_methods_info"]
        )
        assert len(df) == expected_row_count

    def test_raises_on_missing_bioinformatics_methods_info(self):
        """Raises ValueError when bioinformatics_methods_info key is absent."""
        with self.assertRaises(ValueError) as context:
            PMOExporter.export_bioinformatics_methods_info_meta_table(
                {"pmo_header": {}}
            )
        self.assertIn(
            "no bioinformatics_methods_info found",
            str(context.exception),
        )

    def test_empty_bioinformatics_methods_info_list(self):
        """Empty bioinformatics_methods_info list produces an empty DataFrame."""
        pmodata = {"bioinformatics_methods_info": []}
        df = PMOExporter.export_bioinformatics_methods_info_meta_table(pmodata)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0

    def test_build_pmo_sheet_configs_required_sheets_always_present(self):
        """Required sheets are always included regardless of optional fields."""
        configs = PMOExporter._build_pmo_sheet_configs(self.small_example_pmo_data)
        sheet_names = [c.sheet_name for c in configs]
        required_sheets = [
            "PMO Header",
            "Required Panel Targets",
            "Required Panel Info",
            "Required Microhaplotype",
        ]
        for sheet in required_sheets:
            self.assertIn(sheet, sheet_names)

    def test_build_pmo_sheet_configs_optional_sheets_excluded_when_absent(self):
        """Optional sheets are only included when their key exists in the PMO."""
        configs = PMOExporter._build_pmo_sheet_configs(self.small_example_pmo_data)
        sheet_names = [c.sheet_name for c in configs]
        optional_key_sheet_pairs = [
            ("targeted_genomes", "Optional GenomeInfo"),
            ("project_info", "Optional ProjectInfo"),
            ("sequencing_info", "Optional SequencingInfo"),
            ("bioinformatics_methods_info", "Optional Bioinformatics Methods"),
            ("bioinformatics_run_info", "Optional Bioinformatics Run"),
        ]
        for pmo_key, sheet_name in optional_key_sheet_pairs:
            if pmo_key not in self.small_example_pmo_data:
                self.assertNotIn(sheet_name, sheet_names)
            else:
                self.assertIn(sheet_name, sheet_names)

    def test_build_pmo_sheet_configs_all_optional_sheets_present_in_combined(self):
        """All optional sheets appear when the combined (fully-populated) PMO is used."""
        configs = PMOExporter._build_pmo_sheet_configs(self.combined_pmo_data)
        sheet_names = [c.sheet_name for c in configs]
        optional_sheets = [
            "Optional GenomeInfo",
            "Optional ProjectInfo",
            "Optional SequencingInfo",
            "Optional Bioinformatics Methods",
            "Optional Bioinformatics Run",
        ]
        for sheet in optional_sheets:
            self.assertIn(sheet, sheet_names)

    def test_export_to_excel_creates_valid_file_with_expected_sheets(self):
        """export_to_excel writes a valid xlsx file containing all expected sheet names."""
        output_fnp = os.path.join(self.test_dir.name, "test_export.xlsx")
        PMOExporter.export_to_excel(self.small_example_pmo_data, output_fnp)
        self.assertTrue(os.path.exists(output_fnp))
        written_sheets = pd.ExcelFile(output_fnp).sheet_names
        expected_sheets = [
            "PMO Header",
            "Required Panel Targets",
            "Required Panel Info",
            "Required Microhaplotype",
            "Optional Specimen Level",
            "Optional LibrarySampleInfo",
        ]
        for sheet in expected_sheets:
            self.assertIn(sheet, written_sheets)


if __name__ == "__main__":
    unittest.main()
