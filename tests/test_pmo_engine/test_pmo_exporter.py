#!/usr/bin/env python3
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
            self.minimum_pmo_data = json.load(f)

    def tearDown(self):
        self.test_dir.cleanup()

    def test_list_library_sample_names_per_specimen_name(self):
        id_counts = PMOExporter.list_library_sample_names_per_specimen_name(
            self.minimum_pmo_data
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

    def test_extract_alleles_per_sample_table(self):
        allele_data = PMOExporter.extract_alleles_per_sample_table(
            self.combined_pmo_data
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
        self.assertEqual("2898d87133e2e381612f3c0dea70122f", md5sum_of_fnp(output_fnp))

        allele_data_with_seq_reads = PMOExporter.extract_alleles_per_sample_table(
            self.combined_pmo_data,
            additional_microhap_fields=["reads"],
            additional_representative_info_fields=["seq"],
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
        self.assertEqual("744c1c0233066f030881c8b595b9ad5c", md5sum_of_fnp(output_fnp))

        allele_data_with_seq_reads_panel_id_collection_country = (
            PMOExporter.extract_alleles_per_sample_table(
                self.combined_pmo_data,
                additional_microhap_fields=["reads"],
                additional_representative_info_fields=["seq"],
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
        self.assertEqual("c425004244e6af1386b6e7776da76fed", md5sum_of_fnp(output_fnp))

    def test_export_specimen_meta_table(self):
        spec_table = PMOExporter.export_specimen_meta_table(self.minimum_pmo_data)
        spec_table.to_csv(os.path.join(self.test_dir.name, "specimen_meta_table.csv"))
        self.assertEqual(
            "8f94b8b774696e26c4ff6c8086e616a4",
            md5sum_of_fnp(os.path.join(self.test_dir.name, "specimen_meta_table.csv")),
        )

    def test_export_target_info_meta_table(self):
        target_info_table = PMOExporter.export_target_info_meta_table(
            self.minimum_pmo_data
        )
        target_info_table.to_csv(
            os.path.join(self.test_dir.name, "target_info_table.csv")
        )
        self.assertEqual(
            "cb0319482c9da5f9d8b22fba955ce1c8",
            md5sum_of_fnp(os.path.join(self.test_dir.name, "target_info_table.csv")),
        )

    def test_export_panel_info_meta_table(self):
        panel_info_table = PMOExporter.export_panel_info_meta_table(
            self.minimum_pmo_data
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
            self.minimum_pmo_data
        )
        library_sample_table.to_csv(
            os.path.join(self.test_dir.name, "library_sample_table.csv")
        )
        self.assertEqual(
            "7c433a74d215708e9339b5f6dece0bf3",
            md5sum_of_fnp(os.path.join(self.test_dir.name, "library_sample_table.csv")),
        )

    def test_export_sequencing_info_meta_table(self):
        sequencing_info_table = PMOExporter.export_sequencing_info_meta_table(
            self.minimum_pmo_data
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
            self.minimum_pmo_data
        )
        project_info_table.to_csv(
            os.path.join(self.test_dir.name, "project_info_table.csv")
        )
        print(project_info_table)
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
        print(specimen_trable_info_table)
        self.assertEqual(
            "0305350d655184aa385d3d1ddc9b3600",
            md5sum_of_fnp(
                os.path.join(self.test_dir.name, "specimen_trable_info_table.csv")
            ),
        )


if __name__ == "__main__":
    unittest.main()
