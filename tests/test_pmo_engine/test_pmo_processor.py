#!/usr/bin/env python3

import os
import pickle
import unittest
import json

import pandas as pd

from pmotools.pmo_engine.pmo_processor import PMOProcessor
import hashlib
import gzip

def md5sum_of_fnp(filename):
    with open(filename, 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()


class TestPMOProcessor(unittest.TestCase):
    def setUp(self):
        self.working_dir = os.path.dirname(os.path.abspath(__file__))

    def test_list_experiment_sample_ids_per_specimen_id(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example.json")) as f:
            pmo_data = json.load(f)
        id_counts = PMOProcessor.list_experiment_sample_ids_per_specimen_id(pmo_data)
        id_counts_check_data = {
            "specimen_name": ["8025874217", "8025874266"],
            "experiment_sample_name": ["8025874217", "8025874266"],
            "experiment_sample_count": [1, 1]
        }
        id_counts_check_df = pd.DataFrame(id_counts_check_data)
        pd.testing.assert_frame_equal(id_counts_check_df,id_counts_check_df)

        with open(os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json")) as f:
            pmo_data_2 = json.load(f)
        id_counts_2 = PMOProcessor.list_experiment_sample_ids_per_specimen_id(pmo_data_2)
        id_counts_check_data_2 = {
            "specimen_name": ["5tbx", "XUC009"],
            "experiment_sample_name": ["5tbx", "XUC009"],
            "experiment_sample_count": [1, 1]
        }
        id_counts_check_df_2 = pd.DataFrame(id_counts_check_data_2)
        pd.testing.assert_frame_equal(id_counts_check_df_2,id_counts_2)

    def test_count_targets_per_sample(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example.json")) as f:
            pmo_data = json.load(f)
        targets_per_sample_counts = PMOProcessor.count_targets_per_sample(pmo_data)
        targets_per_sample_check_data = {
            "bioinformatics_run_id": [0, 0],
            "experiment_sample_name": ["8025874266", "8025874217"],
            "target_number": [85, 99]
        }
        targets_per_sample_check_df = pd.DataFrame(targets_per_sample_check_data)
        pd.testing.assert_frame_equal(targets_per_sample_check_df,targets_per_sample_counts)

        targets_per_sample_counts_read_count_off1000 = PMOProcessor.count_targets_per_sample(pmo_data, 1000)
        targets_per_sample_read_count_off1000_check_data = {
            "bioinformatics_run_id": [0, 0],
            "experiment_sample_name": ["8025874266", "8025874217"],
            "target_number": [61, 99]
        }
        targets_per_sample_read_count_off1000_check_df = pd.DataFrame(targets_per_sample_read_count_off1000_check_data)
        pd.testing.assert_frame_equal(targets_per_sample_counts_read_count_off1000,targets_per_sample_read_count_off1000_check_df)

        with open(os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json")) as f:
            pmo_data_2 = json.load(f)
        targets_per_sample_counts_2 = PMOProcessor.count_targets_per_sample(pmo_data_2)
        targets_per_sample_check_data_2 = {
            "bioinformatics_run_id": [0, 0],
            "experiment_sample_name": ["XUC009", "5tbx"],
            "target_number": [98, 100]
        }
        targets_per_sample_check_df_2 = pd.DataFrame(targets_per_sample_check_data_2)
        pd.testing.assert_frame_equal(targets_per_sample_check_df_2,targets_per_sample_check_df_2)

        targets_per_sample_counts_2_read_count_off200 = PMOProcessor.count_targets_per_sample(pmo_data_2, 200)
        targets_per_sample_check_data_2_read_count_off200 = {
            "bioinformatics_run_id": [0, 0],
            "experiment_sample_name": ["XUC009", "5tbx"],
            "target_number": [96, 100]
        }
        targets_per_sample_check_df_2_read_count_off200 = pd.DataFrame(targets_per_sample_check_data_2_read_count_off200)
        pd.testing.assert_frame_equal(targets_per_sample_check_df_2_read_count_off200,targets_per_sample_counts_2_read_count_off200)

        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data_combined = json.load(f)
        pmo_data_combined_targets_per_sample = PMOProcessor.count_targets_per_sample(pmo_data_combined)
        pmo_data_combined_targets_per_sample_check = pd.DataFrame([
            {"bioinformatics_run_id": 0, "experiment_sample_name": "8025874266", "target_number": 85},
            {"bioinformatics_run_id": 0, "experiment_sample_name": "8025874217", "target_number": 99},
            {"bioinformatics_run_id": 1, "experiment_sample_name": "XUC009", "target_number": 98},
            {"bioinformatics_run_id": 1, "experiment_sample_name": "5tbx", "target_number": 100}
        ])
        pd.testing.assert_frame_equal(pmo_data_combined_targets_per_sample_check,pmo_data_combined_targets_per_sample)

    def test_count_samples_per_target(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)

        # counts per bio run
        targets_per_sample_counts = PMOProcessor.count_samples_per_target(pmo_data)
        targets_per_sample_counts_check = pd.read_csv(
            os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example_sample_per_target_counts.csv"),
            sep=','
        ).sort_values(by=['bioinformatics_run_id', 'target_name']).reset_index(drop=True)
        pd.testing.assert_frame_equal(targets_per_sample_counts_check, targets_per_sample_counts)

        # collapsing across bio runs
        targets_per_sample_counts_collapsed = PMOProcessor.count_samples_per_target(pmo_data, collapse_across_runs=True)
        targets_per_sample_counts_collapsed_check = pd.read_csv(
            os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example_sample_per_target_counts_collapsed.csv"),
            sep=','
        ).sort_values(by=['target_name']).reset_index(drop=True)
        pd.testing.assert_frame_equal(targets_per_sample_counts_collapsed_check, targets_per_sample_counts_collapsed)

    def test_count_specimen_meta_fields(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        specimen_meta_fields_counts = PMOProcessor.count_specimen_meta_fields(pmo_data)
        specimen_meta_fields_counts_check = pd.DataFrame([
            {"field": "collection_country", "present_in_specimens_count": 4, "total_specimen_count": 4},
            {"field": "collection_date", "present_in_specimens_count": 4, "total_specimen_count": 4},
            {"field": "collector_chief_scientist", "present_in_specimens_count": 4, "total_specimen_count": 4},
            {"field": "geo_admin3", "present_in_specimens_count": 4, "total_specimen_count": 4},
            {"field": "host_taxon_id", "present_in_specimens_count": 4, "total_specimen_count": 4},
            {"field": "lat_lon", "present_in_specimens_count": 2, "total_specimen_count": 4},
            {"field": "parasite_density_info", "present_in_specimens_count": 2, "total_specimen_count": 4},
            {"field": "plate_info", "present_in_specimens_count": 2, "total_specimen_count": 4},
            {"field": "project_name", "present_in_specimens_count": 4, "total_specimen_count": 4},
            {"field": "specimen_collect_device", "present_in_specimens_count": 4, "total_specimen_count": 4},
            {"field": "specimen_name", "present_in_specimens_count": 4, "total_specimen_count": 4},
            {"field": "specimen_store_loc", "present_in_specimens_count": 4, "total_specimen_count": 4},
            {"field": "specimen_taxon_id", "present_in_specimens_count": 4, "total_specimen_count": 4},
        ])
        pd.testing.assert_frame_equal(specimen_meta_fields_counts, specimen_meta_fields_counts_check)

    def test_count_specimen_meta_subfields(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        specimen_meta_subfields_counts = PMOProcessor.count_specimen_meta_subfields(pmo_data, ["collection_country"])

        specimen_meta_subfields_counts_check = pd.DataFrame([
            {"collection_country": "India", "specimens_count": 1, "specimens_freq": 0.25, "total_specimen_count": 4},
            {"collection_country": "Mozambique", "specimens_count": 2, "specimens_freq": 0.50,"total_specimen_count": 4},
            {"collection_country": "Papua New Guinea", "specimens_count": 1, "specimens_freq": 0.25,"total_specimen_count": 4},
        ])
        pd.testing.assert_frame_equal(specimen_meta_subfields_counts, specimen_meta_subfields_counts_check)

    def test_extract_allele_counts_freq_from_pmo(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        allele_counts = PMOProcessor.extract_allele_counts_freq_from_pmo(pmo_data)
        df_string = pd.util.hash_pandas_object(allele_counts).values
        self.assertEqual(hashlib.md5(df_string).hexdigest(), "79fadd12ff7033c23810afa6cd46b4a5")

        allele_counts_collapsed = PMOProcessor.extract_allele_counts_freq_from_pmo(pmo_data, collapse_across_runs=True)
        allele_counts_collapsed_string = pd.util.hash_pandas_object(allele_counts_collapsed).values
        self.assertEqual(hashlib.md5(allele_counts_collapsed_string).hexdigest(), "cb34c7e1357e2e35024a89464b63f06c")


    def test_extract_targets_insert_bed_loc(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        all_target_inserts = PMOProcessor.extract_targets_insert_bed_loc(pmo_data, sort_output=True)
        PMOProcessor.write_bed_locs(all_target_inserts, "all_target_inserts_test1.bed")
        self.assertEqual("b7e477fe327ad7ae85f78ddaa66c313c", md5sum_of_fnp("all_target_inserts_test1.bed"))

        all_target_inserts_8_6_10 = PMOProcessor.extract_targets_insert_bed_loc(pmo_data, [8,6,10,20], sort_output=False)
        PMOProcessor.write_bed_locs(all_target_inserts_8_6_10, "all_target_inserts_test2_8_6_10.bed")
        self.assertEqual("9be1da64a4794489e08eca11c240d879", md5sum_of_fnp("all_target_inserts_test2_8_6_10.bed"))

        all_target_inserts_8_6_10_sorted = PMOProcessor.extract_targets_insert_bed_loc(pmo_data, [8,6,10,20], sort_output=True)
        PMOProcessor.write_bed_locs(all_target_inserts_8_6_10_sorted, "all_target_inserts_test2_8_6_10_sorted.bed")
        self.assertEqual("1831cc6a6f9f2bc4036d1dfad90771a1", md5sum_of_fnp("all_target_inserts_test2_8_6_10_sorted.bed"))

    def test_extract_panels_insert_bed_loc(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        all_target_inserts = PMOProcessor.extract_panels_insert_bed_loc(pmo_data, sort_output=True)
        PMOProcessor.write_bed_locs(all_target_inserts, "all_panel_inserts_test1.bed")
        self.assertEqual("52b1f79a3a89f8265573fa54b5a7ce57", md5sum_of_fnp("all_panel_inserts_test1.bed"))

    def test_extract_alleles_per_sample_table(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        allele_data = PMOProcessor.extract_alleles_per_sample_table(pmo_data).sort_values(by=['bioinformatics_run_id', 'sampleID', 'locus', 'allele'])
        allele_data.to_csv("extracted_alleles_per_sample_table_no_extra_args.csv", index = False)
        self.assertEqual("9efc4b744faa77176c9e224dbb2bf809", md5sum_of_fnp("extracted_alleles_per_sample_table_no_extra_args.csv"))

        allele_data_with_seq_reads = PMOProcessor.extract_alleles_per_sample_table(pmo_data,
                                                                                   additional_microhap_fields = ["reads"],
                                                                                   additional_representative_infos_fields = ["seq"]).sort_values(by=['bioinformatics_run_id', 'sampleID', 'locus', 'allele'])
        allele_data_with_seq_reads.to_csv("extracted_alleles_per_sample_table_no_extra_args_with_seq_reads.csv", index = False)
        self.assertEqual("bd525642ec9ed0a29996304595b43b85", md5sum_of_fnp("extracted_alleles_per_sample_table_no_extra_args_with_seq_reads.csv"))

        allele_data_with_seq_reads_panel_id_collection_country = PMOProcessor.extract_alleles_per_sample_table(pmo_data,
                                                                                   additional_microhap_fields = ["reads"],
                                                                                   additional_representative_infos_fields = ["seq"],
                                                                                   additional_experiment_infos_fields = ["panel_id"],
                                                                                   additional_specimen_info_fields = ["collection_country"],
                                                                                   ).sort_values(by=['bioinformatics_run_id', 'sampleID', 'locus', 'allele'])
        allele_data_with_seq_reads_panel_id_collection_country.to_csv("extracted_alleles_per_sample_table_no_extra_args_with_seq_reads_panel_id_collection_country.csv", index = False)
        self.assertEqual("8409a859490b5b705f828a8579b5bdfc", md5sum_of_fnp("extracted_alleles_per_sample_table_no_extra_args_with_seq_reads_panel_id_collection_country.csv"))




if __name__ == "__main__":
    unittest.main()
