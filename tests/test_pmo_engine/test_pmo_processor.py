#!/usr/bin/env python3

import os
import pickle
import unittest
import json

import pandas as pd

from pmotools.pmo_engine.pmo_checker import PMOChecker
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

    def test_count_specimen_per_meta_fields(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        specimen_meta_fields_counts = PMOProcessor.count_specimen_per_meta_fields(pmo_data)
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

    def test_count_specimen_by_field_value(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        specimen_meta_subfields_counts = PMOProcessor.count_specimen_by_field_value(pmo_data, ["collection_country"])

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
        allele_data = PMOProcessor.extract_alleles_per_sample_table(pmo_data).sort_values(by=['bioinformatics_run_id', 'experiment_sample_name', 'target_name', 'mhap_id'])
        allele_data.to_csv("extracted_alleles_per_sample_table_no_extra_args.csv", index = False)
        self.assertEqual("5e4d462b6de3ace0db54b8d7d2ec7a5e", md5sum_of_fnp("extracted_alleles_per_sample_table_no_extra_args.csv"))
#experiment_sample_name", "target_name", "mhap_id
        allele_data_with_seq_reads = PMOProcessor.extract_alleles_per_sample_table(pmo_data,
                                                                                   additional_microhap_fields = ["reads"],
                                                                                   additional_representative_info_fields = ["seq"]).sort_values(by=['bioinformatics_run_id', 'experiment_sample_name', 'target_name', 'mhap_id'])
        allele_data_with_seq_reads.to_csv("extracted_alleles_per_sample_table_no_extra_args_with_seq_reads.csv", index = False)
        self.assertEqual("72deae7afa21d0df5e7c0e83a18217b0", md5sum_of_fnp("extracted_alleles_per_sample_table_no_extra_args_with_seq_reads.csv"))

        allele_data_with_seq_reads_panel_id_collection_country = PMOProcessor.extract_alleles_per_sample_table(pmo_data,
                                                                                   additional_microhap_fields = ["reads"],
                                                                                   additional_representative_info_fields = ["seq"],
                                                                                   additional_experiment_info_fields = ["panel_id"],
                                                                                   additional_specimen_info_fields = ["collection_country"],
                                                                                   ).sort_values(by=['bioinformatics_run_id', 'experiment_sample_name', 'target_name', 'mhap_id'])
        allele_data_with_seq_reads_panel_id_collection_country.to_csv("extracted_alleles_per_sample_table_no_extra_args_with_seq_reads_panel_id_collection_country.csv", index = False)
        self.assertEqual("a8ccb14976e6bd5d320e6d912ab7fcd2", md5sum_of_fnp("extracted_alleles_per_sample_table_no_extra_args_with_seq_reads_panel_id_collection_country.csv"))

    def test_extract_from_pmo_with_read_filter(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        pmo_data_filtered = PMOProcessor.extract_from_pmo_with_read_filter(pmo_data, 1000)
        with open("read_filter_100_combined_pmo_example.json", "w") as f:
            json.dump(pmo_data_filtered, f)
        self.assertEqual("a624dbce6bf501e4fe970e139c31e575",md5sum_of_fnp("read_filter_100_combined_pmo_example.json"))
        # check pmo extracted against PMO schema
        pmo_jsonschema_fnp = os.path.join(os.path.dirname(os.path.dirname(self.working_dir)),
                                               "etc/portable_microhaplotype_object.schema.json")
        with open(pmo_jsonschema_fnp) as f:
            pmo_jsonschema_data = json.load(f)
        checker = PMOChecker(pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_filtered)

    def test_filter_pmo_by_target_ids(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        pmo_data_select_targets = PMOProcessor.filter_pmo_by_target_ids(pmo_data, {1,10,11,55})
        with open("combined_pmo_target_ids_1_10_11_55.json", "w") as f:
            json.dump(pmo_data_select_targets, f)
        self.assertEqual("a09d7320b95f549b1ddd36384c1a9349",md5sum_of_fnp("combined_pmo_target_ids_1_10_11_55.json"))
        # check pmo extracted against PMO schema
        pmo_jsonschema_fnp = os.path.join(os.path.dirname(os.path.dirname(self.working_dir)),
                                               "etc/portable_microhaplotype_object.schema.json")
        with open(pmo_jsonschema_fnp) as f:
            pmo_jsonschema_data = json.load(f)
        checker = PMOChecker(pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_select_targets)

    def test_filter_pmo_by_target_names(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        pmo_data_select_targets = PMOProcessor.filter_pmo_by_target_names(pmo_data, {"t96", "t80", "t34", "t55"})
        with open("combined_pmo_target_ids_t96_t80_t34_t55.json", "w") as f:
            json.dump(pmo_data_select_targets, f)
        self.assertEqual("af2135e65a0bb0e2bc14da0c18fe21bb",md5sum_of_fnp("combined_pmo_target_ids_t96_t80_t34_t55.json"))
        # check pmo extracted against PMO schema
        pmo_jsonschema_fnp = os.path.join(os.path.dirname(os.path.dirname(self.working_dir)),
                                               "etc/portable_microhaplotype_object.schema.json")
        with open(pmo_jsonschema_fnp) as f:
            pmo_jsonschema_data = json.load(f)
        checker = PMOChecker(pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_select_targets)

    def test_filter_pmo_by_experiment_sample_ids(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        pmo_data_select_targets = PMOProcessor.filter_pmo_by_experiment_sample_ids(pmo_data, {1, 3})
        with open("combined_pmo_example_ids_1_3.json", "w") as f:
            json.dump(pmo_data_select_targets, f)
        self.assertEqual("1961bf5b7f41439191ecf12a296152aa",md5sum_of_fnp("combined_pmo_example_ids_1_3.json"))
        # check pmo extracted against PMO schema
        pmo_jsonschema_fnp = os.path.join(os.path.dirname(os.path.dirname(self.working_dir)),
                                               "etc/portable_microhaplotype_object.schema.json")
        with open(pmo_jsonschema_fnp) as f:
            pmo_jsonschema_data = json.load(f)
        checker = PMOChecker(pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_select_targets)

    def test_filter_pmo_by_experiment_sample_names(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        pmo_data_select_experiment_sample_names = PMOProcessor.filter_pmo_by_experiment_sample_names(pmo_data, {"8025874217", "XUC009"})
        with open("combined_pmo_example_ids_8025874217_XUC009.json", "w") as f:
            json.dump(pmo_data_select_experiment_sample_names, f)
        self.assertEqual("fe440e75e79f6b0348339c799647f9be",md5sum_of_fnp("combined_pmo_example_ids_8025874217_XUC009.json"))
        # check pmo extracted against PMO schema
        pmo_jsonschema_fnp = os.path.join(os.path.dirname(os.path.dirname(self.working_dir)),
                                               "etc/portable_microhaplotype_object.schema.json")
        with open(pmo_jsonschema_fnp) as f:
            pmo_jsonschema_data = json.load(f)
        checker = PMOChecker(pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_select_experiment_sample_names)

    def test_filter_pmo_by_specimen_ids(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        pmo_data_select_targets = PMOProcessor.filter_pmo_by_specimen_ids(pmo_data, {0,2})
        with open("combined_pmo_specimen_ids_0_2.json", "w") as f:
            json.dump(pmo_data_select_targets, f)
        self.assertEqual("0ecc93cfd4730b2fbcf1c815ade61fa9",md5sum_of_fnp("combined_pmo_specimen_ids_0_2.json"))
        # check pmo extracted against PMO schema
        pmo_jsonschema_fnp = os.path.join(os.path.dirname(os.path.dirname(self.working_dir)),
                                               "etc/portable_microhaplotype_object.schema.json")
        with open(pmo_jsonschema_fnp) as f:
            pmo_jsonschema_data = json.load(f)
        checker = PMOChecker(pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_select_targets)

    def test_filter_pmo_by_specimen_names(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        pmo_data_select_targets = PMOProcessor.filter_pmo_by_specimen_names(pmo_data, {"8025874217", "5tbx"})
        with open("combined_pmo_specimen_ids_8025874217_5tbx.json", "w") as f:
            json.dump(pmo_data_select_targets, f)
        self.assertEqual("0ecc93cfd4730b2fbcf1c815ade61fa9",md5sum_of_fnp("combined_pmo_specimen_ids_8025874217_5tbx.json"))
        # check pmo extracted against PMO schema
        pmo_jsonschema_fnp = os.path.join(os.path.dirname(os.path.dirname(self.working_dir)),
                                               "etc/portable_microhaplotype_object.schema.json")
        with open(pmo_jsonschema_fnp) as f:
            pmo_jsonschema_data = json.load(f)
        checker = PMOChecker(pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_select_targets)

    def test_extract_from_pmo_samples_with_meta_groupings(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            pmo_data = json.load(f)
        pmo_data_select_meta, group_counts = PMOProcessor.extract_from_pmo_samples_with_meta_groupings(pmo_data, "collection_country=Mozambique")
        with open("combined_pmo_collection_country_Mozambique.json", "w") as f:
            json.dump(pmo_data_select_meta, f)
        self.assertEqual("7310416dd5354647bb6caa2572461c3b",md5sum_of_fnp("combined_pmo_collection_country_Mozambique.json"))
        # check pmo extracted against PMO schema
        pmo_jsonschema_fnp = os.path.join(os.path.dirname(os.path.dirname(self.working_dir)),
                                               "etc/portable_microhaplotype_object.schema.json")
        with open(pmo_jsonschema_fnp) as f:
            pmo_jsonschema_data = json.load(f)
        checker = PMOChecker(pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_select_meta)


if __name__ == "__main__":
    unittest.main()
