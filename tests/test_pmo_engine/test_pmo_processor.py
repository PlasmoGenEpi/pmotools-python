#!/usr/bin/env python3

import os
import tempfile
import unittest
import json

import pandas as pd

from pmotools.pmo_engine.pmo_checker import PMOChecker
from pmotools.pmo_engine.pmo_processor import PMOProcessor
import hashlib
from pmotools.utils.schema_loader import load_schema


def md5sum_of_fnp(filename):
    with open(filename, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


class TestPMOProcessor(unittest.TestCase):
    def setUp(self):
        self.working_dir = os.path.dirname(os.path.abspath(__file__))
        self.test_dir = tempfile.TemporaryDirectory()
        self.pmo_jsonschema_data = load_schema(
            "portable_microhaplotype_object_v1.0.0.schema.json"
        )
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

    def test_count_targets_per_library_sample(self):
        targets_per_sample_counts = PMOProcessor.count_targets_per_library_sample(
            self.minimum_pmo_data
        )
        targets_per_sample_check_data = {
            "bioinformatics_run_id": [0, 0],
            "library_sample_name": ["8025874266_lib_name", "8025874217_lib_name"],
            "target_number": [85, 99],
        }
        targets_per_sample_check_df = pd.DataFrame(targets_per_sample_check_data)
        pd.testing.assert_frame_equal(
            targets_per_sample_check_df, targets_per_sample_counts
        )

        targets_per_sample_counts_read_count_off1000 = (
            PMOProcessor.count_targets_per_library_sample(self.minimum_pmo_data, 1000)
        )
        targets_per_sample_read_count_off1000_check_data = {
            "bioinformatics_run_id": [0, 0],
            "library_sample_name": ["8025874266_lib_name", "8025874217_lib_name"],
            "target_number": [61, 99],
        }
        targets_per_sample_read_count_off1000_check_df = pd.DataFrame(
            targets_per_sample_read_count_off1000_check_data
        )
        pd.testing.assert_frame_equal(
            targets_per_sample_counts_read_count_off1000,
            targets_per_sample_read_count_off1000_check_df,
        )

        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json"
            )
        ) as f:
            pmo_data_2 = json.load(f)
        targets_per_sample_counts_2 = PMOProcessor.count_targets_per_library_sample(
            pmo_data_2
        )
        targets_per_sample_check_data_2 = {
            "bioinformatics_run_id": [0, 0],
            "library_sample_name": ["XUC009_lib_name", "5tbx_lib_name"],
            "target_number": [98, 100],
        }
        targets_per_sample_check_df_2 = pd.DataFrame(targets_per_sample_check_data_2)
        pd.testing.assert_frame_equal(
            targets_per_sample_counts_2, targets_per_sample_check_df_2
        )

        targets_per_sample_counts_2_read_count_off200 = (
            PMOProcessor.count_targets_per_library_sample(pmo_data_2, 200)
        )
        targets_per_sample_check_data_2_read_count_off200 = {
            "bioinformatics_run_id": [0, 0],
            "library_sample_name": ["XUC009_lib_name", "5tbx_lib_name"],
            "target_number": [96, 100],
        }
        targets_per_sample_check_df_2_read_count_off200 = pd.DataFrame(
            targets_per_sample_check_data_2_read_count_off200
        )
        pd.testing.assert_frame_equal(
            targets_per_sample_check_df_2_read_count_off200,
            targets_per_sample_counts_2_read_count_off200,
        )

        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/combined_pmo_example.json"
            )
        ) as f:
            pmo_data_combined = json.load(f)
        pmo_data_combined_targets_per_sample = (
            PMOProcessor.count_targets_per_library_sample(pmo_data_combined)
        )
        pmo_data_combined_targets_per_sample_check = pd.DataFrame(
            [
                {
                    "bioinformatics_run_id": 0,
                    "library_sample_name": "8025874266_lib_name",
                    "target_number": 85,
                },
                {
                    "bioinformatics_run_id": 0,
                    "library_sample_name": "8025874217_lib_name",
                    "target_number": 99,
                },
                {
                    "bioinformatics_run_id": 1,
                    "library_sample_name": "XUC009_lib_name",
                    "target_number": 98,
                },
                {
                    "bioinformatics_run_id": 1,
                    "library_sample_name": "5tbx_lib_name",
                    "target_number": 100,
                },
            ]
        )
        pd.testing.assert_frame_equal(
            pmo_data_combined_targets_per_sample_check,
            pmo_data_combined_targets_per_sample,
        )

    def test_count_library_samples_per_target(self):
        # counts per bio run
        targets_per_sample_counts = PMOProcessor.count_library_samples_per_target(
            self.combined_pmo_data
        )
        targets_per_sample_counts_check = (
            pd.read_csv(
                os.path.join(
                    os.path.dirname(self.working_dir),
                    "data/combined_pmo_example_sample_per_target_counts.csv",
                ),
                sep=",",
            )
            .sort_values(by=["bioinformatics_run_id", "target_name"])
            .reset_index(drop=True)
        )
        pd.testing.assert_frame_equal(
            targets_per_sample_counts_check, targets_per_sample_counts
        )

        # collapsing across bio runs
        targets_per_sample_counts_collapsed = (
            PMOProcessor.count_library_samples_per_target(
                self.combined_pmo_data, collapse_across_runs=True
            )
        )
        targets_per_sample_counts_collapsed_check = (
            pd.read_csv(
                os.path.join(
                    os.path.dirname(self.working_dir),
                    "data/combined_pmo_example_sample_per_target_counts_collapsed.csv",
                ),
                sep=",",
            )
            .sort_values(by=["target_name"])
            .reset_index(drop=True)
        )
        pd.testing.assert_frame_equal(
            targets_per_sample_counts_collapsed_check,
            targets_per_sample_counts_collapsed,
        )

    def test_count_specimen_per_meta_fields(self):
        specimen_meta_fields_counts = PMOProcessor.count_specimen_per_meta_fields(
            self.combined_pmo_data
        )
        specimen_meta_fields_counts_check = pd.DataFrame(
            [
                {
                    "field": "collection_country",
                    "present_in_specimens_count": 4,
                    "total_specimen_count": 4,
                },
                {
                    "field": "collection_date",
                    "present_in_specimens_count": 4,
                    "total_specimen_count": 4,
                },
                {
                    "field": "geo_admin3",
                    "present_in_specimens_count": 4,
                    "total_specimen_count": 4,
                },
                {
                    "field": "host_taxon_id",
                    "present_in_specimens_count": 4,
                    "total_specimen_count": 4,
                },
                {
                    "field": "lat_lon",
                    "present_in_specimens_count": 2,
                    "total_specimen_count": 4,
                },
                {
                    "field": "parasite_density_info",
                    "present_in_specimens_count": 2,
                    "total_specimen_count": 4,
                },
                {
                    "field": "project_id",
                    "present_in_specimens_count": 4,
                    "total_specimen_count": 4,
                },
                {
                    "field": "specimen_collect_device",
                    "present_in_specimens_count": 4,
                    "total_specimen_count": 4,
                },
                {
                    "field": "specimen_name",
                    "present_in_specimens_count": 4,
                    "total_specimen_count": 4,
                },
                {
                    "field": "specimen_store_loc",
                    "present_in_specimens_count": 4,
                    "total_specimen_count": 4,
                },
                {
                    "field": "specimen_taxon_id",
                    "present_in_specimens_count": 4,
                    "total_specimen_count": 4,
                },
                {
                    "field": "storage_plate_info",
                    "present_in_specimens_count": 2,
                    "total_specimen_count": 4,
                },
            ]
        )
        pd.testing.assert_frame_equal(
            specimen_meta_fields_counts, specimen_meta_fields_counts_check
        )

    def test_count_specimen_by_field_value(self):
        specimen_meta_subfields_counts = PMOProcessor.count_specimen_by_field_value(
            self.combined_pmo_data, ["collection_country"]
        )

        specimen_meta_subfields_counts_check = pd.DataFrame(
            [
                {
                    "collection_country": "India",
                    "specimens_count": 1,
                    "specimens_freq": 0.25,
                    "total_specimen_count": 4,
                },
                {
                    "collection_country": "Mozambique",
                    "specimens_count": 2,
                    "specimens_freq": 0.50,
                    "total_specimen_count": 4,
                },
                {
                    "collection_country": "Papua New Guinea",
                    "specimens_count": 1,
                    "specimens_freq": 0.25,
                    "total_specimen_count": 4,
                },
            ]
        )
        pd.testing.assert_frame_equal(
            specimen_meta_subfields_counts, specimen_meta_subfields_counts_check
        )

    def test_extract_allele_counts_freq_from_pmo(self):
        allele_counts = PMOProcessor.extract_allele_counts_freq_from_pmo(
            self.combined_pmo_data
        )
        df_string = pd.util.hash_pandas_object(allele_counts).values
        self.assertEqual(
            hashlib.md5(df_string).hexdigest(), "79fadd12ff7033c23810afa6cd46b4a5"
        )

        allele_counts_collapsed = PMOProcessor.extract_allele_counts_freq_from_pmo(
            self.combined_pmo_data, collapse_across_runs=True
        )
        allele_counts_collapsed_string = pd.util.hash_pandas_object(
            allele_counts_collapsed
        ).values
        self.assertEqual(
            hashlib.md5(allele_counts_collapsed_string).hexdigest(),
            "cb34c7e1357e2e35024a89464b63f06c",
        )

    def test_extract_from_pmo_with_read_filter(self):
        pmo_data_filtered = PMOProcessor.extract_from_pmo_with_read_filter(
            self.combined_pmo_data, 1000
        )
        output_fnp = os.path.join(
            self.test_dir.name, "read_filter_100_combined_pmo_example.json"
        )
        with open(output_fnp, "w") as f:
            json.dump(pmo_data_filtered, f)
        self.assertEqual("879b1a0c62a77a8fcc910a21d8ec9de6", md5sum_of_fnp(output_fnp))
        checker = PMOChecker(self.pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_filtered)

    def test_filter_pmo_by_target_ids(self):
        pmo_data_select_targets = PMOProcessor.filter_pmo_by_target_ids(
            self.combined_pmo_data, {1, 10, 11, 55}
        )
        output_fnp = os.path.join(
            self.test_dir.name, "combined_pmo_target_ids_1_10_11_55.json"
        )
        with open(output_fnp, "w") as f:
            json.dump(pmo_data_select_targets, f)
        self.assertEqual("d9525c0a067d8a291af8bf89b80bfb42", md5sum_of_fnp(output_fnp))
        # check pmo extracted against PMO schema
        checker = PMOChecker(self.pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_select_targets)

    def test_filter_pmo_by_target_names(self):
        pmo_data_select_targets = PMOProcessor.filter_pmo_by_target_names(
            self.combined_pmo_data, {"t96", "t80", "t34", "t55"}
        )
        output_fnp = os.path.join(
            self.test_dir.name, "combined_pmo_target_ids_t96_t80_t34_t55.json"
        )
        with open(output_fnp, "w") as f:
            json.dump(pmo_data_select_targets, f)
        self.assertEqual("853d0806237d4078abc601f0c7dcdc3c", md5sum_of_fnp(output_fnp))
        # check pmo extracted against PMO schema
        checker = PMOChecker(self.pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_select_targets)

    def test_filter_pmo_by_library_sample_ids(self):
        pmo_data_select_targets = PMOProcessor.filter_pmo_by_library_sample_ids(
            self.combined_pmo_data, {1, 3}
        )
        output_fnp = os.path.join(
            self.test_dir.name, "combined_pmo_example_ids_1_3.json"
        )
        with open(output_fnp, "w") as f:
            json.dump(pmo_data_select_targets, f)
        self.assertEqual("b942daeac359bbe063563f8692986d1c", md5sum_of_fnp(output_fnp))
        # check pmo extracted against PMO schema
        checker = PMOChecker(self.pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_select_targets)

    def test_filter_pmo_by_library_sample_names(self):
        pmo_data_select_library_sample_names = (
            PMOProcessor.filter_pmo_by_library_sample_names(
                self.combined_pmo_data, {"8025874217_lib_name", "XUC009_lib_name"}
            )
        )
        output_fnp = os.path.join(
            self.test_dir.name, "combined_pmo_example_ids_8025874217_XUC009.json"
        )
        with open(output_fnp, "w") as f:
            json.dump(pmo_data_select_library_sample_names, f)
        self.assertEqual("99873ba369fa94722969ab9ce6581287", md5sum_of_fnp(output_fnp))
        # check pmo extracted against PMO schema
        checker = PMOChecker(self.pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_select_library_sample_names)

    def test_filter_pmo_by_specimen_ids(self):
        pmo_data_select_targets = PMOProcessor.filter_pmo_by_specimen_ids(
            self.combined_pmo_data, {0, 2}
        )
        output_fnp = os.path.join(
            self.test_dir.name, "combined_pmo_specimen_ids_0_2.json"
        )
        with open(output_fnp, "w") as f:
            json.dump(pmo_data_select_targets, f)
        self.assertEqual("b5d1d39d508a7c7c1d8ccce65c52bc5d", md5sum_of_fnp(output_fnp))
        # check pmo extracted against PMO schema
        checker = PMOChecker(self.pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_select_targets)

    def test_filter_pmo_by_specimen_names(self):
        pmo_data_select_targets = PMOProcessor.filter_pmo_by_specimen_names(
            self.combined_pmo_data, {"8025874217", "5tbx"}
        )
        output_fnp = os.path.join(
            self.test_dir.name, "combined_pmo_specimen_ids_8025874217_5tbx.json"
        )
        with open(output_fnp, "w") as f:
            json.dump(pmo_data_select_targets, f)
        self.assertEqual("b5d1d39d508a7c7c1d8ccce65c52bc5d", md5sum_of_fnp(output_fnp))
        # check pmo extracted against PMO schema
        checker = PMOChecker(self.pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_select_targets)

    def test_extract_from_pmo_samples_with_meta_groupings(self):
        (
            pmo_data_select_meta,
            group_counts,
        ) = PMOProcessor.extract_from_pmo_samples_with_meta_groupings(
            self.combined_pmo_data, "collection_country=Mozambique"
        )
        output_fnp = os.path.join(
            self.test_dir.name, "combined_pmo_collection_country_Mozambique.json"
        )
        with open(output_fnp, "w") as f:
            json.dump(pmo_data_select_meta, f)
        self.assertEqual("338e838297d79a8038d8d5392a6ca57e", md5sum_of_fnp(output_fnp))
        # check pmo extracted against PMO schema
        checker = PMOChecker(self.pmo_jsonschema_data)
        checker.validate_pmo_json(pmo_data_select_meta)

    def test_get_sorted_bioinformatics_run_names(self):
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/combined_pmo_example.json"
            )
        ) as f:
            pmo_data_combined = json.load(f)
        names = PMOProcessor.get_sorted_bioinformatics_run_names(pmo_data_combined)
        self.assertEqual(["Mozambique2018-SeekDeep", "PathWeaver-Heome1"], names)

    def test_get_sorted_specimen_names(self):
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/combined_pmo_example.json"
            )
        ) as f:
            pmo_data_combined = json.load(f)
        names = PMOProcessor.get_sorted_specimen_names(pmo_data_combined)
        self.assertEqual(["5tbx", "8025874217", "8025874266", "XUC009"], names)

    def test_get_sorted_library_sample_names(self):
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/combined_pmo_example.json"
            )
        ) as f:
            pmo_data_combined = json.load(f)
        names = PMOProcessor.get_sorted_library_sample_names(pmo_data_combined)
        self.assertEqual(
            [
                "5tbx_lib_name",
                "8025874217_lib_name",
                "8025874266_lib_name",
                "XUC009_lib_name",
            ],
            names,
        )

    def test_sorted_get_target_names(self):
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/combined_pmo_example.json"
            )
        ) as f:
            pmo_data_combined = json.load(f)
        names = PMOProcessor.get_sorted_target_names(pmo_data_combined)
        self.assertEqual(
            [
                "t1",
                "t10",
                "t100",
                "t11",
                "t12",
                "t13",
                "t14",
                "t15",
                "t16",
                "t17",
                "t18",
                "t19",
                "t2",
                "t20",
                "t21",
                "t22",
                "t23",
                "t24",
                "t25",
                "t26",
                "t27",
                "t28",
                "t29",
                "t3",
                "t30",
                "t31",
                "t32",
                "t33",
                "t34",
                "t35",
                "t36",
                "t37",
                "t38",
                "t39",
                "t4",
                "t40",
                "t41",
                "t42",
                "t43",
                "t44",
                "t45",
                "t46",
                "t47",
                "t48",
                "t49",
                "t5",
                "t50",
                "t51",
                "t52",
                "t53",
                "t54",
                "t55",
                "t56",
                "t57",
                "t58",
                "t59",
                "t6",
                "t60",
                "t61",
                "t62",
                "t63",
                "t64",
                "t65",
                "t66",
                "t67",
                "t68",
                "t69",
                "t7",
                "t70",
                "t71",
                "t72",
                "t73",
                "t74",
                "t75",
                "t76",
                "t77",
                "t78",
                "t79",
                "t8",
                "t80",
                "t81",
                "t82",
                "t83",
                "t84",
                "t85",
                "t86",
                "t87",
                "t88",
                "t89",
                "t9",
                "t90",
                "t91",
                "t92",
                "t93",
                "t94",
                "t95",
                "t96",
                "t97",
                "t98",
                "t99",
            ],
            names,
        )

    def test_get_sorted_panel_names(self):
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/combined_pmo_example.json"
            )
        ) as f:
            pmo_data_combined = json.load(f)
        names = PMOProcessor.get_sorted_panel_names(pmo_data_combined)
        self.assertEqual(["heomev1"], names)

    def test_get_bioinformatics_run_names(self):
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/combined_pmo_example.json"
            )
        ) as f:
            pmo_data_combined = json.load(f)
        names = PMOProcessor.get_bioinformatics_run_names(pmo_data_combined)
        self.assertEqual(["Mozambique2018-SeekDeep", "PathWeaver-Heome1"], names)

    def test_get_specimen_names(self):
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/combined_pmo_example.json"
            )
        ) as f:
            pmo_data_combined = json.load(f)
        names = PMOProcessor.get_specimen_names(pmo_data_combined)
        self.assertEqual(["8025874217", "8025874266", "5tbx", "XUC009"], names)

    def test_get_library_sample_names(self):
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/combined_pmo_example.json"
            )
        ) as f:
            pmo_data_combined = json.load(f)
        names = PMOProcessor.get_library_sample_names(pmo_data_combined)
        self.assertEqual(
            [
                "8025874217_lib_name",
                "8025874266_lib_name",
                "5tbx_lib_name",
                "XUC009_lib_name",
            ],
            names,
        )

    def test_get_target_names(self):
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/combined_pmo_example.json"
            )
        ) as f:
            pmo_data_combined = json.load(f)
        names = PMOProcessor.get_target_names(pmo_data_combined)
        self.assertEqual(
            [
                "t96",
                "t95",
                "t94",
                "t50",
                "t36",
                "t88",
                "t86",
                "t85",
                "t81",
                "t80",
                "t82",
                "t83",
                "t98",
                "t78",
                "t74",
                "t23",
                "t73",
                "t79",
                "t71",
                "t69",
                "t66",
                "t70",
                "t65",
                "t64",
                "t6",
                "t55",
                "t54",
                "t37",
                "t53",
                "t40",
                "t52",
                "t99",
                "t51",
                "t33",
                "t62",
                "t45",
                "t77",
                "t15",
                "t5",
                "t39",
                "t7",
                "t42",
                "t11",
                "t21",
                "t38",
                "t44",
                "t35",
                "t87",
                "t19",
                "t34",
                "t97",
                "t58",
                "t28",
                "t76",
                "t32",
                "t89",
                "t67",
                "t16",
                "t26",
                "t41",
                "t68",
                "t31",
                "t4",
                "t60",
                "t3",
                "t61",
                "t93",
                "t30",
                "t84",
                "t72",
                "t22",
                "t63",
                "t25",
                "t14",
                "t20",
                "t43",
                "t9",
                "t56",
                "t12",
                "t1",
                "t48",
                "t29",
                "t2",
                "t90",
                "t8",
                "t47",
                "t46",
                "t59",
                "t92",
                "t57",
                "t18",
                "t75",
                "t27",
                "t49",
                "t17",
                "t13",
                "t24",
                "t91",
                "t10",
                "t100",
            ],
            names,
        )

    def test_get_panel_names(self):
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/combined_pmo_example.json"
            )
        ) as f:
            pmo_data_combined = json.load(f)
        names = PMOProcessor.get_panel_names(pmo_data_combined)
        self.assertEqual(["heomev1"], names)

    def test_count_targets_per_panel_single_panel_single_reaction(self):
        """Test count_targets_per_panel with a single panel containing one reaction"""
        # Create a simple PMO with one panel and one reaction
        test_pmo = {
            "panel_info": [
                {
                    "panel_name": "test_panel_1",
                    "reactions": [
                        {
                            "reaction_name": "reaction_1",
                            "panel_targets": [0, 1, 2, 3, 4],
                        }
                    ],
                }
            ]
        }

        result = PMOProcessor.count_targets_per_panel(test_pmo)

        expected_data = {"panel_name": ["test_panel_1"], "panel_target_count": [5]}
        expected_df = pd.DataFrame(expected_data)

        pd.testing.assert_frame_equal(result, expected_df)

    def test_count_targets_per_panel_single_panel_multiple_reactions(self):
        """Test count_targets_per_panel with a single panel containing multiple reactions"""
        # Create a PMO with one panel and multiple reactions
        test_pmo = {
            "panel_info": [
                {
                    "panel_name": "test_panel_1",
                    "reactions": [
                        {"reaction_name": "reaction_1", "panel_targets": [0, 1, 2]},
                        {"reaction_name": "reaction_2", "panel_targets": [3, 4, 5, 6]},
                        {"reaction_name": "reaction_3", "panel_targets": [7, 8]},
                    ],
                }
            ]
        }

        result = PMOProcessor.count_targets_per_panel(test_pmo)

        expected_data = {
            "panel_name": ["test_panel_1"],
            "panel_target_count": [9],  # 3 + 4 + 2 = 9 total targets
        }
        expected_df = pd.DataFrame(expected_data)

        pd.testing.assert_frame_equal(result, expected_df)

    def test_count_targets_per_panel_multiple_panels(self):
        """Test count_targets_per_panel with multiple panels"""
        # Create a PMO with multiple panels
        test_pmo = {
            "panel_info": [
                {
                    "panel_name": "panel_A",
                    "reactions": [
                        {"reaction_name": "reaction_1", "panel_targets": [0, 1, 2, 3]}
                    ],
                },
                {
                    "panel_name": "panel_B",
                    "reactions": [
                        {"reaction_name": "reaction_1", "panel_targets": [4, 5]},
                        {
                            "reaction_name": "reaction_2",
                            "panel_targets": [6, 7, 8, 9, 10],
                        },
                    ],
                },
                {
                    "panel_name": "panel_C",
                    "reactions": [
                        {"reaction_name": "reaction_1", "panel_targets": [11]}
                    ],
                },
            ]
        }

        result = PMOProcessor.count_targets_per_panel(test_pmo)

        expected_data = {
            "panel_name": ["panel_A", "panel_B", "panel_C"],
            "panel_target_count": [4, 7, 1],  # 4, (2+5), 1
        }
        expected_df = pd.DataFrame(expected_data)

        pd.testing.assert_frame_equal(result, expected_df)

    def test_count_targets_per_panel_with_real_data_structure(self):
        """Test count_targets_per_panel with realistic PMO data structure"""
        # Create a realistic PMO structure similar to the real data
        # Based on the real data structure with one panel containing 100 targets
        test_pmo = {
            "panel_info": [
                {
                    "panel_name": "heomev1",
                    "reactions": [
                        {
                            "panel_targets": list(range(100)),  # targets 0-99
                            "reaction_name": "full",
                        }
                    ],
                }
            ]
        }

        result = PMOProcessor.count_targets_per_panel(test_pmo)

        expected_data = {
            "panel_name": ["heomev1"],
            "panel_target_count": [100],  # targets 0-99 = 100 total
        }
        expected_df = pd.DataFrame(expected_data)

        pd.testing.assert_frame_equal(result, expected_df)


if __name__ == "__main__":
    unittest.main()
