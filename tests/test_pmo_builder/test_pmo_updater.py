#!/usr/bin/env python3

import os
import unittest
import pandas as pd
from pmotools.pmo_builder.pmo_updater import PMOUpdater


class TestPMOUpdater(unittest.TestCase):
    def setUp(self):
        self.working_dir = os.path.dirname(os.path.abspath(__file__))

        self.specimen_main_list = [
            {"specimen_name": "specimen_001"},
            {"specimen_name": "specimen_002"},
            {"specimen_name": "specimen_003"},
            {"specimen_name": "specimen_004"},
            {"specimen_name": "specimen_005"},
        ]
        self.specimen_update_list = [
            {
                "specimen_name": "specimen_001",
                "collection_date": "2024-03-15",
                "collection_country": "Uganda",
            },
            {
                "specimen_name": "specimen_002",
                "collection_date": "2024-07-22",
                "collection_country": "Kenya",
            },
            {
                "specimen_name": "specimen_003",
                "collection_date": "2025-01-08",
                "collection_country": "Ethiopia",
            },
            {
                "specimen_name": "specimen_004",
                "collection_date": "2024-11-30",
                "collection_country": "Uganda",
            },
            {
                "specimen_name": "specimen_005",
                "collection_date": "2025-02-14",
                "collection_country": "Kenya",
            },
        ]

    def test_check_if_date_yyyy_mm_or_yyyy_mm_dd(self):
        self.assertFalse(PMOUpdater.check_if_date_yyyy_mm_or_yyyy_mm_dd("2023/11/24"))
        self.assertFalse(PMOUpdater.check_if_date_yyyy_mm_or_yyyy_mm_dd("11-24-2023"))
        self.assertFalse(PMOUpdater.check_if_date_yyyy_mm_or_yyyy_mm_dd("invalid-date"))

        self.assertTrue(PMOUpdater.check_if_date_yyyy_mm_or_yyyy_mm_dd("2023-11-24"))
        self.assertTrue(PMOUpdater.check_if_date_yyyy_mm_or_yyyy_mm_dd("2023-11"))

    def test_update_specimen_meta_with_traveler_info(self):
        test_pmo = {
            "specimen_info": [{"specimen_name": "spec1"}, {"specimen_name": "spec2"}],
        }
        traveler_info = pd.DataFrame(
            {
                "specimen_name": ["spec1", "spec1", "spec2"],
                "travel_country": ["Kenya", "Kenya", "Tanzania"],
                "travel_start_date": ["2024-01", "2024-04", "2024-02-15"],
                "travel_end_date": ["2024-02", "2024-06", "2024-02-27"],
            }
        )

        PMOUpdater.update_specimen_meta_with_traveler_info(test_pmo, traveler_info)
        test_out_pmo = {
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
        self.assertEqual(test_out_pmo, test_pmo)

    def test_update_specimen_meta_with_traveler_info_raises(self):
        test_pmo = {
            "specimen_info": [{"specimen_name": "spec1"}, {"specimen_name": "spec2"}],
        }
        traveler_info = pd.DataFrame(
            {
                "specimen_name": ["spec1", "spec2"],
                "travel_country": ["Kenya", "Tanzania"],
                "travel_start_date": ["24-01", "2024-02"],  # BAD: "24-01"
                "travel_end_date": ["2024-02-05", "2024-03"],
            }
        )

        with self.assertRaises(ValueError):
            PMOUpdater.update_specimen_meta_with_traveler_info(test_pmo, traveler_info)

    def test_update_specimen_meta_with_traveler_info_with_optional(self):
        test_pmo = {
            "specimen_info": [{"specimen_name": "spec1"}, {"specimen_name": "spec2"}],
        }
        traveler_info = pd.DataFrame(
            {
                "specimen_name": ["spec1", "spec2"],
                "travel_country": ["Kenya", "Tanzania"],
                "travel_start_date": ["2024-01", "2024-02"],
                "travel_end_date": ["2024-01-20", "2024-02-15"],
                "bed_net": [0.50, 0.0],
                "admin1": ["Nairobi", "Dar es Salaam"],
                "admin2": ["SubCounty1", "SubCounty2"],
                "admin3": ["Ward1", "Ward2"],
                "latlon": ["-1.2921,36.8219", "-6.7924,39.2083"],
            }
        )

        PMOUpdater.update_specimen_meta_with_traveler_info(
            test_pmo,
            traveler_info,
            bed_net_usage_col="bed_net",
            geo_admin1_col="admin1",
            geo_admin2_col="admin2",
            geo_admin3_col="admin3",
            lat_lon_col="latlon",
        )
        test_out_pmo = {
            "specimen_info": [
                {
                    "specimen_name": "spec1",
                    "travel_out_six_month": [
                        {
                            "travel_country": "Kenya",
                            "travel_start_date": "2024-01",
                            "travel_end_date": "2024-01-20",
                            "bed_net": 0.5,
                            "admin1": "Nairobi",
                            "admin2": "SubCounty1",
                            "admin3": "Ward1",
                            "latlon": "-1.2921,36.8219",
                        }
                    ],
                },
                {
                    "specimen_name": "spec2",
                    "travel_out_six_month": [
                        {
                            "travel_country": "Tanzania",
                            "travel_start_date": "2024-02",
                            "travel_end_date": "2024-02-15",
                            "bed_net": 0.0,
                            "admin1": "Dar es Salaam",
                            "admin2": "SubCounty2",
                            "admin3": "Ward2",
                            "latlon": "-6.7924,39.2083",
                        }
                    ],
                },
            ]
        }
        self.assertEqual(test_out_pmo, test_pmo)

    def test_update_specimen_meta_with_traveler_info_with_optional_replace_old(self):
        test_pmo = {
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
                {"specimen_name": "spec2"},
            ],
        }
        traveler_info = pd.DataFrame(
            {
                "specimen_name": ["spec1", "spec2"],
                "travel_country": ["Kenya", "Tanzania"],
                "travel_start_date": ["2024-01", "2024-02"],
                "travel_end_date": ["2024-01-20", "2024-02-15"],
                "bed_net": [0.50, 0.0],
                "admin1": ["Nairobi", "Dar es Salaam"],
                "admin2": ["SubCounty1", "SubCounty2"],
                "admin3": ["Ward1", "Ward2"],
                "latlon": ["-1.2921,36.8219", "-6.7924,39.2083"],
            }
        )

        PMOUpdater.update_specimen_meta_with_traveler_info(
            test_pmo,
            traveler_info,
            bed_net_usage_col="bed_net",
            geo_admin1_col="admin1",
            geo_admin2_col="admin2",
            geo_admin3_col="admin3",
            lat_lon_col="latlon",
            replace_current_traveler_info=True,
        )
        test_out_pmo = {
            "specimen_info": [
                {
                    "specimen_name": "spec1",
                    "travel_out_six_month": [
                        {
                            "travel_country": "Kenya",
                            "travel_start_date": "2024-01",
                            "travel_end_date": "2024-01-20",
                            "bed_net": 0.5,
                            "admin1": "Nairobi",
                            "admin2": "SubCounty1",
                            "admin3": "Ward1",
                            "latlon": "-1.2921,36.8219",
                        }
                    ],
                },
                {
                    "specimen_name": "spec2",
                    "travel_out_six_month": [
                        {
                            "travel_country": "Tanzania",
                            "travel_start_date": "2024-02",
                            "travel_end_date": "2024-02-15",
                            "bed_net": 0.0,
                            "admin1": "Dar es Salaam",
                            "admin2": "SubCounty2",
                            "admin3": "Ward2",
                            "latlon": "-6.7924,39.2083",
                        }
                    ],
                },
            ]
        }
        self.assertEqual(test_out_pmo, test_pmo)

    # PMOUpdater.merge_dicts_by_key
    def test_merge_dicts_by_key_correct_fields_added(self):
        result = PMOUpdater.merge_dicts_by_key(
            self.specimen_main_list,
            self.specimen_update_list,
            key_field="specimen_name",
        )
        result_map = {r["specimen_name"]: r for r in result}
        self.assertEqual(result_map["specimen_001"]["collection_date"], "2024-03-15")
        self.assertEqual(result_map["specimen_001"]["collection_country"], "Uganda")
        self.assertEqual(result_map["specimen_002"]["collection_date"], "2024-07-22")
        self.assertEqual(result_map["specimen_002"]["collection_country"], "Kenya")
        self.assertEqual(result_map["specimen_003"]["collection_date"], "2025-01-08")
        self.assertEqual(result_map["specimen_003"]["collection_country"], "Ethiopia")

    def test_merge_dicts_by_key_all_specimens_present_in_result(self):
        result = PMOUpdater.merge_dicts_by_key(
            self.specimen_main_list,
            self.specimen_update_list,
            key_field="specimen_name",
        )
        result_names = {r["specimen_name"] for r in result}
        expected_names = {r["specimen_name"] for r in self.specimen_main_list}
        self.assertEqual(result_names, expected_names)

    def test_merge_dicts_by_key_does_not_mutate_main_list(self):
        import copy

        main_copy = copy.deepcopy(self.specimen_main_list)
        PMOUpdater.merge_dicts_by_key(
            self.specimen_main_list,
            self.specimen_update_list,
            key_field="specimen_name",
        )
        self.assertEqual(self.specimen_main_list, main_copy)

    def test_merge_dicts_by_key_does_not_mutate_update_list(self):
        import copy

        update_copy = copy.deepcopy(self.specimen_update_list)
        PMOUpdater.merge_dicts_by_key(
            self.specimen_main_list,
            self.specimen_update_list,
            key_field="specimen_name",
        )
        self.assertEqual(self.specimen_update_list, update_copy)

    def test_merge_dicts_by_key_partial_update_only_updates_provided_specimens(self):
        # update_list covering only a subset of main_list should leave others untouched."""
        partial_update = [
            {
                "specimen_name": "specimen_001",
                "collection_date": "2024-03-15",
                "collection_country": "Uganda",
            },
            {
                "specimen_name": "specimen_003",
                "collection_date": "2025-01-08",
                "collection_country": "Ethiopia",
            },
        ]
        result = PMOUpdater.merge_dicts_by_key(
            self.specimen_main_list, partial_update, key_field="specimen_name"
        )
        result_map = {r["specimen_name"]: r for r in result}
        self.assertIn("collection_date", result_map["specimen_001"])
        self.assertIn("collection_date", result_map["specimen_003"])
        self.assertNotIn("collection_date", result_map["specimen_002"])
        self.assertNotIn("collection_date", result_map["specimen_004"])
        self.assertNotIn("collection_date", result_map["specimen_005"])

    # PMOUpdater.merge_dicts_by_key, testing replacement
    def test_merge_dicts_by_key_replace_true_overwrites_existing_field(self):
        main_with_existing = [
            {"specimen_name": "specimen_001", "collection_country": "Uganda"},
            {"specimen_name": "specimen_002", "collection_country": "Kenya"},
        ]
        update = [
            {"specimen_name": "specimen_001", "collection_country": "Ethiopia"},
            {"specimen_name": "specimen_002", "collection_country": "Uganda"},
        ]
        result = PMOUpdater.merge_dicts_by_key(
            main_with_existing, update, key_field="specimen_name", replace=True
        )
        result_map = {r["specimen_name"]: r for r in result}
        self.assertEqual(result_map["specimen_001"]["collection_country"], "Ethiopia")
        self.assertEqual(result_map["specimen_002"]["collection_country"], "Uganda")

    # PMOUpdater.merge_dicts_by_key test ignoring fields

    def test_merge_dicts_by_key_ignore_fields_not_added(self):
        result = PMOUpdater.merge_dicts_by_key(
            self.specimen_main_list,
            self.specimen_update_list,
            key_field="specimen_name",
            ignore_fields=["collection_country"],
        )
        result_map = {r["specimen_name"]: r for r in result}
        self.assertIn("collection_date", result_map["specimen_001"])
        self.assertNotIn("collection_country", result_map["specimen_001"])
        self.assertNotIn("collection_country", result_map["specimen_003"])

    def test_merge_dicts_by_key_ignore_fields_does_not_affect_other_fields(self):
        result = PMOUpdater.merge_dicts_by_key(
            self.specimen_main_list,
            self.specimen_update_list,
            key_field="specimen_name",
            ignore_fields=["collection_country"],
        )
        result_map = {r["specimen_name"]: r for r in result}
        self.assertEqual(result_map["specimen_002"]["collection_date"], "2024-07-22")

    # PMOUpdater.merge_dicts_by_key testing for expected errors
    def test_merge_dicts_by_key_replace_false_raises_on_existing_field(self):
        # replace=False should raise ValueError when update would overwrite an existing value
        main_with_existing = [
            {"specimen_name": "specimen_001", "collection_country": "Uganda"},
            {"specimen_name": "specimen_002"},
        ]
        update = [
            {"specimen_name": "specimen_001", "collection_country": "Kenya"},
        ]
        with self.assertRaises(ValueError) as context:
            PMOUpdater.merge_dicts_by_key(
                main_with_existing, update, key_field="specimen_name", replace=False
            )
        self.assertIn("collection_country", str(context.exception))

    def test_merge_dicts_by_key_missing_key_field_in_main_raises(self):
        main_missing_key = [
            {"specimen_name": "specimen_001"},
            {"NOT_specimen_name": "specimen_002"},  # missing key
        ]
        with self.assertRaises(KeyError) as context:
            PMOUpdater.merge_dicts_by_key(
                main_missing_key,
                self.specimen_update_list[:1],
                key_field="specimen_name",
            )
        self.assertIn("main_list", str(context.exception))

    def test_merge_dicts_by_key_missing_key_field_in_update_raises(self):
        update_missing_key = [
            {
                "specimen_name": "specimen_001",
                "collection_date": "2024-03-15",
                "collection_country": "Uganda",
            },
            {
                "NOT_specimen_name": "specimen_002",
                "collection_date": "2024-07-22",
            },  # missing key
        ]
        with self.assertRaises(KeyError) as context:
            PMOUpdater.merge_dicts_by_key(
                self.specimen_main_list, update_missing_key, key_field="specimen_name"
            )
        self.assertIn("update_list", str(context.exception))

    def test_merge_dicts_by_key_duplicate_in_main_raises(self):
        main_with_dupes = [
            {"specimen_name": "specimen_001"},
            {"specimen_name": "specimen_001"},  # duplicate
            {"specimen_name": "specimen_002"},
        ]
        with self.assertRaises(ValueError) as context:
            PMOUpdater.merge_dicts_by_key(
                main_with_dupes,
                self.specimen_update_list[:1],
                key_field="specimen_name",
            )
        self.assertIn("specimen_001", str(context.exception))

    def test_merge_dicts_by_key_duplicate_in_update_raises(self):
        update_with_dupes = [
            {
                "specimen_name": "specimen_001",
                "collection_date": "2024-03-15",
                "collection_country": "Uganda",
            },
            {
                "specimen_name": "specimen_001",
                "collection_date": "2024-05-10",
                "collection_country": "Kenya",
            },  # duplicate
        ]
        with self.assertRaises(ValueError) as context:
            PMOUpdater.merge_dicts_by_key(
                self.specimen_main_list, update_with_dupes, key_field="specimen_name"
            )
        self.assertIn("specimen_001", str(context.exception))

    def test_merge_dicts_by_key_update_key_not_in_main_raises(self):
        update_with_unknown = [
            {
                "specimen_name": "specimen_001",
                "collection_date": "2024-03-15",
                "collection_country": "Uganda",
            },
            {
                "specimen_name": "specimen_UNKNOWN",
                "collection_date": "2024-06-01",
                "collection_country": "Kenya",
            },
        ]
        with self.assertRaises(KeyError) as context:
            PMOUpdater.merge_dicts_by_key(
                self.specimen_main_list, update_with_unknown, key_field="specimen_name"
            )
        self.assertIn("specimen_UNKNOWN", str(context.exception))

    # PMOUpdater.merge_dicts_by_key edge cases

    def test_merge_dicts_by_key_empty_update_list_returns_main_unchanged(self):
        result = PMOUpdater.merge_dicts_by_key(
            self.specimen_main_list, [], key_field="specimen_name"
        )
        result_names = sorted(r["specimen_name"] for r in result)
        main_names = sorted(r["specimen_name"] for r in self.specimen_main_list)
        self.assertListEqual(result_names, main_names)
        for r in result:
            self.assertNotIn("collection_date", r)
            self.assertNotIn("collection_country", r)

    def test_merge_dicts_by_key_empty_main_and_update_returns_empty(self):
        result = PMOUpdater.merge_dicts_by_key([], [], key_field="specimen_name")
        self.assertListEqual(result, [])


if __name__ == "__main__":
    unittest.main()
