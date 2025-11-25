#!/usr/bin/env python3

import os
import unittest
import pandas as pd
from pmotools.pmo_builder.pmo_updater import PMOUpdater


class TestPMOUpdater(unittest.TestCase):
    def setUp(self):
        self.working_dir = os.path.dirname(os.path.abspath(__file__))

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


if __name__ == "__main__":
    unittest.main()
