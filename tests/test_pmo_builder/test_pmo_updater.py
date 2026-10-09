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
        self.assertFalse(PMOUpdater._check_if_date_yyyy_mm_or_yyyy_mm_dd("2023/11/24"))
        self.assertFalse(PMOUpdater._check_if_date_yyyy_mm_or_yyyy_mm_dd("11-24-2023"))
        self.assertFalse(
            PMOUpdater._check_if_date_yyyy_mm_or_yyyy_mm_dd("invalid-date")
        )

        self.assertTrue(PMOUpdater._check_if_date_yyyy_mm_or_yyyy_mm_dd("2023-11-24"))
        self.assertTrue(PMOUpdater._check_if_date_yyyy_mm_or_yyyy_mm_dd("2023-11"))

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
                            "bed_net_usage": 0.5,
                            "geo_admin1": "Nairobi",
                            "geo_admin2": "SubCounty1",
                            "geo_admin3": "Ward1",
                            "lat_lon": "-1.2921,36.8219",
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
                            "bed_net_usage": 0.0,
                            "geo_admin1": "Dar es Salaam",
                            "geo_admin2": "SubCounty2",
                            "geo_admin3": "Ward2",
                            "lat_lon": "-6.7924,39.2083",
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
                            "bed_net_usage": 0.5,
                            "geo_admin1": "Nairobi",
                            "geo_admin2": "SubCounty1",
                            "geo_admin3": "Ward1",
                            "lat_lon": "-1.2921,36.8219",
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
                            "bed_net_usage": 0.0,
                            "geo_admin1": "Dar es Salaam",
                            "geo_admin2": "SubCounty2",
                            "geo_admin3": "Ward2",
                            "lat_lon": "-6.7924,39.2083",
                        }
                    ],
                },
            ]
        }
        self.assertEqual(test_out_pmo, test_pmo)

    def test_update_specimen_meta_with_traveler_info_skips_blank_optional(self):
        test_pmo = {"specimen_info": [{"specimen_name": "spec1"}]}
        traveler_info = pd.DataFrame(
            {
                "specimen_name": ["spec1"],
                "travel_country": ["Kenya"],
                "travel_start_date": ["2024-01"],
                "travel_end_date": ["2024-02"],
                "bed_net": [None],
                "admin1": [""],
            }
        )
        PMOUpdater.update_specimen_meta_with_traveler_info(
            test_pmo,
            traveler_info,
            bed_net_usage_col="bed_net",
            geo_admin1_col="admin1",
        )
        self.assertEqual(
            test_pmo["specimen_info"][0]["travel_out_six_month"],
            [
                {
                    "travel_country": "Kenya",
                    "travel_start_date": "2024-01",
                    "travel_end_date": "2024-02",
                }
            ],
        )

    def test_update_target_info_with_markers_of_interest(self):
        test_pmo = {
            "target_info": [{"target_name": "t1"}, {"target_name": "t2"}],
        }
        markers_info = pd.DataFrame(
            {
                "target_name": ["t1", "t1"],
                "chrom": ["Pf3D7_04_v3", "Pf3D7_04_v3"],
                "start": [748238, 748409],
                "end": [748239, 748410],
                "ref_seq": ["A", ""],
                "associations": ["SP resistance, dhfr", None],
            }
        )
        PMOUpdater.update_target_info_with_markers_of_interest(
            test_pmo,
            markers_info,
            ref_seq_col="ref_seq",
            associations_col="associations",
        )
        self.assertEqual(
            test_pmo["target_info"][0]["markers_of_interest"],
            [
                {
                    "marker_location": {
                        "genome_id": 0,
                        "chrom": "Pf3D7_04_v3",
                        "start": 748238,
                        "end": 748239,
                        "ref_seq": "A",
                    },
                    "associations": ["SP resistance", "dhfr"],
                },
                {
                    "marker_location": {
                        "genome_id": 0,
                        "chrom": "Pf3D7_04_v3",
                        "start": 748409,
                        "end": 748410,
                    }
                },
            ],
        )
        self.assertNotIn("markers_of_interest", test_pmo["target_info"][1])

    def test_update_target_info_with_markers_of_interest_replace(self):
        old_marker = {
            "marker_location": {"genome_id": 0, "chrom": "c", "start": 1, "end": 2}
        }
        markers_info = pd.DataFrame(
            {"target_name": ["t1"], "chrom": ["c"], "start": [5], "end": [6]}
        )
        for replace, expected_count in ((False, 2), (True, 1)):
            test_pmo = {
                "target_info": [
                    {"target_name": "t1", "markers_of_interest": [dict(old_marker)]}
                ]
            }
            PMOUpdater.update_target_info_with_markers_of_interest(
                test_pmo, markers_info, replace_current_markers=replace
            )
            self.assertEqual(
                len(test_pmo["target_info"][0]["markers_of_interest"]),
                expected_count,
            )

    def test_update_target_info_with_markers_of_interest_raises(self):
        test_pmo = {"target_info": [{"target_name": "t1"}]}
        unknown_target = pd.DataFrame(
            {"target_name": ["t9"], "chrom": ["c"], "start": [1], "end": [2]}
        )
        with self.assertRaises(ValueError):
            PMOUpdater.update_target_info_with_markers_of_interest(
                test_pmo, unknown_target
            )
        blank_start = pd.DataFrame(
            {"target_name": ["t1"], "chrom": ["c"], "start": [None], "end": [2]}
        )
        with self.assertRaises(ValueError):
            PMOUpdater.update_target_info_with_markers_of_interest(
                test_pmo, blank_start
            )
        with self.assertRaises(ValueError):
            PMOUpdater.update_target_info_with_markers_of_interest(
                test_pmo, unknown_target, chrom_col="not_a_column"
            )

    def _rep_mhap_pmo(self):
        return {
            "target_info": [{"target_name": "t1"}, {"target_name": "t2"}],
            "representative_microhaplotypes": {
                "targets": [
                    {
                        "target_id": 1,
                        "microhaplotypes": [{"seq": "ACGT"}, {"seq": "ACTT"}],
                    }
                ]
            },
        }

    def test_update_representative_microhaplotypes_with_seq_variants(self):
        test_pmo = self._rep_mhap_pmo()
        seq_variants_info = pd.DataFrame(
            {
                "target_name": ["t2", "t2"],
                "seq": ["ACTT", "ACTT"],
                "chrom": ["c", "c"],
                "start": [2, 3],
                "end": [3, 4],
                "ref_seq": ["G", "T"],
                "alt_seq": ["T", None],
            }
        )
        PMOUpdater.update_representative_microhaplotypes_with_seq_variants(
            test_pmo,
            seq_variants_info,
            ref_seq_col="ref_seq",
            alt_seq_col="alt_seq",
        )
        mhaps = test_pmo["representative_microhaplotypes"]["targets"][0][
            "microhaplotypes"
        ]
        self.assertNotIn("associated_seq_variants", mhaps[0])
        self.assertEqual(
            mhaps[1]["associated_seq_variants"],
            [
                {
                    "genome_id": 0,
                    "chrom": "c",
                    "start": 2,
                    "end": 3,
                    "ref_seq": "G",
                    "alt_seq": "T",
                },
                {"genome_id": 0, "chrom": "c", "start": 3, "end": 4, "ref_seq": "T"},
            ],
        )

    def test_update_representative_microhaplotypes_with_seq_variants_raises(self):
        for target_name, seq in (("t2", "GGGG"), ("t1", "ACGT")):
            seq_variants_info = pd.DataFrame(
                {
                    "target_name": [target_name],
                    "seq": [seq],
                    "chrom": ["c"],
                    "start": [1],
                    "end": [2],
                }
            )
            with self.assertRaises(ValueError):
                PMOUpdater.update_representative_microhaplotypes_with_seq_variants(
                    self._rep_mhap_pmo(), seq_variants_info
                )

    def test_update_representative_microhaplotypes_with_protein_variants(self):
        test_pmo = self._rep_mhap_pmo()
        protein_variants_info = pd.DataFrame(
            {
                "target_name": ["t2", "t2"],
                "seq": ["ACGT", "ACGT"],
                "transcript": ["PF3D7_0417200.1", "PF3D7_0417200.1"],
                "protein_start": [50, 107],
                "protein_end": [51, 108],
                "protein_ref": ["C", "S"],
                "protein_alt": ["R", "N"],
                "gene": ["dhfr", ""],
                "codon_chrom": ["Pf3D7_04_v3", None],
                "codon_start": [748238, None],
                "codon_end": [748241, None],
            }
        )
        PMOUpdater.update_representative_microhaplotypes_with_protein_variants(
            test_pmo,
            protein_variants_info,
            protein_ref_seq_col="protein_ref",
            protein_alt_seq_col="protein_alt",
            gene_name_col="gene",
            codon_chrom_col="codon_chrom",
            codon_start_col="codon_start",
            codon_end_col="codon_end",
        )
        mhaps = test_pmo["representative_microhaplotypes"]["targets"][0][
            "microhaplotypes"
        ]
        self.assertEqual(
            mhaps[0]["associated_protein_variants"],
            [
                {
                    "protein_location": {
                        "genome_id": 0,
                        "chrom": "PF3D7_0417200.1",
                        "start": 50,
                        "end": 51,
                        "ref_seq": "C",
                        "alt_seq": "R",
                    },
                    "gene_name": "dhfr",
                    "codon_genomic_location": {
                        "genome_id": 0,
                        "chrom": "Pf3D7_04_v3",
                        "start": 748238,
                        "end": 748241,
                    },
                },
                {
                    "protein_location": {
                        "genome_id": 0,
                        "chrom": "PF3D7_0417200.1",
                        "start": 107,
                        "end": 108,
                        "ref_seq": "S",
                        "alt_seq": "N",
                    },
                },
            ],
        )
        self.assertNotIn("associated_protein_variants", mhaps[1])

    def test_update_representative_microhaplotypes_with_protein_variants_partial_codon_cols(
        self,
    ):
        protein_variants_info = pd.DataFrame(
            {
                "target_name": ["t2"],
                "seq": ["ACGT"],
                "transcript": ["tx"],
                "protein_start": [1],
                "protein_end": [2],
                "codon_chrom": ["c"],
            }
        )
        with self.assertRaises(ValueError):
            PMOUpdater.update_representative_microhaplotypes_with_protein_variants(
                self._rep_mhap_pmo(),
                protein_variants_info,
                codon_chrom_col="codon_chrom",
            )

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
