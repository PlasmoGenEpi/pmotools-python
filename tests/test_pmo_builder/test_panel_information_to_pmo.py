import unittest
import pandas as pd
from unittest.mock import patch

from pmotools.pmo_builder.panel_information_to_pmo import (
    PMOPanelBuilder,
    check_genome_info,
    merge_panel_info_dicts,
    panel_info_table_to_pmo,
)


class TestPanelInformationToPMO(unittest.TestCase):
    def setUp(self):
        self.min_target_table = pd.DataFrame(
            {
                "target_name": ["target1", "target2", "target3"],
                "fwd_primer": ["CTT", "CTA", "TTG"],
                "rev_primer": ["GTT", "TGG", "ATT"],
            }
        )
        self.min_target_json = [
            {
                "target_name": "target1",
                "forward_primer": {"seq": "CTT"},
                "reverse_primer": {"seq": "GTT"},
            },
            {
                "target_name": "target2",
                "forward_primer": {"seq": "CTA"},
                "reverse_primer": {"seq": "TGG"},
            },
            {
                "target_name": "target3",
                "forward_primer": {"seq": "TTG"},
                "reverse_primer": {"seq": "ATT"},
            },
        ]
        self.genome_info = {
            "name": "3D7",
            "url": "somthing.com",
            "genome_version": 1,
            "taxon_id": 1,
        }
        self.min_builder = PMOPanelBuilder(
            self.min_target_table,
            "test_panel",
            self.genome_info,
            "target_name",
            "fwd_primer",
            "rev_primer",
        )

        min_target_table_with_locations = self.min_target_table.copy()
        min_target_table_with_locations["target_start"] = 1
        min_target_table_with_locations["target_end"] = 5
        min_target_table_with_locations["insert_start"] = 2
        min_target_table_with_locations["insert_end"] = 3
        min_target_table_with_locations["chrom"] = "chrom1"
        self.builder_with_locations = PMOPanelBuilder(
            min_target_table_with_locations,
            "test_panel",
            self.genome_info,
            "target_name",
            "fwd_primer",
            "rev_primer",
            forward_primers_start_col="target_start",
            forward_primers_end_col="insert_start",
            reverse_primers_start_col="insert_end",
            reverse_primers_end_col="target_end",
            insert_start_col="insert_start",
            insert_end_col="insert_end",
            chrom_col="chrom",
        )

        min_target_table_with_locations_missing = self.min_target_table.copy()
        min_target_table_with_locations_missing["target_start"] = [1, None, 1]
        min_target_table_with_locations_missing["target_end"] = None
        min_target_table_with_locations_missing["insert_start"] = [2, None, None]
        min_target_table_with_locations_missing["insert_end"] = [3, 3, 3]
        min_target_table_with_locations_missing["chrom"] = "chrom1"

        self.builder_with_missing_locations = PMOPanelBuilder(
            min_target_table_with_locations_missing,
            "test_panel",
            self.genome_info,
            "target_name",
            "fwd_primer",
            "rev_primer",
            forward_primers_start_col="target_start",
            forward_primers_end_col="insert_start",
            reverse_primers_start_col="insert_end",
            reverse_primers_end_col="target_end",
            insert_start_col="insert_start",
            insert_end_col="insert_end",
            chrom_col="chrom",
        )

    def test_check_targets_are_unique_passes(self):
        self.min_builder.check_targets_are_unique()

    def test_check_targets_are_unique_fails(self):
        new_row = pd.DataFrame(
            {"target_name": ["target1"], "fwd_primer": ["TTT"], "rev_primer": ["GGG"]}
        )
        duplicate_target_table = pd.concat([self.min_target_table, new_row])
        builder = PMOPanelBuilder(
            duplicate_target_table,
            "test_panel",
            self.genome_info,
            "target_name",
            "fwd_primer",
            "rev_primer",
        )

        with self.assertRaises(ValueError) as context:
            builder.check_targets_are_unique()
        self.assertEqual(
            str(context.exception),
            "The following target_ids are duplicated: ['target1']",
        )

    def test_check_unique_target_info_passes(self):
        self.min_builder.check_unique_target_info(["rev_primer", "fwd_primer"])

    def test_check_unique_target_info_fails(self):
        new_row = pd.DataFrame(
            {"target_name": ["target4"], "fwd_primer": ["CTT"], "rev_primer": ["GTT"]}
        )
        duplicate_target_primers_table = pd.concat([self.min_target_table, new_row])
        builder = PMOPanelBuilder(
            duplicate_target_primers_table,
            "test_panel",
            self.genome_info,
            "target_name",
            "fwd_primer",
            "rev_primer",
        )

        with self.assertRaises(ValueError) as context:
            builder.check_unique_target_info(["rev_primer", "fwd_primer"])
        self.assertEqual(
            str(context.exception),
            "The following targets have duplicated information:\ntargets: target1, target4 → rev_primer=GTT, fwd_primer=CTT",
        )

    def test_summarise_targets_missing_optional_info_nothing_missing(self):
        (
            missing_insert_loc,
            missing_fwd_primer_loc,
            missing_rev_primer_loc,
        ) = self.builder_with_locations.summarise_targets_missing_optional_info()
        self.assertEqual(missing_insert_loc, [])
        self.assertEqual(missing_fwd_primer_loc, [])
        self.assertEqual(missing_rev_primer_loc, [])

    def test_summarise_targets_missing_optional_info_with_missing(self):
        with self.assertWarns(UserWarning):
            (
                missing_insert_loc,
                missing_fwd_primer_loc,
                missing_rev_primer_loc,
            ) = self.builder_with_missing_locations.summarise_targets_missing_optional_info()
        self.assertEqual(missing_insert_loc, ["target2", "target3"])
        self.assertEqual(missing_fwd_primer_loc, ["target2", "target3"])
        self.assertEqual(missing_rev_primer_loc, ["target1", "target2", "target3"])

    def test_check_genome_info_passes(self):
        check_genome_info(self.genome_info)

    def test_check_genome_info_passes_with_list(self):
        genome_info_list = [
            {
                "name": "3D7",
                "url": "somthing.com",
                "genome_version": 1,
                "taxon_id": 1,
            },
            {
                "name": "DD2",
                "url": "somthingelse.com",
                "genome_version": 2,
                "taxon_id": 2,
            },
        ]
        check_genome_info(genome_info_list)

    def test_check_genome_info_fails_with_string(self):
        with self.assertRaises(TypeError) as context:
            check_genome_info("genome_info")
        self.assertIn(
            "genome_info must be a dict or list",
            str(context.exception),
        )

    def test_check_genome_info_fails_with_missing_key(self):
        with self.assertRaises(ValueError):
            check_genome_info({"name": "3D7", "taxon_id": 1})

    def test_check_genome_info_fails_with_empty_list(self):
        with self.assertRaises(ValueError) as context:
            check_genome_info([])
        self.assertEqual(str(context.exception), "genome_info list cannot be empty")

    def test_check_genome_info_fails_with_list_containing_non_dict(self):
        with self.assertRaises(TypeError) as context:
            check_genome_info([self.genome_info, "not_a_dict"])
        self.assertIn("genome_info[1] must be a dict", str(context.exception))

    def test_check_genome_info_fails_with_list_containing_invalid_dict(self):
        with self.assertRaises(ValueError) as context:
            check_genome_info(
                [
                    self.genome_info,
                    {"name": "3D7", "taxon_id": 1},  # missing required keys
                ]
            )
        self.assertIn("genome_info[1] missing required keys", str(context.exception))

    def test_build_panel_info(self):
        panel_info = self.min_builder.build_panel_info(self.min_target_json)
        expected_panel_info = {
            "panel_name": "test_panel",
            "reactions": [{"reaction_name": "1", "panel_targets": [0, 1, 2]}],
        }
        self.assertEqual(panel_info, expected_panel_info)

    def test_merge_panel_info_dicts_no_overlap(self):
        target_info_b = [
            {
                "target_name": "target4",
                "forward_primer": {"seq": "CTA"},
                "reverse_primer": {"seq": "TGG"},
            },
            {
                "target_name": "target5",
                "forward_primer": {"seq": "TTG"},
                "reverse_primer": {"seq": "ATT"},
            },
        ]

        panel_info_dict_a = {
            "panel_info": [
                {
                    "panel_name": "test_panel1",
                    "reactions": [{"reaction_name": "1", "panel_targets": [0, 1, 2]}],
                }
            ],
            "targeted_genomes": [self.genome_info],
            "target_info": self.min_target_json,
        }

        panel_info_dict_b = {
            "panel_info": [
                {
                    "panel_name": "test_panel2",
                    "reactions": [{"reaction_name": "1", "panel_targets": [0, 1]}],
                }
            ],
            "targeted_genomes": [self.genome_info],
            "target_info": target_info_b,
        }

        merged = merge_panel_info_dicts([panel_info_dict_a, panel_info_dict_b])

        expected_merged = {
            "panel_info": [
                {
                    "panel_name": "test_panel1",
                    "reactions": [{"reaction_name": "1", "panel_targets": [0, 1, 2]}],
                },
                {
                    "panel_name": "test_panel2",
                    "reactions": [{"reaction_name": "1", "panel_targets": [3, 4]}],
                },
            ],
            "targeted_genomes": [self.genome_info],
            "target_info": self.min_target_json + target_info_b,
        }

        self.assertEqual(merged, expected_merged)

    def test_merge_panel_info_dicts_with_overlap(self):
        target_info_b = [
            {
                "target_name": "target2",
                "forward_primer": {"seq": "CTA"},
                "reverse_primer": {"seq": "TGG"},
            },
            {
                "target_name": "target5",
                "forward_primer": {"seq": "TTG"},
                "reverse_primer": {"seq": "ATT"},
            },
        ]

        panel_info_dict_a = {
            "panel_info": [
                {
                    "panel_name": "test_panel1",
                    "reactions": [{"reaction_name": "1", "panel_targets": [0, 1, 2]}],
                }
            ],
            "targeted_genomes": [self.genome_info],
            "target_info": self.min_target_json,
        }
        panel_info_dict_b = {
            "panel_info": [
                {
                    "panel_name": "test_panel2",
                    "reactions": [{"reaction_name": "1", "panel_targets": [0, 1]}],
                }
            ],
            "targeted_genomes": [self.genome_info],
            "target_info": target_info_b,
        }

        merged = merge_panel_info_dicts([panel_info_dict_a, panel_info_dict_b])

        expected_merged = {
            "panel_info": [
                {
                    "panel_name": "test_panel1",
                    "reactions": [{"reaction_name": "1", "panel_targets": [0, 1, 2]}],
                },
                {
                    "panel_name": "test_panel2",
                    "reactions": [{"reaction_name": "1", "panel_targets": [1, 3]}],
                },
            ],
            "targeted_genomes": [self.genome_info],
            "target_info": self.min_target_json
            + [
                {
                    "target_name": "target5",
                    "forward_primer": {"seq": "TTG"},
                    "reverse_primer": {"seq": "ATT"},
                },
            ],
        }

        self.assertEqual(merged, expected_merged)

    def test_build_panel_info_multi_reaction(self):
        target_table_with_reactions = self.min_target_table
        target_table_with_reactions["reaction"] = [
            "reaction1",
            "reaction2",
            "reaction2",
        ]
        builder = PMOPanelBuilder(
            target_table_with_reactions,
            "test_panel",
            self.genome_info,
            "target_name",
            "fwd_primer",
            "rev_primer",
            reaction_name_col="reaction",
        )
        panel_info = builder.build_panel_info(self.min_target_json)
        expected_panel_info = {
            "panel_name": "test_panel",
            "reactions": [
                {"reaction_name": "reaction1", "panel_targets": [0]},
                {"reaction_name": "reaction2", "panel_targets": [1, 2]},
            ],
        }
        self.assertEqual(panel_info, expected_panel_info)

    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_targets_are_unique"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_unique_target_info"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.summarise_targets_missing_optional_info"
    )
    def test_create_targets_dict_min_info(
        self,
        mock_summarise_targets_missing_optional_info,
        mock_check_unique_target_info,
        mock_check_targets_are_unique,
    ):
        mock_summarise_targets_missing_optional_info.return_value = [], [], []
        target_info = self.min_builder.create_targets_dict()
        self.assertEqual(self.min_target_json, target_info)

    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_targets_are_unique"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_unique_target_info"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.summarise_targets_missing_optional_info"
    )
    def test_create_targets_dict_full_info(
        self,
        mock_summarise_targets_missing_optional_info,
        mock_check_unique_target_info,
        mock_check_targets_are_unique,
    ):
        mock_summarise_targets_missing_optional_info.return_value = [], [], []
        target_info = self.builder_with_locations.create_targets_dict()
        expected_json = [
            {
                "target_name": "target1",
                "insert_location": {
                    "genome_id": 0,
                    "chrom": "chrom1",
                    "start": 2,
                    "end": 3,
                },
                "forward_primer": {
                    "seq": "CTT",
                    "location": {
                        "genome_id": 0,
                        "chrom": "chrom1",
                        "end": 1,
                        "start": 2,
                    },
                },
                "reverse_primer": {
                    "seq": "GTT",
                    "location": {
                        "genome_id": 0,
                        "chrom": "chrom1",
                        "start": 3,
                        "end": 5,
                    },
                },
            },
            {
                "target_name": "target2",
                "insert_location": {
                    "genome_id": 0,
                    "chrom": "chrom1",
                    "start": 2,
                    "end": 3,
                },
                "forward_primer": {
                    "seq": "CTA",
                    "location": {
                        "genome_id": 0,
                        "chrom": "chrom1",
                        "end": 1,
                        "start": 2,
                    },
                },
                "reverse_primer": {
                    "seq": "TGG",
                    "location": {
                        "genome_id": 0,
                        "chrom": "chrom1",
                        "start": 3,
                        "end": 5,
                    },
                },
            },
            {
                "target_name": "target3",
                "insert_location": {
                    "genome_id": 0,
                    "chrom": "chrom1",
                    "start": 2,
                    "end": 3,
                },
                "forward_primer": {
                    "seq": "TTG",
                    "location": {
                        "genome_id": 0,
                        "chrom": "chrom1",
                        "end": 1,
                        "start": 2,
                    },
                },
                "reverse_primer": {
                    "seq": "ATT",
                    "location": {
                        "genome_id": 0,
                        "chrom": "chrom1",
                        "start": 3,
                        "end": 5,
                    },
                },
            },
        ]
        self.assertEqual(expected_json, target_info)

    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_targets_are_unique"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_unique_target_info"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.summarise_targets_missing_optional_info"
    )
    def test_create_targets_dict_missing_info(
        self,
        mock_summarise_targets_missing_optional_info,
        mock_check_unique_target_info,
        mock_check_targets_are_unique,
    ):
        missing_insert_loc = ["target2", "target3"]
        missing_fwd_primer_loc = ["target2", "target3"]
        missing_rev_primer_loc = ["target1", "target2", "target3"]

        mock_summarise_targets_missing_optional_info.return_value = (
            missing_insert_loc,
            missing_fwd_primer_loc,
            missing_rev_primer_loc,
        )
        target_info = self.builder_with_locations.create_targets_dict()
        expected_json = [
            {
                "target_name": "target1",
                "insert_location": {
                    "genome_id": 0,
                    "chrom": "chrom1",
                    "start": 2,
                    "end": 3,
                },
                "forward_primer": {
                    "seq": "CTT",
                    "location": {
                        "genome_id": 0,
                        "chrom": "chrom1",
                        "end": 1,
                        "start": 2,
                    },
                },
                "reverse_primer": {"seq": "GTT"},
            },
            {
                "target_name": "target2",
                "forward_primer": {"seq": "CTA"},
                "reverse_primer": {"seq": "TGG"},
            },
            {
                "target_name": "target3",
                "forward_primer": {"seq": "TTG"},
                "reverse_primer": {"seq": "ATT"},
            },
        ]
        self.assertEqual(expected_json, target_info)

    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_targets_are_unique"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_unique_target_info"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.summarise_targets_missing_optional_info"
    )
    def test_create_targets_dict_additional_info(
        self,
        mock_summarise_targets_missing_optional_info,
        mock_check_unique_target_info,
        mock_check_targets_are_unique,
    ):
        target_table = self.min_target_table.copy()
        target_table["extra_col"] = "my target description"
        target_table["special_col"] = "special target"
        builder = PMOPanelBuilder(
            target_table,
            "test_panel",
            self.genome_info,
            additional_target_info_cols=["extra_col", "special_col"],
        )
        mock_summarise_targets_missing_optional_info.return_value = [], [], []
        target_info = builder.create_targets_dict()
        expected_json = [
            {
                "target_name": "target1",
                "forward_primer": {"seq": "CTT"},
                "reverse_primer": {"seq": "GTT"},
                "extra_col": "my target description",
                "special_col": "special target",
            },
            {
                "target_name": "target2",
                "forward_primer": {"seq": "CTA"},
                "reverse_primer": {"seq": "TGG"},
                "extra_col": "my target description",
                "special_col": "special target",
            },
            {
                "target_name": "target3",
                "forward_primer": {"seq": "TTG"},
                "reverse_primer": {"seq": "ATT"},
                "extra_col": "my target description",
                "special_col": "special target",
            },
        ]
        self.assertEqual(expected_json, target_info)

    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_targets_are_unique"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_unique_target_info"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.summarise_targets_missing_optional_info"
    )
    def test_panel_info_table_to_pmo_with_dict_genome_info(
        self,
        mock_summarise_targets_missing_optional_info,
        mock_check_unique_target_info,
        mock_check_targets_are_unique,
    ):
        """Test panel_info_table_to_pmo with dict genome_info (should be converted to list)"""
        mock_summarise_targets_missing_optional_info.return_value = [], [], []
        result = panel_info_table_to_pmo(
            self.min_target_table, "test_panel", self.genome_info
        )

        # Check that genome_info dict was converted to list
        self.assertIsInstance(result["targeted_genomes"], list)
        self.assertEqual(len(result["targeted_genomes"]), 1)
        self.assertEqual(result["targeted_genomes"][0], self.genome_info)
        self.assertIn("panel_info", result)
        self.assertIn("target_info", result)

    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_targets_are_unique"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_unique_target_info"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.summarise_targets_missing_optional_info"
    )
    def test_panel_info_table_to_pmo_with_list_genome_info(
        self,
        mock_summarise_targets_missing_optional_info,
        mock_check_unique_target_info,
        mock_check_targets_are_unique,
    ):
        """Test panel_info_table_to_pmo with list genome_info"""
        mock_summarise_targets_missing_optional_info.return_value = [], [], []
        genome_info_list = [
            {
                "name": "3D7",
                "url": "somthing.com",
                "genome_version": 1,
                "taxon_id": 1,
            },
            {
                "name": "DD2",
                "url": "somthingelse.com",
                "genome_version": 2,
                "taxon_id": 2,
            },
        ]
        result = panel_info_table_to_pmo(
            self.min_target_table, "test_panel", genome_info_list
        )

        # Check that genome_info list is used directly
        self.assertIsInstance(result["targeted_genomes"], list)
        self.assertEqual(len(result["targeted_genomes"]), 2)
        self.assertEqual(result["targeted_genomes"], genome_info_list)
        self.assertIn("panel_info", result)
        self.assertIn("target_info", result)

    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_targets_are_unique"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_unique_target_info"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.summarise_targets_missing_optional_info"
    )
    def test_create_targets_dict_with_genome_id_col(
        self,
        mock_summarise_targets_missing_optional_info,
        mock_check_unique_target_info,
        mock_check_targets_are_unique,
    ):
        """Test create_targets_dict with genome_id_col parameter"""
        mock_summarise_targets_missing_optional_info.return_value = [], [], []
        target_table = self.min_target_table.copy()
        target_table["target_start"] = [1, 2, 3]
        target_table["target_end"] = [5, 6, 7]
        target_table["insert_start"] = [2, 3, 4]
        target_table["insert_end"] = [3, 4, 5]
        target_table["chrom"] = "chrom1"
        target_table["genome_id"] = [
            0,
            1,
            0,
        ]  # Different genome_ids for different targets

        builder = PMOPanelBuilder(
            target_table,
            "test_panel",
            self.genome_info,
            "target_name",
            "fwd_primer",
            "rev_primer",
            forward_primers_start_col="target_start",
            forward_primers_end_col="insert_start",
            reverse_primers_start_col="insert_end",
            reverse_primers_end_col="target_end",
            insert_start_col="insert_start",
            insert_end_col="insert_end",
            chrom_col="chrom",
        )

        target_info = builder.create_targets_dict(genome_id_col="genome_id")

        # Check that genome_id values come from the column
        self.assertEqual(target_info[0]["insert_location"]["genome_id"], 0)
        self.assertEqual(target_info[0]["forward_primer"]["location"]["genome_id"], 0)
        self.assertEqual(target_info[0]["reverse_primer"]["location"]["genome_id"], 0)

        self.assertEqual(target_info[1]["insert_location"]["genome_id"], 1)
        self.assertEqual(target_info[1]["forward_primer"]["location"]["genome_id"], 1)
        self.assertEqual(target_info[1]["reverse_primer"]["location"]["genome_id"], 1)

        self.assertEqual(target_info[2]["insert_location"]["genome_id"], 0)
        self.assertEqual(target_info[2]["forward_primer"]["location"]["genome_id"], 0)
        self.assertEqual(target_info[2]["reverse_primer"]["location"]["genome_id"], 0)

    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_targets_are_unique"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_unique_target_info"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.summarise_targets_missing_optional_info"
    )
    def test_create_targets_dict_without_genome_id_col(
        self,
        mock_summarise_targets_missing_optional_info,
        mock_check_unique_target_info,
        mock_check_targets_are_unique,
    ):
        """Test create_targets_dict without genome_id_col (should default to 0)"""
        mock_summarise_targets_missing_optional_info.return_value = [], [], []
        target_table = self.min_target_table.copy()
        target_table["target_start"] = [1, 2, 3]
        target_table["target_end"] = [5, 6, 7]
        target_table["insert_start"] = [2, 3, 4]
        target_table["insert_end"] = [3, 4, 5]
        target_table["chrom"] = "chrom1"

        builder = PMOPanelBuilder(
            target_table,
            "test_panel",
            self.genome_info,
            "target_name",
            "fwd_primer",
            "rev_primer",
            forward_primers_start_col="target_start",
            forward_primers_end_col="insert_start",
            reverse_primers_start_col="insert_end",
            reverse_primers_end_col="target_end",
            insert_start_col="insert_start",
            insert_end_col="insert_end",
            chrom_col="chrom",
        )

        target_info = builder.create_targets_dict()

        # Check that genome_id defaults to 0 when genome_id_col is not provided
        for target in target_info:
            self.assertEqual(target["insert_location"]["genome_id"], 0)
            self.assertEqual(target["forward_primer"]["location"]["genome_id"], 0)
            self.assertEqual(target["reverse_primer"]["location"]["genome_id"], 0)

    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_targets_are_unique"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.check_unique_target_info"
    )
    @patch(
        "pmotools.pmo_builder.panel_information_to_pmo.PMOPanelBuilder.summarise_targets_missing_optional_info"
    )
    def test_panel_info_table_to_pmo_with_genome_id_col(
        self,
        mock_summarise_targets_missing_optional_info,
        mock_check_unique_target_info,
        mock_check_targets_are_unique,
    ):
        """Test panel_info_table_to_pmo with genome_id_col parameter"""
        mock_summarise_targets_missing_optional_info.return_value = [], [], []
        target_table = self.min_target_table.copy()
        target_table["target_start"] = [1, 2, 3]
        target_table["target_end"] = [5, 6, 7]
        target_table["insert_start"] = [2, 3, 4]
        target_table["insert_end"] = [3, 4, 5]
        target_table["chrom"] = "chrom1"
        target_table["genome_id"] = [0, 1, 0]

        result = panel_info_table_to_pmo(
            target_table,
            "test_panel",
            self.genome_info,
            forward_primers_start_col="target_start",
            forward_primers_end_col="insert_start",
            reverse_primers_start_col="insert_end",
            reverse_primers_end_col="target_end",
            insert_start_col="insert_start",
            insert_end_col="insert_end",
            chrom_col="chrom",
            genome_id_col="genome_id",
        )

        # Check that genome_id values come from the column in the result
        target_info = result["target_info"]
        self.assertEqual(target_info[0]["insert_location"]["genome_id"], 0)
        self.assertEqual(target_info[1]["insert_location"]["genome_id"], 1)
        self.assertEqual(target_info[2]["insert_location"]["genome_id"], 0)
        self.assertEqual(target_info[0]["forward_primer"]["location"]["genome_id"], 0)
        self.assertEqual(target_info[1]["forward_primer"]["location"]["genome_id"], 1)
        self.assertEqual(target_info[2]["forward_primer"]["location"]["genome_id"], 0)


if __name__ == "__main__":
    unittest.main()
