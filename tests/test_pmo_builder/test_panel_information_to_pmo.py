import unittest
import pandas as pd
from unittest.mock import patch

from pmotools.pmo_builder.panel_information_to_pmo import (
    PMOPanelBuilder,
    check_genome_info,
    merge_panel_info_dicts,
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

    def test_check_genome_info_fails_with_string(self):
        with self.assertRaises(TypeError) as context:
            check_genome_info("genome_info")
        self.assertEqual(
            str(context.exception), "genome_info must be a dict, but got str"
        )

    def test_check_genome_info_fails_with_missing_key(self):
        with self.assertRaises(ValueError):
            check_genome_info({"name": "3D7", "taxon_id": 1})

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

    def test_merge_panel_info_dicts_with_overlap_no_reaction(self):
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


if __name__ == "__main__":
    unittest.main()
