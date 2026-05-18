import unittest
from unittest.mock import patch
from datetime import date


from pmotools.pmo_builder.merge_to_pmo import (
    _report_missing_IDs,
    _generate_pmo_header,
    _replace_key_with_id,
    _make_lookup,
    merge_to_pmo,
)


class TestMergeToPMO(unittest.TestCase):
    def setUp(self):
        self.ref_list = [{"name": "name1"}, {"name": "name2"}, {"name": "name3"}]
        self.pmo_header_v1_0_0 = {
            "pmo_version": "1.0.0",
            "creation_date": "2025-07-22",
            "generation_method": {
                "program_name": "pmotools-python",
                "program_version": "1.0.0",
            },
        }
        self.pmo_header_v1_1_0 = {
            "pmo_version": "1.1.0",
            "creation_date": "2025-07-22",
            "generation_method": {
                "program_name": "pmotools-python",
                "program_version": "1.1.0",
            },
        }

    def test_report_missing_IDs_passes(self):
        _report_missing_IDs(
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
        )

    def test_report_missing_IDs_fails_correctly(self):
        with self.assertRaises(ValueError) as context:
            _report_missing_IDs(
                ["something"],
                [],
                [],
                [],
                ["something else", "something else2"],
                [],
                [],
                [],
                [],
                [],
            )
        self.assertEqual(
            "The following fields were found in one table and not another:\nProject names in Specimen Info not in Project Info: ['something']\nTarget names in Representative Microhaplotypes not in Target Info: ['something else', 'something else2']\n",
            str(context.exception),
        )

    def test_merge_to_pmo_with_read_counts_by_stage(self):
        """Test merge_to_pmo with read_counts_by_stage_info."""
        # Sample data
        specimen_info = [{"specimen_name": "spec1", "project_name": "proj1"}]
        library_sample_info = [
            {
                "library_sample_name": "lib1",
                "specimen_name": "spec1",
                "sequencing_info_name": "seq1",
                "panel_name": "panel1",
            }
        ]
        sequencing_info = [{"sequencing_info_name": "seq1"}]
        panel_info = {
            "panel_info": [{"panel_name": "panel1"}],
            "target_info": [{"target_name": "target1"}],
        }
        mhap_info = {
            "representative_microhaplotypes": {
                "targets": [{"target_name": "target1", "microhaplotypes": []}]
            },
            "detected_microhaplotypes": [],
        }
        bioinfo_method_info = []
        bioinfo_run_info = [{"bioinformatics_run_name": "run1"}]
        project_info = [{"project_name": "proj1"}]

        # Read counts by stage data
        read_counts_by_stage_info = [
            {
                "bioinformatics_run_name": "run1",
                "read_counts_by_library_sample_by_stage": [
                    {
                        "library_sample_name": "lib1",
                        "total_raw_count": 1000,
                        "read_counts_for_targets": [
                            {
                                "target_name": "target1",
                                "stages": [
                                    {"stage": "demultiplexed", "read_count": 100}
                                ],
                            }
                        ],
                    }
                ],
            }
        ]

        result = merge_to_pmo(
            specimen_info=specimen_info,
            library_sample_info=library_sample_info,
            sequencing_info=sequencing_info,
            panel_and_target_info=panel_info,
            mhap_info=mhap_info,
            bioinfo_method_info=bioinfo_method_info,
            bioinfo_run_info=bioinfo_run_info,
            project_info=project_info,
            read_counts_by_stage_info=read_counts_by_stage_info,
        )

        # Check that read_counts_by_stage is included
        self.assertIn("read_counts_by_stage", result)
        self.assertEqual(len(result["read_counts_by_stage"]), 1)

        # Check that names have been replaced with IDs
        read_counts_run = result["read_counts_by_stage"][0]
        self.assertIn("bioinformatics_run_id", read_counts_run)
        self.assertNotIn("bioinformatics_run_name", read_counts_run)

        # Check library sample names have been replaced
        library_sample = read_counts_run["read_counts_by_library_sample_by_stage"][0]
        self.assertIn("library_sample_id", library_sample)
        self.assertNotIn("library_sample_name", library_sample)

    # need to test raise on not giving spec or libr and giving more than 1 panel

    def test_merge_to_pmo_minimum_info(self):
        """Test merge_to_pmo with only target and panel info and mhaps detected."""

        panel_info = {
            "panel_info": [{"panel_name": "panel1"}],
            "target_info": [{"target_name": "target1"}],
        }
        mhap_info = {
            "representative_microhaplotypes": {
                "targets": [{"target_name": "target1", "microhaplotypes": []}]
            },
            "detected_microhaplotypes": [
                {
                    "library_samples": [
                        {"library_sample_name": "samp_A", "target_results": []},
                        {"library_sample_name": "samp_B", "target_results": []},
                    ]
                }
            ],
        }

        result = merge_to_pmo(
            panel_and_target_info=panel_info,
            mhap_info=mhap_info,
        )

        # test that library_sample_info and specimen_info are created and contain the library_sample names from detected_microhaplotypes
        self.assertIn("library_sample_info", result)
        self.assertIn("specimen_info", result)
        self.assertEqual(
            sorted(
                [
                    result["library_sample_info"][0]["library_sample_name"],
                    result["library_sample_info"][1]["library_sample_name"],
                ]
            ),
            sorted(["samp_A", "samp_B"]),
        )
        self.assertEqual(
            sorted(
                [
                    result["specimen_info"][0]["specimen_name"],
                    result["specimen_info"][1]["specimen_name"],
                ]
            ),
            sorted(["samp_A", "samp_B"]),
        )

    def test_merge_to_pmo_minimum_info_plus_library_sample(self):
        """Test merge_to_pmo with only target and panel info and mhaps detected."""

        panel_info = {
            "panel_info": [{"panel_name": "panel1"}, {"panel_name": "panel2"}],
            "target_info": [{"target_name": "target1"}],
        }
        library_sample_info = [
            {"library_sample_name": "samp_A", "panel_name": "panel1"},
            {"library_sample_name": "samp_B", "panel_name": "panel2"},
        ]
        mhap_info = {
            "representative_microhaplotypes": {
                "targets": [{"target_name": "target1", "microhaplotypes": []}]
            },
            "detected_microhaplotypes": [
                {
                    "library_samples": [
                        {"library_sample_name": "samp_A", "target_results": []},
                        {"library_sample_name": "samp_B", "target_results": []},
                    ]
                }
            ],
        }

        result = merge_to_pmo(
            panel_and_target_info=panel_info,
            mhap_info=mhap_info,
            library_sample_info=library_sample_info,
        )

        # test that specimen_info was created and contain the library_sample names from library_sample_info
        self.assertIn("specimen_info", result)
        self.assertEqual(
            sorted(
                [
                    result["library_sample_info"][0]["library_sample_name"],
                    result["library_sample_info"][1]["library_sample_name"],
                ]
            ),
            sorted(["samp_A", "samp_B"]),
        )
        self.assertEqual(
            sorted(
                [
                    result["specimen_info"][0]["specimen_name"],
                    result["specimen_info"][1]["specimen_name"],
                ]
            ),
            sorted(["samp_A", "samp_B"]),
        )

    def test_merge_to_pmo_minimum_info_raise_spec_but_no_lib(self):
        panel_info = {
            "panel_info": [{"panel_name": "panel1"}, {"panel_name": "panel2"}],
            "target_info": [{"target_name": "target1"}],
        }
        specimen_info = [{"specimen_name": "samp_A"}, {"specimen_name": "samp_B"}]
        mhap_info = {
            "representative_microhaplotypes": {
                "targets": [{"target_name": "target1", "microhaplotypes": []}]
            },
            "detected_microhaplotypes": [
                {
                    "library_samples": [
                        {"library_sample_name": "samp_A", "target_results": []},
                        {"library_sample_name": "samp_B", "target_results": []},
                    ]
                }
            ],
        }
        with self.assertRaises(ValueError) as context:
            merge_to_pmo(
                panel_and_target_info=panel_info,
                mhap_info=mhap_info,
                specimen_info=specimen_info,
            )
        self.assertIn(
            "library_sample_info must be provided if specimen_info is provided",
            str(context.exception),
        )

    def test_merge_to_pmo_minimum_info_raise_multi_panel_no_lib(self):
        panel_info = {
            "panel_info": [{"panel_name": "panel1"}, {"panel_name": "panel2"}],
            "target_info": [{"target_name": "target1"}],
        }
        mhap_info = {
            "representative_microhaplotypes": {
                "targets": [{"target_name": "target1", "microhaplotypes": []}]
            },
            "detected_microhaplotypes": [
                {
                    "library_samples": [
                        {"library_sample_name": "samp_A", "target_results": []},
                        {"library_sample_name": "samp_B", "target_results": []},
                    ]
                }
            ],
        }
        with self.assertRaises(Exception) as context:
            merge_to_pmo(panel_and_target_info=panel_info, mhap_info=mhap_info)
        self.assertIn(
            "If providing only microhaplotypes, but not specimen_info or library_sample_info have to have only 1 panel to default to, found",
            str(context.exception),
        )

    @patch("pmotools.pmo_builder.merge_to_pmo.date")
    def test_generate_pmo_header(self, mock_date):
        mock_date.today.return_value = date(2025, 7, 22)
        mock_date.side_effect = lambda *args, **kwargs: date(*args, **kwargs)
        actual = _generate_pmo_header("1.0.0")
        # expected = {'pmo_version': '1.0.0', 'creation_date': '2025-07-22', 'generation_method': {
        #     'program_name': 'pmotools-python', 'program_version': '1.0.0'}}
        self.assertEqual(actual, self.pmo_header_v1_0_0)

    def test_replace_key_with_id(self):
        test_target_list = [
            {"name": "name1"},
            {"name": "name2"},
            {"name": "name1"},
            {"name": "name2"},
            {"name": "name3"},
        ]
        actual = _replace_key_with_id(test_target_list, self.ref_list, "name", "ID")
        expected = [{"ID": 0}, {"ID": 1}, {"ID": 0}, {"ID": 1}, {"ID": 2}]
        self.assertEqual(actual, [])
        self.assertEqual(test_target_list, expected)

    def test_replace_key_with_id_with_lookup(self):
        lookup = {"name1": 5, "name2": 7, "name3": 49}
        test_target_list = [
            {"name": "name1"},
            {"name": "name2"},
            {"name": "name1"},
            {"name": "name2"},
            {"name": "name3"},
        ]
        actual = _replace_key_with_id(
            test_target_list, self.ref_list, "name", "ID", lookup
        )
        expected = [{"ID": 5}, {"ID": 7}, {"ID": 5}, {"ID": 7}, {"ID": 49}]
        self.assertEqual(actual, [])
        self.assertEqual(test_target_list, expected)

    def test_replace_key_with_id_with_missing(self):
        lookup = {"name1": 5, "name2": 7}
        test_target_list = [
            {"name": "name1"},
            {"name": "name2"},
            {"name": "name1"},
            {"name": "name2"},
            {"name": "name3"},
        ]
        actual = _replace_key_with_id(
            test_target_list, self.ref_list, "name", "ID", lookup
        )
        expected = [{"ID": 5}, {"ID": 7}, {"ID": 5}, {"ID": 7}, {"ID": None}]
        self.assertEqual(actual, ["name3"])
        self.assertEqual(test_target_list, expected)

    def test_make_lookup(self):
        actual = _make_lookup(self.ref_list, "name")
        expected = {"name1": 0, "name2": 1, "name3": 2}
        self.assertEqual(expected, actual)

    @patch("pmotools.pmo_builder.merge_to_pmo._replace_names_with_IDs")
    @patch("pmotools.pmo_builder.merge_to_pmo._generate_pmo_header")
    def test_merge_to_pmo(self, mock_generate_pmo_header, _):
        mock_generate_pmo_header.return_value = self.pmo_header_v1_0_0
        actual = merge_to_pmo(
            specimen_info=[{"specimens": "specinfo"}],
            library_sample_info=[{"library_samples": "library_samples"}],
            sequencing_info=[{"sequencing": "sequencing"}],
            panel_and_target_info={
                "panel_info": ["panels"],
                "target_info": ["targets"],
            },
            mhap_info={
                "representative_microhaplotypes": ["mhap_seqs"],
                "detected_microhaplotypes": ["mhaps for sample"],
            },
            bioinfo_method_info=[{"bioinfo_methods": "bioinfo_methods"}],
            bioinfo_run_info=[{"bioinfo_runs": "bioinfo_runs"}],
            project_info=[{"projects": "projects"}],
        )

        expected = {
            "pmo_header": self.pmo_header_v1_0_0,
            "library_sample_info": [{"library_samples": "library_samples"}],
            "specimen_info": [{"specimens": "specinfo"}],
            "sequencing_info": [{"sequencing": "sequencing"}],
            "bioinformatics_methods_info": [{"bioinfo_methods": "bioinfo_methods"}],
            "bioinformatics_run_info": [{"bioinfo_runs": "bioinfo_runs"}],
            "project_info": [{"projects": "projects"}],
            "panel_info": ["panels"],
            "target_info": ["targets"],
            "representative_microhaplotypes": ["mhap_seqs"],
            "detected_microhaplotypes": ["mhaps for sample"],
        }
        self.assertEqual(expected, actual)


if __name__ == "__main__":
    unittest.main()
