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
        self.pmo_header = {
            "pmo_version": "0.1.0",
            "creation_date": "2025-07-22",
            "generation_method": {
                "program_name": "pmotools-python",
                "program_version": "0.1.0",
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
            )
        self.assertEqual(
            "The following fields were found in one table and not another:\nProject names in Specimen Info not in Project Info: ['something']\nTarget names in Representative Microhaplotypes not in Target Info: ['something else', 'something else2']\n",
            str(context.exception),
        )

    @patch("pmotools.pmo_builder.merge_to_pmo.date")
    def test_generate_pmo_header(self, mock_date):
        mock_date.today.return_value = date(2025, 7, 22)
        mock_date.side_effect = lambda *args, **kwargs: date(*args, **kwargs)
        actual = _generate_pmo_header()
        # expected = {'pmo_version': '1.0.0', 'creation_date': '2025-07-22', 'generation_method': {
        #     'program_name': 'pmotools-python', 'program_version': '1.0.0'}}
        self.assertEqual(actual, self.pmo_header)

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
        mock_generate_pmo_header.return_value = self.pmo_header
        actual = merge_to_pmo(
            [{"specimens": "specinfo"}],
            [{"library_samples": "library_samples"}],
            [{"sequencing": "sequencing"}],
            {"panel_info": ["panels"], "target_info": ["targets"]},
            {
                "representative_microhaplotypes": ["mhap_seqs"],
                "detected_microhaplotypes": ["mhaps for sample"],
            },
            [{"bioinfo_methods": "bioinfo_methods"}],
            [{"bioinfo_runs": "bioinfo_runs"}],
            [{"projects": "projects"}],
        )

        expected = {
            "pmo_header": self.pmo_header,
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
