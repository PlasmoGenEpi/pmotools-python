#!/usr/bin/env python3
import json
import os
import unittest

from pmotools.pmo_engine.pmo_checker import PMOChecker
from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.utils.schema_loader import load_schema


class TestPMOReader(unittest.TestCase):
    def setUp(self):
        self.working_dir = os.path.dirname(os.path.abspath(__file__))

    def test_read_in_pmo(self):
        PMOReader.read_in_pmo(
            os.path.join(
                os.path.dirname(self.working_dir), "data/minimum_pmo_example.json"
            )
        )
        PMOReader.read_in_pmo(
            os.path.join(
                os.path.dirname(self.working_dir), "data/minimum_pmo_example.json.gz"
            )
        )
        PMOReader.read_in_pmo(
            os.path.join(
                os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json"
            )
        )
        PMOReader.read_in_pmo(
            os.path.join(
                os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json.gz"
            )
        )

    def test_read_in_pmos(self):
        PMOReader.read_in_pmos(
            [
                os.path.join(
                    os.path.dirname(self.working_dir), "data/minimum_pmo_example.json"
                ),
                os.path.join(
                    os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json"
                ),
            ]
        )

        PMOReader.read_in_pmos(
            [
                os.path.join(
                    os.path.dirname(self.working_dir),
                    "data/minimum_pmo_example.json.gz",
                ),
                os.path.join(
                    os.path.dirname(self.working_dir),
                    "data/minimum_pmo_example_2.json.gz",
                ),
            ]
        )

        PMOReader.read_in_pmos(
            [
                os.path.join(
                    os.path.dirname(self.working_dir),
                    "data/minimum_pmo_example.json.gz",
                ),
                os.path.join(
                    os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json"
                ),
            ]
        )

    def test_combine_multiple_pmos(self):
        pmo_data_list = PMOReader.read_in_pmos(
            [
                os.path.join(
                    os.path.dirname(self.working_dir), "data/minimum_pmo_example.json"
                ),
                os.path.join(
                    os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json"
                ),
            ]
        )
        combined_pmo = PMOReader.combine_multiple_pmos(pmo_data_list)
        # validate with schema
        pmo_jsonschema_data = load_schema(
            "portable_microhaplotype_object_v1.0.0.schema.json"
        )
        checker = PMOChecker(pmo_jsonschema_data)
        checker.validate_pmo_json(combined_pmo)
        # check against expected
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/combined_pmo_example.json"
            )
        ) as f:
            expected_pmo = json.load(f)
        # remove the pmo_header as the generation date will be new each time it's created
        expected_pmo.pop("pmo_header")
        combined_pmo.pop("pmo_header")

        self.assertEqual(expected_pmo, combined_pmo)

    def test_combine_multiple_pmos_fail_dup_specimen_names(self):
        # the two files below have same specimen_names but have different meta so will fail when trying to combine
        pmo_data_list = PMOReader.read_in_pmos(
            [
                os.path.join(
                    os.path.dirname(self.working_dir),
                    "data/minimum_pmo_example_2_for_spec_dup_testing.json",
                ),
                os.path.join(
                    os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json"
                ),
            ]
        )
        self.assertRaises(Exception, PMOReader.combine_multiple_pmos, pmo_data_list)

    def test_combine_multiple_pmos_fail_dup_library_sample_names(self):
        # the two files below have same library sample names so will fail for duplicated library_sample_names
        pmo_data_list_2 = PMOReader.read_in_pmos(
            [
                os.path.join(
                    os.path.dirname(self.working_dir),
                    "data/minimum_pmo_example_2_for_library_sample_dup_testing.json",
                ),
                os.path.join(
                    os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json"
                ),
            ]
        )
        self.assertRaises(Exception, PMOReader.combine_multiple_pmos, pmo_data_list_2)

    def test_combine_multiple_pmos_fail_for_combine_only_one_file(self):
        # will fail for only having 1 PMO
        pmo_data_list_2 = PMOReader.read_in_pmos(
            [
                os.path.join(
                    os.path.dirname(self.working_dir),
                    "data/minimum_pmo_example_2_for_library_sample_dup_testing.json",
                )
            ]
        )
        self.assertRaises(Exception, PMOReader.combine_multiple_pmos, pmo_data_list_2)


if __name__ == "__main__":
    unittest.main()
