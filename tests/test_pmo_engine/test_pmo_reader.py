#!/usr/bin/env python3
import json
import os
import unittest

from pmotools.pmo_engine.pmo_checker import PMOChecker
from pmotools.pmo_engine.pmo_reader import PMOReader


class TestPMOReader(unittest.TestCase):
    def setUp(self):
        self.working_dir = os.path.dirname(os.path.abspath(__file__))

    def test_read_in_pmo(self):
        pmo_data = PMOReader.read_in_pmo(os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example.json"))
        pmo_data_gz = PMOReader.read_in_pmo(os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example.json.gz"))
        pmo_data_2 = PMOReader.read_in_pmo(os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json"))
        pmo_data_2_gz = PMOReader.read_in_pmo(os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json.gz"))

    def test_read_in_pmos(self):
        pmo_data_list = PMOReader.read_in_pmos([os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example.json"),
                                            os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json")])

        pmo_data_list_gz = PMOReader.read_in_pmos([os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example.json.gz"),
                                            os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json.gz")])

        pmo_data_list_mix = PMOReader.read_in_pmos([os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example.json.gz"),
                                            os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json")])

    def test_combine_multiple_pmos(self):
        pmo_data_list = PMOReader.read_in_pmos([
            os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example.json"),
            os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json")
        ])
        combined_pmo = PMOReader.combine_multiple_pmos(pmo_data_list)
        # validate with schema
        pmo_jsonschema_fnp = os.path.join(os.path.dirname(os.path.dirname(self.working_dir)),
                                               "etc/portable_microhaplotype_object.schema.json")
        with open(pmo_jsonschema_fnp) as f:
            pmo_jsonschema_data = json.load(f)
        checker = PMOChecker(pmo_jsonschema_data)
        checker.validate_pmo_json(combined_pmo)
        # check against expected
        with open(os.path.join(os.path.dirname(self.working_dir), "data/combined_pmo_example.json")) as f:
            expected_pmo = json.load(f)
        self.assertEqual(expected_pmo, combined_pmo)



if __name__ == "__main__":
    unittest.main()
