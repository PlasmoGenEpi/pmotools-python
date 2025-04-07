#!/usr/bin/env python3

import os
import unittest
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



if __name__ == "__main__":
    unittest.main()
