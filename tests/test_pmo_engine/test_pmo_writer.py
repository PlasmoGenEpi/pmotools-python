#!/usr/bin/env python3

import os
import tempfile
import unittest
import json
from pmotools.pmo_engine.pmo_writer import PMOWriter
import hashlib
import gzip


class TestPMOWriter(unittest.TestCase):
    def setUp(self):
        self.working_dir = os.path.dirname(os.path.abspath(__file__))
        self.test_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.test_dir.cleanup()

    def test_add_pmo_extension_as_needed(self):
        # the output fnps to test
        output_fnp_1 = "out_pmo"
        output_fnp_2 = "out_pmo.json"
        output_fnp_3 = "out_pmo.json.gz"

        # output_fnp_1
        self.assertEqual(
            "out_pmo.json", PMOWriter.add_pmo_extension_as_needed(output_fnp_1, False)
        )
        self.assertEqual(
            "out_pmo.json.gz", PMOWriter.add_pmo_extension_as_needed(output_fnp_1, True)
        )
        # output_fnp_2
        self.assertEqual(
            "out_pmo.json", PMOWriter.add_pmo_extension_as_needed(output_fnp_2, False)
        )
        self.assertEqual(
            "out_pmo.json.gz", PMOWriter.add_pmo_extension_as_needed(output_fnp_2, True)
        )
        # output_fnp_3
        self.assertEqual(
            "out_pmo.json.gz", PMOWriter.add_pmo_extension_as_needed(output_fnp_3, True)
        )

    def test_write_out_pmo(self):
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/minimum_pmo_example.json"
            )
        ) as f:
            pmo_data = json.load(f)
        output_fnp = os.path.join(self.test_dir.name, "out_pmo.json")
        PMOWriter.write_out_pmo(pmo_data, output_fnp, True)
        with open(output_fnp, "rb") as file_to_check:
            md5_returned = hashlib.md5(file_to_check.read()).hexdigest()
        self.assertEqual("f56b922855f471346376e6d928894e4d", md5_returned)

    def test_write_out_pmo_gzip(self):
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/minimum_pmo_example.json"
            )
        ) as f:
            pmo_data = json.load(f)
        output_fnp = os.path.join(self.test_dir.name, "out_pmo.json.gz")
        PMOWriter.write_out_pmo(pmo_data, output_fnp, True)
        hash_md5 = hashlib.md5()
        with gzip.open(output_fnp, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        self.assertEqual("f56b922855f471346376e6d928894e4d", hash_md5.hexdigest())

    def test_write_out_pmo_fail_overwrite(self):
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/minimum_pmo_example.json"
            )
        ) as f:
            pmo_data = json.load(f)
        output_fnp = os.path.join(self.test_dir.name, "out_pmo.json")
        f = open(output_fnp, mode="w")
        f.close()
        self.assertRaises(
            Exception, PMOWriter.write_out_pmo, pmo_data, output_fnp, False
        )


if __name__ == "__main__":
    unittest.main()
