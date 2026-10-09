#!/usr/bin/env python3
import json
import os
import unittest

from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.pmo_engine.pmo_checker import PMOChecker, load_schema


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

    def test_combine_multiple_pmos_v1_1_0(self):
        pmo_data_list = PMOReader.read_in_pmos(
            [
                os.path.join(
                    os.path.dirname(self.working_dir),
                    "data/minimum_fields_pmo_example1.json.gz",
                ),
                os.path.join(
                    os.path.dirname(self.working_dir),
                    "data/minimum_fields_pmo_example2.json.gz",
                ),
            ]
        )
        combined_pmo = PMOReader.combine_multiple_pmos(pmo_data_list)
        # validate with schema
        pmo_jsonschema_data = load_schema(
            "portable_microhaplotype_object_v1.1.0.schema.json"
        )
        checker = PMOChecker(pmo_jsonschema_data)
        checker.validate_pmo_json(combined_pmo)
        # check against expected
        with open(
            os.path.join(
                os.path.dirname(self.working_dir),
                "data/combined_pmo_minimum_fields_example.json",
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

    def test_combine_multiple_pmos_remaps_new_rep_target_id(self):
        # a representative-microhaplotype target newly added from a
        # non-first PMO must have its target_id remapped to the combined
        # target_info (previously it kept its local target_id), this test
        # ensures this is tested as it was previously untested (2026-06-29)
        import pandas as pd
        from pmotools.pmo_builder.mhap_table_to_pmo import mhap_table_to_pmo
        from pmotools.pmo_builder.panel_information_to_pmo import (
            panel_info_table_to_pmo,
        )
        from pmotools.pmo_builder.merge_to_pmo import merge_to_pmo

        def build(libs, targets, seqs, panel_name, panel_targets):
            mhap_info = mhap_table_to_pmo(
                pd.DataFrame(
                    {
                        "library_sample_name": libs,
                        "target_name": targets,
                        "seq": seqs,
                        "reads": [10] * len(libs),
                    }
                )
            )
            panel_info = panel_info_table_to_pmo(
                pd.DataFrame(
                    {
                        "target_name": panel_targets,
                        "fwd_primer": [
                            "A" * (i + 4) for i in range(len(panel_targets))
                        ],
                        "rev_primer": [
                            "C" * (i + 4) for i in range(len(panel_targets))
                        ],
                    }
                ),
                panel_name,
            )
            return merge_to_pmo(mhap_info=mhap_info, panel_target_info=panel_info)

        pmo_a = build(
            ["S1", "S2"], ["t1", "t2"], ["AAA", "GGG"], "panelA", ["t1", "t2"]
        )
        # pmo_b introduces a brand-new target t3
        pmo_b = build(
            ["S3", "S4"], ["t1", "t3"], ["AAA", "CCC"], "panelB", ["t1", "t3"]
        )

        combined = PMOReader.combine_multiple_pmos([pmo_a, pmo_b])
        self.assertEqual(
            [t["target_name"] for t in combined["target_info"]],
            ["t1", "t2", "t3"],
        )
        # every rep target's target_id resolves to the right combined target
        for rep in combined["representative_microhaplotypes"]["targets"]:
            self.assertLess(rep["target_id"], len(combined["target_info"]))
        t3_reps = [
            rep
            for rep in combined["representative_microhaplotypes"]["targets"]
            if combined["target_info"][rep["target_id"]]["target_name"] == "t3"
        ]
        self.assertEqual(len(t3_reps), 1)
        # validate against schema (minimal builder PMO targets v1.1.0)
        checker = PMOChecker(
            load_schema("portable_microhaplotype_object_v1.1.0.schema.json")
        )
        checker.validate_pmo_json(combined)

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
