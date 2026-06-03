#!/usr/bin/env python3
import gzip
import os
import unittest
import json
from jsonschema import ValidationError
from pmotools.pmo_engine.pmo_checker import PMOChecker, load_schema


class TestPMOChecker(unittest.TestCase):
    def setUp(self):
        self.working_dir = os.path.dirname(os.path.abspath(__file__))
        self.pmo_jsonschema_data_v1_0_0 = load_schema(
            "portable_microhaplotype_object_v1.0.0.schema.json"
        )
        self.checker_v1_0_0 = PMOChecker(self.pmo_jsonschema_data_v1_0_0)
        self.pmo_jsonschema_data_v1_1_0 = load_schema(
            "portable_microhaplotype_object_v1.1.0.schema.json"
        )
        self.checker_v1_1_0 = PMOChecker(self.pmo_jsonschema_data_v1_1_0)

        self.pmo_required_base_fields = [
            "library_sample_info",
            "specimen_info",
            "sequencing_info",
            "panel_info",
            "target_info",
            "targeted_genomes",
            "representative_microhaplotypes",
            "bioinformatics_methods_info",
            "bioinformatics_run_info",
            "detected_microhaplotypes",
            "project_info" "pmo_header",
        ]

        self.specimen_required_fields_v1_0_0 = [
            "specimen_name",
            "specimen_taxon_id",
            "host_taxon_id",
            "collection_date",
            "collection_country",
            "project_id",
        ]

    def test_pmo_checker_check_for_required_base_fields(self):
        pmo_test_object_v1_0_0 = {
            "library_sample_info": [],
            "specimen_info": [],
            "sequencing_info": [],
            "panel_info": [],
            "target_info": [],
            "targeted_genomes": [],
            "representative_microhaplotypes": {},
            "bioinformatics_methods_info": [],
            "bioinformatics_run_info": [],
            "detected_microhaplotypes": [],
            "project_info": [],
            "pmo_header": {},
        }
        self.checker_v1_0_0.check_for_required_base_fields(pmo_test_object_v1_0_0)

    def test_pmo_checker_check_for_required_base_fields_fail(self):
        pmo_test_object = {
            "library_sample_info": [],
            "specimen_info": [],
            "sequencing_info": [],
            "panel_info": [],
            "target_info": [],
            "targeted_genome": [],
            "microhaplotypes_info": {},
            "bioinformatics_methods_infos": [],
            "bioinformatics_run_info": [],
            "microhaplotypes_detected": [],
            "pmo_headers": {},
        }
        self.assertRaises(
            Exception,
            self.checker_v1_0_0.check_for_required_base_fields,
            pmo_test_object,
        )

    def test_pmo_checker_get_required_fields_for_pmo_class(self):
        self.assertEqual(
            self.checker_v1_0_0.get_required_fields_for_pmo_class("SpecimenInfo"),
            self.specimen_required_fields_v1_0_0,
        )

    def test_pmo_checker_validate_pmo_json(self):
        # minimum_pmo_example.json and minimum_pmo_example_2.json are version 1.0.0 which required more fields and
        # version 1.1.0 requires less fields but still has the same fields defined so should still pass that validation
        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/minimum_pmo_example.json"
            )
        ) as f:
            pmo_data = json.load(f)
        self.checker_v1_0_0.validate_pmo_json(pmo_data)
        self.checker_v1_1_0.validate_pmo_json(pmo_data)

        with open(
            os.path.join(
                os.path.dirname(self.working_dir), "data/minimum_pmo_example_2.json"
            )
        ) as f:
            pmo_data2 = json.load(f)
        self.checker_v1_0_0.validate_pmo_json(pmo_data2)
        self.checker_v1_1_0.validate_pmo_json(pmo_data2)

        with gzip.open(
            os.path.join(
                os.path.dirname(self.working_dir),
                "data/minimum_fields_pmo_example1.json.gz",
            ),
            "rt",
        ) as f:
            pmo_data3 = json.load(f)
        self.checker_v1_1_0.validate_pmo_json(pmo_data3)

        with gzip.open(
            os.path.join(
                os.path.dirname(self.working_dir),
                "data/minimum_fields_pmo_example2.json.gz",
            ),
            "rt",
        ) as f:
            pmo_data4 = json.load(f)
        self.checker_v1_1_0.validate_pmo_json(pmo_data4)

    def test_pmo_checker_validate_pmo_json_fail(self):
        with open(
            os.path.join(
                os.path.dirname(self.working_dir),
                "data/minimum_pmo_example_bad_format.json",
            )
        ) as f:
            pmo_data = json.load(f)
        self.assertRaises(
            ValidationError, self.checker_v1_0_0.validate_pmo_json, pmo_data
        )

        with gzip.open(
            os.path.join(
                os.path.dirname(self.working_dir),
                "data/minimum_fields_pmo_example2.json.gz",
            ),
            "rt",
        ) as f:
            pmo_data4 = json.load(f)
        self.assertRaises(
            ValidationError, self.checker_v1_0_0.validate_pmo_json, pmo_data4
        )


if __name__ == "__main__":
    unittest.main()
