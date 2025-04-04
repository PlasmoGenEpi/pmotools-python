#!/usr/bin/env python3

import os
import unittest
import json
from jsonschema import ValidationError
from pmotools.pmo_engine.pmo_checker import PMOChecker


class TestPMOChecker(unittest.TestCase):
    def setUp(self):
        self.working_dir = os.path.dirname(os.path.abspath(__file__))
        self.pmo_jsonschema_fnp = os.path.join(os.path.dirname(os.path.dirname(self.working_dir)), "etc/portable_microhaplotype_object.schema.json")
        
        with open(self.pmo_jsonschema_fnp) as f: 
            self.pmo_jsonschema_data = json.load(f)
        self.checker = PMOChecker(self.pmo_jsonschema_data)
        self.pmo_required_base_fields = ["experiment_info",
        "specimen_info",
        "sequencing_info",
        "panel_info",
        "target_info",
        "targeted_genomes",
        "microhaplotypes_info",
        "bioinformatics_methods_info",
        "bioinformatics_run_info",
        "microhaplotypes_detected",
        "pmo_header"]

        self.specimen_required_fields = [
                "specimen_name",
                "specimen_taxon_id",
                "host_taxon_id",
                "collection_date",
                "collection_country",
                "project_name"
            ]

    def test_pmo_checker_check_for_required_base_fields(self):

        pmo_test_object = {
            "experiment_info":[],
            "specimen_info":[],
            "sequencing_info":[],
            "panel_info":[],
            "target_info":[],
            "targeted_genomes":[],
            "microhaplotypes_info":{},
            "bioinformatics_methods_info":[],
            "bioinformatics_run_info":[],
            "microhaplotypes_detected":[],
            "pmo_header":{}
        }
        self.checker.check_for_required_base_fields(pmo_test_object)

    def test_pmo_checker_check_for_required_base_fields_fail(self):
        pmo_test_object = {
            "experiment_info":[],
            "specimen_info":[],
            "sequencing_info":[],
            "panel_info":[],
            "target_info":[],
            "targeted_genome":[],
            "microhaplotypes_info":{},
            "bioinformatics_methods_infos":[],
            "bioinformatics_run_info":[],
            "microhaplotypes_detected":[],
            "pmo_headers":{}
        }
        self.assertRaises(Exception, self.checker.check_for_required_base_fields, pmo_test_object)

    def test_pmo_checker_get_required_fields_for_pmo_class(self):
        self.assertEqual(self.checker.get_required_fields_for_pmo_class("SpecimenInfo"), self.specimen_required_fields)

    def test_pmo_checker_validate_pmo_json(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example.json")) as f:
            pmo_data = json.load(f)
        self.checker.validate_pmo_json(pmo_data)

    def test_pmo_checker_validate_pmo_json_fail(self):
        with open(os.path.join(os.path.dirname(self.working_dir), "data/minimum_pmo_example_bad_format.json")) as f:
            pmo_data = json.load(f)
        self.assertRaises(ValidationError, self.checker.validate_pmo_json, pmo_data)

if __name__ == "__main__":
    unittest.main()
