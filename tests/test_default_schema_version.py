#!/usr/bin/env python3
import inspect
import os
import unittest

from pmotools import __schema_version__
from pmotools.pmo_engine.pmo_exporter import PMOExporter
from pmotools.scripts.pmo_to_tables import extract_allele_table
from pmotools.scripts.pmo_utils import validate_pmo


class TestDefaultSchemaVersion(unittest.TestCase):
    """The default schema must follow the schema version, not the package version."""

    def setUp(self):
        self.expected_schema_fname = (
            f"portable_microhaplotype_object_v{__schema_version__}.schema.json"
        )

    def assert_default_schema(self, fnp):
        self.assertEqual(os.path.basename(fnp), self.expected_schema_fname)
        self.assertTrue(os.path.exists(fnp), f"{fnp} does not exist")

    def test_extract_alleles_per_sample_table_default_schema(self):
        default = (
            inspect.signature(PMOExporter.extract_alleles_per_sample_table)
            .parameters["jsonschema_fnp"]
            .default
        )
        self.assert_default_schema(default)

    def test_extract_allele_table_cli_default_schema(self):
        self.assert_default_schema(
            extract_allele_table.get_parser().get_default("jsonschema")
        )

    def test_validate_pmo_cli_default_schema_version(self):
        self.assertEqual(
            validate_pmo.get_parser_validate_pmo().get_default("jsonschema_version"),
            __schema_version__,
        )


if __name__ == "__main__":
    unittest.main()
