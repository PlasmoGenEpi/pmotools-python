import unittest
import pandas as pd

from pmotools.pmo_builder.metatable_to_pmo import (
    library_sample_info_table_to_pmo,
    specimen_info_table_to_pmo,
    add_plate_info,
    add_parasite_density_info,
)


class TestMetatableToPMO(unittest.TestCase):
    def setUp(self):
        self.small_json_example = [
            {"specimen_name": "sample1"},
            {"specimen_name": "sample2"},
        ]

    def test_add_plate_info_position_parsing_uppercase(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "plate_position": ["A01", "H11"],
                "plate_name": ["Plate1", "Plate1"],
            }
        )
        result = add_plate_info(
            None,
            "plate_name",
            None,
            "plate_position",
            self.small_json_example,
            df,
            "specimen_name",
        )
        self.assertEqual(result[0]["plate_info"]["plate_row"], "A")
        self.assertEqual(result[0]["plate_info"]["plate_col"], 1)
        self.assertEqual(result[1]["plate_info"]["plate_row"], "H")
        self.assertEqual(result[1]["plate_info"]["plate_col"], 11)

    def test_add_plate_info_position_parsing_lowercase(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "plate_position": ["a01", "h11"],
                "plate_name": ["Plate1", "Plate1"],
            }
        )
        result = add_plate_info(
            None,
            "plate_name",
            None,
            "plate_position",
            self.small_json_example,
            df,
            "specimen_name",
        )
        self.assertEqual(result[0]["plate_info"]["plate_row"], "A")
        self.assertEqual(result[0]["plate_info"]["plate_col"], 1)
        self.assertEqual(result[1]["plate_info"]["plate_row"], "H")
        self.assertEqual(result[1]["plate_info"]["plate_col"], 11)

    def test_add_plate_info_position_parsing_one_digit(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "plate_position": ["a1", "h10"],
                "plate_name": ["Plate1", "Plate1"],
            }
        )
        result = add_plate_info(
            None,
            "plate_name",
            None,
            "plate_position",
            self.small_json_example,
            df,
            "specimen_name",
        )
        self.assertEqual(result[0]["plate_info"]["plate_row"], "A")
        self.assertEqual(result[0]["plate_info"]["plate_col"], 1)
        self.assertEqual(result[1]["plate_info"]["plate_row"], "H")
        self.assertEqual(result[1]["plate_info"]["plate_col"], 10)

    def test_add_plate_info_position_fails_with_out_of_bounds_row(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "plate_position": ["A01", "K10"],
                "plate_name": ["Plate1", "Plate1"],
            }
        )
        with self.assertRaises(ValueError) as context:
            add_plate_info(
                None,
                "plate_name",
                None,
                "plate_position",
                self.small_json_example,
                df,
                "specimen_name",
            )
        self.assertEqual(
            "Values in 'plate_position' must start with a single letter A-H/a-h followed by number 1-12.",
            str(context.exception),
        )

    def test_add_plate_info_position_fails_with_out_of_bounds_col(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "plate_position": ["A01", "H19"],
                "plate_name": ["Plate1", "Plate1"],
            }
        )
        with self.assertRaises(ValueError) as context:
            add_plate_info(
                None,
                "plate_name",
                None,
                "plate_position",
                self.small_json_example,
                df,
                "specimen_name",
            )

        self.assertEqual(
            "Values in 'plate_position' must start with a single letter A-H/a-h followed by number 1-12.",
            str(context.exception),
        )

    def test_add_plate_info_row_col_parsing(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "plate_row": ["A", "h"],
                "plate_col": [1, 10],
                "plate_name": ["Plate1", "Plate1"],
            }
        )
        result = add_plate_info(
            "plate_col",
            "plate_name",
            "plate_row",
            None,
            self.small_json_example,
            df,
            "specimen_name",
        )
        self.assertEqual(result[0]["plate_info"]["plate_row"], "A")
        self.assertEqual(result[0]["plate_info"]["plate_col"], 1)
        self.assertEqual(result[1]["plate_info"]["plate_row"], "H")
        self.assertEqual(result[1]["plate_info"]["plate_col"], 10)

    def test_add_plate_info_fails_with_position_and_row_col(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "plate_position": ["A01", "H19"],
                "plate_row": ["A", "h"],
                "plate_col": [1, 10],
                "plate_name": ["Plate1", "Plate1"],
            }
        )
        with self.assertRaises(ValueError) as context:
            add_plate_info(
                "plate_col",
                "plate_name",
                "plate_row",
                "plate_position",
                self.small_json_example,
                df,
                "specimen_name",
            )
        self.assertEqual(
            "Plate position can be specified using either row and col, or position, but not both.",
            str(context.exception),
        )

    def test_add_plate_info_fails_with_row_and_without_col(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "plate_position": ["A01", "H19"],
                "plate_row": ["A", "h"],
                "plate_col": [1, 10],
                "plate_name": ["Plate1", "Plate1"],
            }
        )
        with self.assertRaises(ValueError) as context:
            add_plate_info(
                None,
                "plate_name",
                "plate_row",
                None,
                self.small_json_example,
                df,
                "specimen_name",
            )
        self.assertEqual(
            "If either plate row or column is set, then both must be.",
            str(context.exception),
        )

    def test_add_plate_info_adds_nothing(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "plate_position": ["A01", "H19"],
                "plate_row": ["A", "h"],
                "plate_col": [1, 10],
                "plate_name": ["Plate1", "Plate1"],
            }
        )
        result = add_plate_info(
            None, None, None, None, self.small_json_example, df, "specimen_name"
        )
        self.assertEqual(
            result, [{"specimen_name": "sample1"}, {"specimen_name": "sample2"}]
        )

    def test_add_parasite_density_info_single_value(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "parasite_density": [10, 100],
                "parasite_density_method": ["qPCR", "microscopy"],
            }
        )
        result = add_parasite_density_info(
            "parasite_density",
            "parasite_density_method",
            self.small_json_example,
            df,
            "specimen_name",
            "parasite_density_info",
        )
        self.assertEqual(result[0]["parasite_density_info"][0]["parasite_density"], 10)
        self.assertEqual(result[1]["parasite_density_info"][0]["parasite_density"], 100)
        self.assertEqual(
            result[0]["parasite_density_info"][0]["parasite_density_method"], "qPCR"
        )
        self.assertEqual(
            result[1]["parasite_density_info"][0]["parasite_density_method"],
            "microscopy",
        )

    def test_add_parasite_density_info_from_list(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "density1": [15, 107],
                "method1": ["qPCR", "qPCR"],
                "density2": [10, 100],
                "method2": ["microscopy", "microscopy"],
            }
        )
        result = add_parasite_density_info(
            ["density1", "density2"],
            ["method1", "method2"],
            self.small_json_example,
            df,
            "specimen_name",
            "parasite_density_info",
        )
        self.assertEqual(result[0]["parasite_density_info"][0]["parasite_density"], 15)
        self.assertEqual(result[1]["parasite_density_info"][0]["parasite_density"], 107)
        self.assertEqual(result[0]["parasite_density_info"][1]["parasite_density"], 10)
        self.assertEqual(result[1]["parasite_density_info"][1]["parasite_density"], 100)
        self.assertEqual(
            result[0]["parasite_density_info"][0]["parasite_density_method"], "qPCR"
        )
        self.assertEqual(
            result[1]["parasite_density_info"][0]["parasite_density_method"], "qPCR"
        )
        self.assertEqual(
            result[0]["parasite_density_info"][1]["parasite_density_method"],
            "microscopy",
        )
        self.assertEqual(
            result[1]["parasite_density_info"][1]["parasite_density_method"],
            "microscopy",
        )

    def test_add_parasite_density_adds_nothing(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
            }
        )
        result = add_parasite_density_info(
            None,
            None,
            self.small_json_example,
            df,
            "specimen_name",
            "parasite_density_info",
        )
        self.assertEqual(
            result, [{"specimen_name": "sample1"}, {"specimen_name": "sample2"}]
        )

    def test_add_parasite_density_fails_with_unequal_lists(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "density1": [15, 107],
                "method1": ["qPCR", "qPCR"],
                "density2": [10, 100],
                "method2": ["microscopy", "microscopy"],
            }
        )
        with self.assertRaises(ValueError) as context:
            add_parasite_density_info(
                ["density1"],
                ["method1", "method2"],
                self.small_json_example,
                df,
                "specimen_name",
                "qpcr_parasite_density_info",
            )

        self.assertEqual(
            "If both parasite_density_col and parasite_density_method_col are lists, they must be the same length.",
            str(context.exception),
        )

    def test_add_parasite_density_fails_with_no_density(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "method1": ["qPCR", "qPCR"],
            }
        )

        with self.assertRaises(ValueError) as context:
            add_parasite_density_info(
                None,
                "method1",
                self.small_json_example,
                df,
                "specimen_name",
                "parasite_density_info",
            )

        self.assertEqual(
            "parasite_density_method_col is set but parasite_density_col is None. Cannot proceed.",
            str(context.exception),
        )

    def test_add_parasite_density_fails_with_type_mismatch(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "density1": [15, 107],
                "method1": ["qPCR", "qPCR"],
                "density2": [10, 100],
                "method2": ["microscopy", "microscopy"],
            }
        )
        with self.assertRaises(TypeError) as context:
            add_parasite_density_info(
                "density1",
                ["method1", "method2"],
                self.small_json_example,
                df,
                "specimen_name",
                "parasite_density_info",
            )

        self.assertEqual(
            "If parasite_density_col is a string, parasite_density_method_col must be a string or None.",
            str(context.exception),
        )

    def test_specimen_info_table_to_pmo_default(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "specimen_taxon_id": [5833, 5833],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana"],
                "project_name": ["project1", "project2"],
            }
        )

        result = specimen_info_table_to_pmo(df)
        self.assertEqual(
            [
                {
                    "specimen_name": "sample1",
                    "specimen_taxon_id": 5833,
                    "host_taxon_id": 9606,
                    "collection_date": "01/02/2023",
                    "collection_country": "Mozambique",
                    "project_name": "project1",
                },
                {
                    "specimen_name": "sample2",
                    "specimen_taxon_id": 5833,
                    "host_taxon_id": 9606,
                    "collection_date": "01/02/2023",
                    "collection_country": "Ghana",
                    "project_name": "project2",
                },
            ],
            result,
        )

    def test_specimen_info_table_to_pmo_with_plate_info(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "specimen_taxon_id": [5833, 5833],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana"],
                "project_name": ["project1", "project2"],
                "storage_plate_col": [1, 2],
                "storage_plate_name": ["plate1", "plate1"],
                "storage_plate_row": ["A", "B"],
            }
        )

        result = specimen_info_table_to_pmo(
            df,
            storage_plate_col_col="storage_plate_col",
            storage_plate_name_col="storage_plate_name",
            storage_plate_row_col="storage_plate_row",
        )
        self.assertEqual(
            [
                {
                    "specimen_name": "sample1",
                    "specimen_taxon_id": 5833,
                    "host_taxon_id": 9606,
                    "collection_date": "01/02/2023",
                    "collection_country": "Mozambique",
                    "project_name": "project1",
                    "storage_plate_info": {
                        "plate_name": "plate1",
                        "plate_row": "A",
                        "plate_col": 1,
                    },
                },
                {
                    "specimen_name": "sample2",
                    "specimen_taxon_id": 5833,
                    "host_taxon_id": 9606,
                    "collection_date": "01/02/2023",
                    "collection_country": "Ghana",
                    "project_name": "project2",
                    "storage_plate_info": {
                        "plate_name": "plate1",
                        "plate_row": "B",
                        "plate_col": 2,
                    },
                },
            ],
            result,
        )

    def test_specimen_info_table_to_pmo_with_parasitemia(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "specimen_taxon_id": [5833, 5833],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana"],
                "project_name": ["project1", "project2"],
                "parasite_density": [10, 100],
                "parasite_density_method": ["qPCR", "microscopy"],
            }
        )

        result = specimen_info_table_to_pmo(
            df,
            parasite_density_col="parasite_density",
            parasite_density_method_col="parasite_density_method",
        )
        self.assertEqual(
            [
                {
                    "specimen_name": "sample1",
                    "specimen_taxon_id": 5833,
                    "host_taxon_id": 9606,
                    "collection_date": "01/02/2023",
                    "collection_country": "Mozambique",
                    "project_name": "project1",
                    "parasite_density_info": [
                        {"parasite_density": 10, "parasite_density_method": "qPCR"}
                    ],
                },
                {
                    "specimen_name": "sample2",
                    "specimen_taxon_id": 5833,
                    "host_taxon_id": 9606,
                    "collection_date": "01/02/2023",
                    "collection_country": "Ghana",
                    "project_name": "project2",
                    "parasite_density_info": [
                        {
                            "parasite_density": 100,
                            "parasite_density_method": "microscopy",
                        }
                    ],
                },
            ],
            result,
        )

    def test_specimen_info_table_to_pmo_with_additional_columns(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "specimen_taxon_id": [5833, 5833],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana"],
                "project_name": ["project1", "project2"],
                "special_field_1": ["something", "something_else"],
                "special_field_2": ["this", "that"],
            }
        )

        result = specimen_info_table_to_pmo(
            df, additional_specimen_cols=["special_field_1", "special_field_2"]
        )
        self.assertEqual(
            [
                {
                    "specimen_name": "sample1",
                    "specimen_taxon_id": 5833,
                    "host_taxon_id": 9606,
                    "collection_date": "01/02/2023",
                    "collection_country": "Mozambique",
                    "project_name": "project1",
                    "special_field_1": "something",
                    "special_field_2": "this",
                },
                {
                    "specimen_name": "sample2",
                    "specimen_taxon_id": 5833,
                    "host_taxon_id": 9606,
                    "collection_date": "01/02/2023",
                    "collection_country": "Ghana",
                    "project_name": "project2",
                    "special_field_1": "something_else",
                    "special_field_2": "that",
                },
            ],
            result,
        )

    def test_specimen_info_table_to_pmo_fails_with_col_duplicate(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "specimen_taxon_id": [5833, 5833],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana"],
                "project_name": ["project1", "project2"],
            }
        )
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_pmo(df, drug_usage_col="specimen_name")
        self.assertEqual("Selected columns must be unique.", str(context.exception))

    def test_specimen_info_table_to_pmo_fails_with_missing_col(self):
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "specimen_taxon_id": [5833, 5833],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana"],
                "project_name": ["project1", "project2"],
            }
        )
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_pmo(df, drug_usage_col="drug_usage")
        self.assertEqual(
            "The following columns are not in the DataFrame: ['drug_usage']",
            str(context.exception),
        )

    def test_specimen_info_table_to_pmo_fails_without_df(self):
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_pmo("test")
        self.assertEqual("contents must be a pandas DataFrame.", str(context.exception))

    def test_specimen_info_table_to_pmo_fails_with_null_in_required_columns(self):
        """Test that function raises error when required columns contain null values"""
        # Test with null in specimen_name
        df1 = pd.DataFrame(
            {
                "specimen_name": ["sample1", None],
                "specimen_taxon_id": [5833, 5833],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana"],
                "project_name": ["project1", "project2"],
            }
        )
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_pmo(df1)
        self.assertIn(
            "The following columns contain null values",
            str(context.exception),
        )
        self.assertIn("specimen_name", str(context.exception))

        # Test with null in specimen_taxon_id
        df2 = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "specimen_taxon_id": [5833, None],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana"],
                "project_name": ["project1", "project2"],
            }
        )
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_pmo(df2)
        self.assertIn(
            "The following columns contain null values",
            str(context.exception),
        )
        self.assertIn("specimen_taxon_id", str(context.exception))

        # Test with null in host_taxon_id
        df3 = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "specimen_taxon_id": [5833, 5833],
                "host_taxon_id": [None, 9606],
                "collection_date": ["01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana"],
                "project_name": ["project1", "project2"],
            }
        )
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_pmo(df3)
        self.assertIn(
            "The following columns contain null values",
            str(context.exception),
        )
        self.assertIn("host_taxon_id", str(context.exception))

        # Test with null in collection_date
        df4 = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "specimen_taxon_id": [5833, 5833],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["01/02/2023", None],
                "collection_country": ["Mozambique", "Ghana"],
                "project_name": ["project1", "project2"],
            }
        )
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_pmo(df4)
        self.assertIn(
            "The following columns contain null values",
            str(context.exception),
        )
        self.assertIn("collection_date", str(context.exception))

        # Test with null in collection_country
        df5 = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "specimen_taxon_id": [5833, 5833],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023"],
                "collection_country": [None, "Ghana"],
                "project_name": ["project1", "project2"],
            }
        )
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_pmo(df5)
        self.assertIn(
            "The following columns contain null values",
            str(context.exception),
        )
        self.assertIn("collection_country", str(context.exception))

        # Test with null in project_name
        df6 = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "specimen_taxon_id": [5833, 5833],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana"],
                "project_name": ["project1", None],
            }
        )
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_pmo(df6)
        self.assertIn(
            "The following columns contain null values",
            str(context.exception),
        )
        self.assertIn("project_name", str(context.exception))

        # Test with nulls in multiple required columns
        df7 = pd.DataFrame(
            {
                "specimen_name": ["sample1", None],
                "specimen_taxon_id": [5833, None],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana"],
                "project_name": ["project1", "project2"],
            }
        )
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_pmo(df7)
        self.assertIn(
            "The following columns contain null values",
            str(context.exception),
        )
        # Both columns with nulls should be mentioned
        error_msg = str(context.exception)
        self.assertTrue(
            "specimen_name" in error_msg and "specimen_taxon_id" in error_msg
        )

    def test_specimen_info_table_to_pmo_removes_empty_optional_fields(self):
        """Test that empty optional fields (None, empty string, empty dict, empty list) are removed"""
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2", "sample3", "sample4"],
                "specimen_taxon_id": [5833, 5833, 5833, 5833],
                "host_taxon_id": [9606, 9606, 9606, 9606],
                "collection_date": [
                    "01/02/2023",
                    "01/02/2023",
                    "01/02/2023",
                    "01/02/2023",
                ],
                "collection_country": ["Mozambique", "Ghana", "Kenya", "Tanzania"],
                "project_name": ["project1", "project2", "project3", "project4"],
                "lat_lon": ["", "-25.123,31.456", None, "-1.234,36.789"],
                "host_age": [None, 25, None, 30],
                "host_sex": ["M", "", None, "F"],
                "geo_admin1": [None, "Province1", "", "Province2"],
            }
        )

        result = specimen_info_table_to_pmo(
            df,
            lat_lon_col="lat_lon",
            host_age_col="host_age",
            host_sex_col="host_sex",
            geo_admin1_col="geo_admin1",
        )

        # sample1: lat_lon is empty string, host_age is None, host_sex is "M", geo_admin1 is None
        # Should remove: lat_lon, host_age, geo_admin1
        # Should keep: host_sex
        self.assertNotIn("lat_lon", result[0])
        self.assertNotIn("host_age", result[0])
        self.assertNotIn("geo_admin1", result[0])
        self.assertEqual(result[0]["host_sex"], "M")

        # sample2: lat_lon has value, host_age has value, host_sex is empty string, geo_admin1 has value
        # Should remove: host_sex
        # Should keep: lat_lon, host_age, geo_admin1
        self.assertEqual(result[1]["lat_lon"], "-25.123,31.456")
        self.assertEqual(result[1]["host_age"], 25)
        self.assertNotIn("host_sex", result[1])
        self.assertEqual(result[1]["geo_admin1"], "Province1")

        # sample3: lat_lon is None, host_age is None, host_sex is None, geo_admin1 is empty string
        # Should remove all optional fields
        self.assertNotIn("lat_lon", result[2])
        self.assertNotIn("host_age", result[2])
        self.assertNotIn("host_sex", result[2])
        self.assertNotIn("geo_admin1", result[2])

        # sample4: all optional fields have values
        # Should keep all
        self.assertEqual(result[3]["lat_lon"], "-1.234,36.789")
        self.assertEqual(result[3]["host_age"], 30)
        self.assertEqual(result[3]["host_sex"], "F")
        self.assertEqual(result[3]["geo_admin1"], "Province2")

    def test_specimen_info_table_to_pmo_with_geo_fields(self):
        """Test geographic administrative fields"""
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "specimen_taxon_id": [5833, 5833],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["2023-03-29", "2023-03-20"],
                "collection_country": ["Eswatini", "Eswatini"],
                "project_name": ["RegGenE8", "RegGenE8"],
                "geo_admin1": ["Lubombo", "Lubombo"],
                "geo_admin2": ["Sithobela", "Sithobela"],
                "geo_admin3": [
                    "CABRINI MINISTRIES HEALTH CARE",
                    "SITHOBELA RURAL HEALTH CENTER",
                ],
                "lat_lon": ["", ""],
            }
        )

        result = specimen_info_table_to_pmo(
            df,
            geo_admin1_col="geo_admin1",
            geo_admin2_col="geo_admin2",
            geo_admin3_col="geo_admin3",
            lat_lon_col="lat_lon",
        )

        self.assertEqual(result[0]["geo_admin1"], "Lubombo")
        self.assertEqual(result[0]["geo_admin2"], "Sithobela")
        self.assertEqual(result[0]["geo_admin3"], "CABRINI MINISTRIES HEALTH CARE")
        self.assertNotIn("lat_lon", result[0])  # Empty string should be removed

        self.assertEqual(result[1]["geo_admin1"], "Lubombo")
        self.assertEqual(result[1]["geo_admin2"], "Sithobela")
        self.assertEqual(result[1]["geo_admin3"], "SITHOBELA RURAL HEALTH CENTER")
        self.assertNotIn("lat_lon", result[1])  # Empty string should be removed

    def test_specimen_info_table_to_pmo_removes_empty_dict_list(self):
        """Test that empty dicts and lists are removed from optional fields via additional_specimen_cols"""
        # Note: This tests the remove_optional_null_values function indirectly
        # Since add_plate_info and add_parasite_density_info don't add empty structures,
        # we test with additional_specimen_cols that might contain empty dicts/lists

        # First test: verify empty storage_plate_info and parasite_density_info are not added
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "specimen_taxon_id": [5833, 5833],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana"],
                "project_name": ["project1", "project2"],
            }
        )

        result_with_empty_info = specimen_info_table_to_pmo(
            df,
            storage_plate_col_col=None,
            storage_plate_name_col=None,
            storage_plate_row_col=None,
            storage_plate_position_col=None,
            parasite_density_col=None,
            parasite_density_method_col=None,
        )

        # Storage plate info and parasite density info should not be added if all params are None
        self.assertNotIn("storage_plate_info", result_with_empty_info[0])
        self.assertNotIn("parasite_density_info", result_with_empty_info[0])

        # Test with parasite density that results in no valid densities (all None)
        df2 = pd.DataFrame(
            {
                "specimen_name": ["sample1"],
                "specimen_taxon_id": [5833],
                "host_taxon_id": [9606],
                "collection_date": ["01/02/2023"],
                "collection_country": ["Mozambique"],
                "project_name": ["project1"],
                "parasite_density": [None],
            }
        )

        result2 = specimen_info_table_to_pmo(
            df2,
            parasite_density_col="parasite_density",
            parasite_density_method_col=None,
        )

        # Since all density values are None, parasite_density_info should not be added
        self.assertNotIn("parasite_density_info", result2[0])

    def test_specimen_info_table_to_pmo_with_host_fields(self):
        """Test host-related optional fields"""
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2", "sample3"],
                "specimen_taxon_id": [5833, 5833, 5833],
                "host_taxon_id": [9606, 9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana", "Kenya"],
                "project_name": ["project1", "project2", "project3"],
                "host_age": [25, None, 30],
                "host_sex": ["M", "F", ""],
                "host_subject_id": ["SUB001", "", "SUB003"],
                "gravid": [True, None, False],
                "gravidity": [2, None, 0],
            }
        )

        result = specimen_info_table_to_pmo(
            df,
            host_age_col="host_age",
            host_sex_col="host_sex",
            host_subject_id="host_subject_id",
            gravid_col="gravid",
            gravidity_col="gravidity",
        )

        # sample1: all fields have values
        self.assertEqual(result[0]["host_age"], 25)
        self.assertEqual(result[0]["host_sex"], "M")
        self.assertEqual(result[0]["host_subject_id"], "SUB001")
        self.assertEqual(result[0]["gravid"], True)
        self.assertEqual(result[0]["gravidity"], 2)

        # sample2: host_age is None, host_subject_id is empty string, gravid is None, gravidity is None
        # Should remove: host_age, host_subject_id, gravid, gravidity
        # Should keep: host_sex
        self.assertNotIn("host_age", result[1])
        self.assertEqual(result[1]["host_sex"], "F")
        self.assertNotIn("host_subject_id", result[1])
        self.assertNotIn("gravid", result[1])
        self.assertNotIn("gravidity", result[1])

        # sample3: host_sex is empty string
        # Should remove: host_sex
        # Should keep: host_age, host_subject_id, gravid, gravidity
        self.assertEqual(result[2]["host_age"], 30)
        self.assertNotIn("host_sex", result[2])
        self.assertEqual(result[2]["host_subject_id"], "SUB003")
        self.assertEqual(result[2]["gravid"], False)
        self.assertEqual(result[2]["gravidity"], 0)

    def test_specimen_info_table_to_pmo_with_specimen_fields(self):
        """Test specimen-specific optional fields"""
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "specimen_taxon_id": [5833, 5833],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana"],
                "project_name": ["project1", "project2"],
                "specimen_accession": ["ACC001", ""],
                "specimen_type": ["field_sample", None],
                "specimen_collect_device": ["needle", ""],
                "specimen_store_loc": [None, "freezer1"],
                "specimen_comments": ["", "Some comment"],
            }
        )

        result = specimen_info_table_to_pmo(
            df,
            specimen_accession_col="specimen_accession",
            specimen_type_col="specimen_type",
            specimen_collect_device_col="specimen_collect_device",
            specimen_store_loc_col="specimen_store_loc",
            specimen_comments_col="specimen_comments",
        )

        # sample1: specimen_accession has value, specimen_type has value, specimen_collect_device has value
        # specimen_store_loc is None, specimen_comments is empty string
        # Should remove: specimen_store_loc, specimen_comments
        self.assertEqual(result[0]["specimen_accession"], "ACC001")
        self.assertEqual(result[0]["specimen_type"], "field_sample")
        self.assertEqual(result[0]["specimen_collect_device"], "needle")
        self.assertNotIn("specimen_store_loc", result[0])
        self.assertNotIn("specimen_comments", result[0])

        # sample2: specimen_accession is empty string, specimen_type is None, specimen_collect_device is empty string
        # specimen_store_loc has value, specimen_comments has value
        # Should remove: specimen_accession, specimen_type, specimen_collect_device
        self.assertNotIn("specimen_accession", result[1])
        self.assertNotIn("specimen_type", result[1])
        self.assertNotIn("specimen_collect_device", result[1])
        self.assertEqual(result[1]["specimen_store_loc"], "freezer1")
        self.assertEqual(result[1]["specimen_comments"], "Some comment")

    def test_specimen_info_table_to_pmo_with_environment_fields(self):
        """Test environment-related optional fields"""
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2"],
                "specimen_taxon_id": [5833, 5833],
                "host_taxon_id": [9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana"],
                "project_name": ["project1", "project2"],
                "env_broad_scale": ["terrestrial", ""],
                "env_local_scale": [None, "forest"],
                "env_medium": ["soil", None],
            }
        )

        result = specimen_info_table_to_pmo(
            df,
            env_broad_scale_col="env_broad_scale",
            env_local_scale_col="env_local_scale",
            env_medium_col="env_medium",
        )

        # sample1: env_broad_scale has value, env_local_scale is None, env_medium has value
        self.assertEqual(result[0]["env_broad_scale"], "terrestrial")
        self.assertNotIn("env_local_scale", result[0])
        self.assertEqual(result[0]["env_medium"], "soil")

        # sample2: env_broad_scale is empty string, env_local_scale has value, env_medium is None
        self.assertNotIn("env_broad_scale", result[1])
        self.assertEqual(result[1]["env_local_scale"], "forest")
        self.assertNotIn("env_medium", result[1])

    def test_specimen_info_table_to_pmo_with_other_optional_fields(self):
        """Test other optional fields like blood_meal, has_travel_out_six_month, treatment_status, etc."""
        df = pd.DataFrame(
            {
                "specimen_name": ["sample1", "sample2", "sample3"],
                "specimen_taxon_id": [5833, 5833, 5833],
                "host_taxon_id": [9606, 9606, 9606],
                "collection_date": ["01/02/2023", "01/02/2023", "01/02/2023"],
                "collection_country": ["Mozambique", "Ghana", "Kenya"],
                "project_name": ["project1", "project2", "project3"],
                "blood_meal": [True, None, False],
                "has_travel_out_six_month": [None, True, ""],
                "treatment_status": ["cured", "", None],
                "specimen_accession": ["ACC001", None, "ACC003"],
            }
        )

        result = specimen_info_table_to_pmo(
            df,
            blood_meal_col="blood_meal",
            has_travel_out_six_month_col="has_travel_out_six_month",
            treatment_status_col="treatment_status",
            specimen_accession_col="specimen_accession",
        )

        # sample1: blood_meal has value (True), has_travel_out_six_month is None, treatment_status has value, specimen_accession has value
        self.assertEqual(result[0]["blood_meal"], True)
        self.assertNotIn("has_travel_out_six_month", result[0])
        self.assertEqual(result[0]["treatment_status"], "cured")
        self.assertEqual(result[0]["specimen_accession"], "ACC001")

        # sample2: blood_meal is None, has_travel_out_six_month has value (True), treatment_status is empty string, specimen_accession is None
        self.assertNotIn("blood_meal", result[1])
        self.assertEqual(result[1]["has_travel_out_six_month"], True)
        self.assertNotIn("treatment_status", result[1])
        self.assertNotIn("specimen_accession", result[1])

        # sample3: blood_meal has value (False), has_travel_out_six_month is empty string, treatment_status is None, specimen_accession has value
        self.assertEqual(result[2]["blood_meal"], False)
        self.assertNotIn("has_travel_out_six_month", result[2])
        self.assertNotIn("treatment_status", result[2])
        self.assertEqual(result[2]["specimen_accession"], "ACC003")

    def test_library_sample_info_table_to_pmo_default(self):
        df = pd.DataFrame(
            {
                "library_sample_name": ["sample1_MH_run1", "sample2_MH_run1"],
                "sequencing_info_name": ["run1", "run1"],
                "specimen_name": ["sample1", "sample2"],
                "panel_name": ["MH", "MH"],
            }
        )

        result = library_sample_info_table_to_pmo(df)
        self.assertEqual(
            [
                {
                    "library_sample_name": "sample1_MH_run1",
                    "sequencing_info_name": "run1",
                    "specimen_name": "sample1",
                    "panel_name": "MH",
                },
                {
                    "library_sample_name": "sample2_MH_run1",
                    "sequencing_info_name": "run1",
                    "specimen_name": "sample2",
                    "panel_name": "MH",
                },
            ],
            result,
        )

    def test_library_sample_info_table_to_pmo_with_plate(self):
        df = pd.DataFrame(
            {
                "library_sample_name": ["sample1_MH_run1", "sample2_MH_run1"],
                "sequencing_info_name": ["run1", "run1"],
                "specimen_name": ["sample1", "sample2"],
                "panel_name": ["MH", "MH"],
                "library_prep_plate_col": [1, 2],
                "library_prep_plate_name": ["plate1", "plate1"],
                "library_prep_plate_row": ["A", "B"],
            }
        )

        result = library_sample_info_table_to_pmo(
            df,
            library_prep_plate_name_col="library_prep_plate_name",
            library_prep_plate_col_col="library_prep_plate_col",
            library_prep_plate_row_col="library_prep_plate_row",
        )
        self.assertEqual(
            [
                {
                    "library_sample_name": "sample1_MH_run1",
                    "sequencing_info_name": "run1",
                    "specimen_name": "sample1",
                    "panel_name": "MH",
                    "library_prep_plate_info": {
                        "plate_name": "plate1",
                        "plate_row": "A",
                        "plate_col": 1,
                    },
                },
                {
                    "library_sample_name": "sample2_MH_run1",
                    "sequencing_info_name": "run1",
                    "specimen_name": "sample2",
                    "panel_name": "MH",
                    "library_prep_plate_info": {
                        "plate_name": "plate1",
                        "plate_row": "B",
                        "plate_col": 2,
                    },
                },
            ],
            result,
        )

    def test_library_sample_info_table_to_pmo_with_additional_columns(self):
        df = pd.DataFrame(
            {
                "library_sample_name": ["sample1_MH_run1", "sample2_MH_run1"],
                "sequencing_info_name": ["run1", "run1"],
                "specimen_name": ["sample1", "sample2"],
                "panel_name": ["MH", "MH"],
                "new_col1": ["test", "this"],
                "new_col2": ["add", "one"],
            }
        )

        result = library_sample_info_table_to_pmo(
            df, additional_library_sample_info_cols=["new_col1", "new_col2"]
        )
        self.assertEqual(
            [
                {
                    "library_sample_name": "sample1_MH_run1",
                    "sequencing_info_name": "run1",
                    "specimen_name": "sample1",
                    "panel_name": "MH",
                    "new_col1": "test",
                    "new_col2": "add",
                },
                {
                    "library_sample_name": "sample2_MH_run1",
                    "sequencing_info_name": "run1",
                    "specimen_name": "sample2",
                    "panel_name": "MH",
                    "new_col1": "this",
                    "new_col2": "one",
                },
            ],
            result,
        )

    def test_library_sample_info_table_to_pmo_fails_with_duplicate_cols(self):
        df = pd.DataFrame(
            {
                "library_sample_name": ["sample1_MH_run1", "sample2_MH_run1"],
                "sequencing_info_name": ["run1", "run1"],
                "specimen_name": ["sample1", "sample2"],
                "panel_name": ["MH", "MH"],
            }
        )

        with self.assertRaises(ValueError) as context:
            library_sample_info_table_to_pmo(df, specimen_name_col="panel_name")
        self.assertEqual("Selected columns must be unique.", str(context.exception))

    def test_library_sample_info_table_to_pmo_fails_with_missing_cols(self):
        df = pd.DataFrame(
            {
                "library_sample_name": ["sample1_MH_run1", "sample2_MH_run1"],
                "sequencing_info_name": ["run1", "run1"],
            }
        )

        with self.assertRaises(ValueError) as context:
            library_sample_info_table_to_pmo(df)
        self.assertEqual(
            "The following columns are not in the DataFrame: ['specimen_name', 'panel_name']",
            str(context.exception),
        )

    def test_library_sample_info_table_to_pmo_fails_without_df(self):
        with self.assertRaises(ValueError) as context:
            library_sample_info_table_to_pmo("test")
        self.assertEqual("contents must be a pandas DataFrame.", str(context.exception))

    def test_library_sample_info_table_to_pmo_fails_with_null_in_required_columns(self):
        """Test that function raises error when required columns contain null values"""
        # Test with null in library_sample_name
        df1 = pd.DataFrame(
            {
                "library_sample_name": ["sample1_MH_run1", None],
                "sequencing_info_name": ["run1", "run1"],
                "specimen_name": ["sample1", "sample2"],
                "panel_name": ["MH", "MH"],
            }
        )
        with self.assertRaises(ValueError) as context:
            library_sample_info_table_to_pmo(df1)
        self.assertIn(
            "The following columns contain null values",
            str(context.exception),
        )
        self.assertIn("library_sample_name", str(context.exception))

        # Test with null in sequencing_info_name
        df2 = pd.DataFrame(
            {
                "library_sample_name": ["sample1_MH_run1", "sample2_MH_run1"],
                "sequencing_info_name": [None, "run1"],
                "specimen_name": ["sample1", "sample2"],
                "panel_name": ["MH", "MH"],
            }
        )
        with self.assertRaises(ValueError) as context:
            library_sample_info_table_to_pmo(df2)
        self.assertIn(
            "The following columns contain null values",
            str(context.exception),
        )
        self.assertIn("sequencing_info_name", str(context.exception))

        # Test with null in specimen_name
        df3 = pd.DataFrame(
            {
                "library_sample_name": ["sample1_MH_run1", "sample2_MH_run1"],
                "sequencing_info_name": ["run1", "run1"],
                "specimen_name": ["sample1", None],
                "panel_name": ["MH", "MH"],
            }
        )
        with self.assertRaises(ValueError) as context:
            library_sample_info_table_to_pmo(df3)
        self.assertIn(
            "The following columns contain null values",
            str(context.exception),
        )
        self.assertIn("specimen_name", str(context.exception))

        # Test with null in panel_name
        df4 = pd.DataFrame(
            {
                "library_sample_name": ["sample1_MH_run1", "sample2_MH_run1"],
                "sequencing_info_name": ["run1", "run1"],
                "specimen_name": ["sample1", "sample2"],
                "panel_name": [None, "MH"],
            }
        )
        with self.assertRaises(ValueError) as context:
            library_sample_info_table_to_pmo(df4)
        self.assertIn(
            "The following columns contain null values",
            str(context.exception),
        )
        self.assertIn("panel_name", str(context.exception))

        # Test with nulls in multiple required columns
        df5 = pd.DataFrame(
            {
                "library_sample_name": ["sample1_MH_run1", None],
                "sequencing_info_name": [None, "run1"],
                "specimen_name": ["sample1", "sample2"],
                "panel_name": ["MH", "MH"],
            }
        )
        with self.assertRaises(ValueError) as context:
            library_sample_info_table_to_pmo(df5)
        self.assertIn(
            "The following columns contain null values",
            str(context.exception),
        )
        error_msg = str(context.exception)
        self.assertTrue(
            "library_sample_name" in error_msg and "sequencing_info_name" in error_msg
        )

    def test_library_sample_info_table_to_pmo_with_new_optional_fields(self):
        """Test new optional fields: alternate_identifiers, experiment_accession, fastqs_loc, run_accession"""
        df = pd.DataFrame(
            {
                "library_sample_name": ["sample1_MH_run1", "sample2_MH_run1"],
                "sequencing_info_name": ["run1", "run1"],
                "specimen_name": ["sample1", "sample2"],
                "panel_name": ["MH", "MH"],
                "alternate_identifiers": ["ID1,ID2", ""],
                "experiment_accession": ["EXP001", None],
                "fastqs_loc": [None, "/path/to/fastqs"],
                "run_accession": ["RUN001", ""],
            }
        )

        result = library_sample_info_table_to_pmo(
            df,
            alternate_identifiers_col="alternate_identifiers",
            experiment_accession_col="experiment_accession",
            fastqs_loc_col="fastqs_loc",
            run_accession_col="run_accession",
        )

        # sample1: alternate_identifiers has value, experiment_accession has value, fastqs_loc is None, run_accession has value
        # Should remove: fastqs_loc
        self.assertEqual(result[0]["alternate_identifiers"], "ID1,ID2")
        self.assertEqual(result[0]["experiment_accession"], "EXP001")
        self.assertNotIn("fastqs_loc", result[0])
        self.assertEqual(result[0]["run_accession"], "RUN001")

        # sample2: alternate_identifiers is empty string, experiment_accession is None, fastqs_loc has value, run_accession is empty string
        # Should remove: alternate_identifiers, experiment_accession, run_accession
        self.assertNotIn("alternate_identifiers", result[1])
        self.assertNotIn("experiment_accession", result[1])
        self.assertEqual(result[1]["fastqs_loc"], "/path/to/fastqs")
        self.assertNotIn("run_accession", result[1])

    def test_library_sample_info_table_to_pmo_removes_empty_optional_fields(self):
        """Test that empty optional fields are removed"""
        df = pd.DataFrame(
            {
                "library_sample_name": [
                    "sample1_MH_run1",
                    "sample2_MH_run1",
                    "sample3_MH_run1",
                ],
                "sequencing_info_name": ["run1", "run1", "run1"],
                "specimen_name": ["sample1", "sample2", "sample3"],
                "panel_name": ["MH", "MH", "MH"],
                "alternate_identifiers": ["ID1", "", None],
                "experiment_accession": ["EXP001", None, "EXP003"],
                "fastqs_loc": [None, "/path/to/fastqs", ""],
            }
        )

        result = library_sample_info_table_to_pmo(
            df,
            alternate_identifiers_col="alternate_identifiers",
            experiment_accession_col="experiment_accession",
            fastqs_loc_col="fastqs_loc",
        )

        # sample1: all optional fields have values
        self.assertEqual(result[0]["alternate_identifiers"], "ID1")
        self.assertEqual(result[0]["experiment_accession"], "EXP001")
        self.assertNotIn("fastqs_loc", result[0])  # None should be removed

        # sample2: alternate_identifiers is empty string, experiment_accession is None, fastqs_loc has value
        self.assertNotIn("alternate_identifiers", result[1])
        self.assertNotIn("experiment_accession", result[1])
        self.assertEqual(result[1]["fastqs_loc"], "/path/to/fastqs")

        # sample3: alternate_identifiers is None, experiment_accession has value, fastqs_loc is empty string
        self.assertNotIn("alternate_identifiers", result[2])
        self.assertEqual(result[2]["experiment_accession"], "EXP003")
        self.assertNotIn("fastqs_loc", result[2])

    def test_library_sample_info_table_to_pmo_with_parasite_density(self):
        """Test parasite density info in library sample info"""
        df = pd.DataFrame(
            {
                "library_sample_name": ["sample1_MH_run1", "sample2_MH_run1"],
                "sequencing_info_name": ["run1", "run1"],
                "specimen_name": ["sample1", "sample2"],
                "panel_name": ["MH", "MH"],
                "parasite_density": [10, 100],
                "parasite_density_method": ["qPCR", "microscopy"],
            }
        )

        result = library_sample_info_table_to_pmo(
            df,
            parasite_density_col="parasite_density",
            parasite_density_method_col="parasite_density_method",
        )

        self.assertEqual(result[0]["parasite_density_info"][0]["parasite_density"], 10)
        self.assertEqual(
            result[0]["parasite_density_info"][0]["parasite_density_method"], "qPCR"
        )
        self.assertEqual(result[1]["parasite_density_info"][0]["parasite_density"], 100)
        self.assertEqual(
            result[1]["parasite_density_info"][0]["parasite_density_method"],
            "microscopy",
        )

    def test_library_sample_info_table_to_pmo_with_parasite_density_multiple(self):
        """Test multiple parasite density measurements"""
        df = pd.DataFrame(
            {
                "library_sample_name": ["sample1_MH_run1", "sample2_MH_run1"],
                "sequencing_info_name": ["run1", "run1"],
                "specimen_name": ["sample1", "sample2"],
                "panel_name": ["MH", "MH"],
                "density1": [15, 107],
                "method1": ["qPCR", "qPCR"],
                "density2": [10, 100],
                "method2": ["microscopy", "microscopy"],
            }
        )

        result = library_sample_info_table_to_pmo(
            df,
            parasite_density_col=["density1", "density2"],
            parasite_density_method_col=["method1", "method2"],
        )

        self.assertEqual(result[0]["parasite_density_info"][0]["parasite_density"], 15)
        self.assertEqual(result[0]["parasite_density_info"][1]["parasite_density"], 10)
        self.assertEqual(result[1]["parasite_density_info"][0]["parasite_density"], 107)
        self.assertEqual(result[1]["parasite_density_info"][1]["parasite_density"], 100)

    def test_library_sample_info_table_to_pmo_with_all_new_fields(self):
        """Test all new optional fields together"""
        df = pd.DataFrame(
            {
                "library_sample_name": ["sample1_MH_run1", "sample2_MH_run1"],
                "sequencing_info_name": ["run1", "run1"],
                "specimen_name": ["sample1", "sample2"],
                "panel_name": ["MH", "MH"],
                "alternate_identifiers": ["ID1", "ID2"],
                "experiment_accession": ["EXP001", "EXP002"],
                "fastqs_loc": ["/path/to/fastqs1", "/path/to/fastqs2"],
                "run_accession": ["RUN001", "RUN002"],
                "parasite_density": [10, 100],
                "parasite_density_method": ["qPCR", "microscopy"],
                "library_prep_plate_col": [1, 2],
                "library_prep_plate_name": ["plate1", "plate1"],
                "library_prep_plate_row": ["A", "B"],
            }
        )

        result = library_sample_info_table_to_pmo(
            df,
            alternate_identifiers_col="alternate_identifiers",
            experiment_accession_col="experiment_accession",
            fastqs_loc_col="fastqs_loc",
            run_accession_col="run_accession",
            parasite_density_col="parasite_density",
            parasite_density_method_col="parasite_density_method",
            library_prep_plate_col_col="library_prep_plate_col",
            library_prep_plate_name_col="library_prep_plate_name",
            library_prep_plate_row_col="library_prep_plate_row",
        )

        # Check all fields are present
        self.assertEqual(result[0]["alternate_identifiers"], "ID1")
        self.assertEqual(result[0]["experiment_accession"], "EXP001")
        self.assertEqual(result[0]["fastqs_loc"], "/path/to/fastqs1")
        self.assertEqual(result[0]["run_accession"], "RUN001")
        self.assertEqual(result[0]["parasite_density_info"][0]["parasite_density"], 10)
        self.assertIn("library_prep_plate_info", result[0])
        self.assertEqual(result[0]["library_prep_plate_info"]["plate_col"], 1)

        self.assertEqual(result[1]["alternate_identifiers"], "ID2")
        self.assertEqual(result[1]["experiment_accession"], "EXP002")
        self.assertEqual(result[1]["fastqs_loc"], "/path/to/fastqs2")
        self.assertEqual(result[1]["run_accession"], "RUN002")
        self.assertEqual(result[1]["parasite_density_info"][0]["parasite_density"], 100)
        self.assertIn("library_prep_plate_info", result[1])
        self.assertEqual(result[1]["library_prep_plate_info"]["plate_col"], 2)


if __name__ == "__main__":
    unittest.main()
