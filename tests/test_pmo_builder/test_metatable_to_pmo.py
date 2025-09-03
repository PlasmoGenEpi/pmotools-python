import unittest
import pandas as pd

from pmotools.pmo_builder.metatable_to_pmo import *


class TestMetatableToPMO(unittest.TestCase):

    def setUp(self):
        self.small_json_example = [
            {'specimen_name': 'sample1'}, {'specimen_name': 'sample2'}]

    def test_add_plate_info_position_parsing_uppercase(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'plate_position': ['A01', 'H11'],
            'plate_name': ['Plate1', 'Plate1']
        })
        result = add_plate_info(None, 'plate_name', None,
                                'plate_position', self.small_json_example, df, 'specimen_name')
        self.assertEqual(result[0]['plate_info']['plate_row'], 'A')
        self.assertEqual(result[0]['plate_info']['plate_col'], 1)
        self.assertEqual(result[1]['plate_info']['plate_row'], 'H')
        self.assertEqual(result[1]['plate_info']['plate_col'], 11)

    def test_add_plate_info_position_parsing_lowercase(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'plate_position': ['a01', 'h11'],
            'plate_name': ['Plate1', 'Plate1']
        })
        result = add_plate_info(None, 'plate_name', None,
                                'plate_position', self.small_json_example, df, 'specimen_name')
        self.assertEqual(result[0]['plate_info']['plate_row'], 'A')
        self.assertEqual(result[0]['plate_info']['plate_col'], 1)
        self.assertEqual(result[1]['plate_info']['plate_row'], 'H')
        self.assertEqual(result[1]['plate_info']['plate_col'], 11)

    def test_add_plate_info_position_parsing_one_digit(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'plate_position': ['a1', 'h10'],
            'plate_name': ['Plate1', 'Plate1']
        })
        result = add_plate_info(None, 'plate_name', None,
                                'plate_position', self.small_json_example, df, 'specimen_name')
        self.assertEqual(result[0]['plate_info']['plate_row'], 'A')
        self.assertEqual(result[0]['plate_info']['plate_col'], 1)
        self.assertEqual(result[1]['plate_info']['plate_row'], 'H')
        self.assertEqual(result[1]['plate_info']['plate_col'], 10)

    def test_add_plate_info_position_fails_with_out_of_bounds_row(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'plate_position': ['A01', 'K10'],
            'plate_name': ['Plate1', 'Plate1']
        })
        with self.assertRaises(ValueError) as context:
            add_plate_info(None, 'plate_name', None,
                           'plate_position', self.small_json_example, df, 'specimen_name')
        self.assertEqual(
            "Values in 'plate_position' must start with a single letter A-H/a-h followed by number 1-12.", str(context.exception))

    def test_add_plate_info_position_fails_with_out_of_bounds_col(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'plate_position': ['A01', 'H19'],
            'plate_name': ['Plate1', 'Plate1']
        })
        with self.assertRaises(ValueError) as context:
            add_plate_info(None, 'plate_name', None,
                           'plate_position', self.small_json_example, df, 'specimen_name')

        self.assertEqual(
            "Values in 'plate_position' must start with a single letter A-H/a-h followed by number 1-12.", str(context.exception))

    def test_add_plate_info_row_col_parsing(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'plate_row': ['A', 'h'],
            'plate_col': [1, 10],
            'plate_name': ['Plate1', 'Plate1']
        })
        result = add_plate_info('plate_col', 'plate_name', 'plate_row',
                                None, self.small_json_example, df, 'specimen_name')
        self.assertEqual(result[0]['plate_info']['plate_row'], 'A')
        self.assertEqual(result[0]['plate_info']['plate_col'], 1)
        self.assertEqual(result[1]['plate_info']['plate_row'], 'H')
        self.assertEqual(result[1]['plate_info']['plate_col'], 10)

    def test_add_plate_info_fails_with_position_and_row_col(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'plate_position': ['A01', 'H19'],
            'plate_row': ['A', 'h'],
            'plate_col': [1, 10],
            'plate_name': ['Plate1', 'Plate1']
        })
        with self.assertRaises(ValueError) as context:
            add_plate_info('plate_col', 'plate_name', 'plate_row',
                           'plate_position', self.small_json_example, df, 'specimen_name')
        self.assertEqual(
            "Plate position can be specified using either row and col, or position, but not both.", str(context.exception))

    def test_add_plate_info_fails_with_row_and_without_col(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'plate_position': ['A01', 'H19'],
            'plate_row': ['A', 'h'],
            'plate_col': [1, 10],
            'plate_name': ['Plate1', 'Plate1']
        })
        with self.assertRaises(ValueError) as context:
            add_plate_info(None, 'plate_name', 'plate_row',
                           None, self.small_json_example, df, 'specimen_name')
        self.assertEqual(
            "If either plate row or column is set, then both must be.", str(context.exception))

    def test_add_plate_info_adds_nothing(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'plate_position': ['A01', 'H19'],
            'plate_row': ['A', 'h'],
            'plate_col': [1, 10],
            'plate_name': ['Plate1', 'Plate1']
        })
        result = add_plate_info(None, None, None,
                                None, self.small_json_example, df, 'specimen_name')
        self.assertEqual(result, [{'specimen_name': 'sample1'},  {
                         'specimen_name': 'sample2'}])

    def test_add_parasite_density_info_single_value(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'parasite_density': [10, 100],
            'parasite_density_method': ['qPCR', 'microscopy']
        })
        result = add_parasite_density_info(
            'parasite_density', 'parasite_density_method', self.small_json_example, df, 'specimen_name',
            "parasite_density_info")
        self.assertEqual(result[0]['parasite_density_info']
                         [0]['parasite_density'], 10)
        self.assertEqual(
            result[1]['parasite_density_info'][0]['parasite_density'], 100)
        self.assertEqual(
            result[0]['parasite_density_info'][0]['parasite_density_method'], 'qPCR')
        self.assertEqual(result[1]['parasite_density_info'][0]
                         ['parasite_density_method'], 'microscopy')

    def test_add_parasite_density_info_from_list(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'density1': [15, 107],
            'method1': ['qPCR', 'qPCR'],
            'density2': [10, 100],
            'method2': ['microscopy', 'microscopy']
        })
        result = add_parasite_density_info(
            ['density1', 'density2'], ['method1',
                                       'method2'], self.small_json_example, df, 'specimen_name',
            "parasite_density_info")
        self.assertEqual(
            result[0]['parasite_density_info'][0]['parasite_density'], 15)
        self.assertEqual(result[1]['parasite_density_info']
                         [0]['parasite_density'], 107)
        self.assertEqual(
            result[0]['parasite_density_info'][1]['parasite_density'], 10)
        self.assertEqual(
            result[1]['parasite_density_info'][1]['parasite_density'], 100)
        self.assertEqual(
            result[0]['parasite_density_info'][0]['parasite_density_method'], 'qPCR')
        self.assertEqual(result[1]['parasite_density_info'][0]
                         ['parasite_density_method'], 'qPCR')
        self.assertEqual(
            result[0]['parasite_density_info'][1]['parasite_density_method'], 'microscopy')
        self.assertEqual(result[1]['parasite_density_info'][1]
                         ['parasite_density_method'], 'microscopy')

    def test_add_parasite_density_adds_nothing(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
        })
        result = add_parasite_density_info(
            None, None, self.small_json_example, df, 'specimen_name',
            "parasite_density_info")
        self.assertEqual(result, [{'specimen_name': 'sample1'}, {
                         'specimen_name': 'sample2'}])

    def test_add_parasite_density_fails_with_unequal_lists(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'density1': [15, 107],
            'method1': ['qPCR', 'qPCR'],
            'density2': [10, 100],
            'method2': ['microscopy', 'microscopy']
        })
        with self.assertRaises(ValueError) as context:
            add_parasite_density_info(
                ['density1'], ['method1',
                               'method2'], self.small_json_example, df, 'specimen_name',
                "qpcr_parasite_density_info")

        self.assertEqual(
            "If both parasite_density_col and parasite_density_method_col are lists, they must be the same length.", str(context.exception))

    def test_add_parasite_density_fails_with_no_density(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'method1': ['qPCR', 'qPCR'],
        })

        with self.assertRaises(ValueError) as context:
            add_parasite_density_info(
                None, 'method1', self.small_json_example, df, 'specimen_name',
                "parasite_density_info")

        self.assertEqual(
            "parasite_density_method_col is set but parasite_density_col is None. Cannot proceed.", str(context.exception))

    def test_add_parasite_density_fails_with_type_mismatch(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'density1': [15, 107],
            'method1': ['qPCR', 'qPCR'],
            'density2': [10, 100],
            'method2': ['microscopy', 'microscopy']
        })
        with self.assertRaises(TypeError) as context:
            add_parasite_density_info(
                'density1', [
                    'method1', 'method2'], self.small_json_example, df, 'specimen_name',
                "parasite_density_info")

        self.assertEqual(
            "If parasite_density_col is a string, parasite_density_method_col must be a string or None.", str(context.exception))

    def test_specimen_info_table_to_pmo_default(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'specimen_taxon_id': [5833, 5833],
            'host_taxon_id': [9606, 9606],
            'collection_date': ['01/02/2023', '01/02/2023'],
            'collection_country': ['Mozambique', 'Ghana'],
            'project_name': ['project1', 'project2']
        })

        result = specimen_info_table_to_pmo(df)
        self.assertEqual([{'specimen_name': 'sample1',
                           'specimen_taxon_id': 5833,
                           'host_taxon_id': 9606,
                           'collection_date': '01/02/2023',
                           'collection_country': 'Mozambique',
                           'project_name': 'project1'},
                          {'specimen_name': 'sample2',
                           'specimen_taxon_id': 5833,
                           'host_taxon_id': 9606,
                           'collection_date': '01/02/2023',
                           'collection_country': 'Ghana',
                           'project_name': 'project2'}], result)

    def test_specimen_info_table_to_pmo_with_plate_info(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'specimen_taxon_id': [5833, 5833],
            'host_taxon_id': [9606, 9606],
            'collection_date': ['01/02/2023', '01/02/2023'],
            'collection_country': ['Mozambique', 'Ghana'],
            'project_name': ['project1', 'project2'],
            'storage_plate_col': [1, 2],
            'storage_plate_name': ['plate1', 'plate1'],
            'storage_plate_row': ['A', 'B']
        })

        result = specimen_info_table_to_pmo(
            df, storage_plate_col_col='storage_plate_col', storage_plate_name_col='storage_plate_name', storage_plate_row_col='storage_plate_row')
        self.assertEqual([{'specimen_name': 'sample1',
                           'specimen_taxon_id': 5833,
                           'host_taxon_id': 9606,
                           'collection_date': '01/02/2023',
                           'collection_country': 'Mozambique',
                           'project_name': 'project1',
                           'storage_plate_info': {
                               'plate_name': 'plate1',
                               'plate_row': 'A',
                               'plate_col': 1
                           }},
                          {'specimen_name': 'sample2',
                           'specimen_taxon_id': 5833,
                           'host_taxon_id': 9606,
                           'collection_date': '01/02/2023',
                           'collection_country': 'Ghana',
                           'project_name': 'project2',
                           'storage_plate_info': {
                               'plate_name': 'plate1',
                               'plate_row': 'B',
                               'plate_col': 2
                           }}], result)

    def test_specimen_info_table_to_pmo_with_parasitemia(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'specimen_taxon_id': [5833, 5833],
            'host_taxon_id': [9606, 9606],
            'collection_date': ['01/02/2023', '01/02/2023'],
            'collection_country': ['Mozambique', 'Ghana'],
            'project_name': ['project1', 'project2'],
            'parasite_density': [10, 100],
            'parasite_density_method': ['qPCR', 'microscopy']
        })

        result = specimen_info_table_to_pmo(
            df, parasite_density_col='parasite_density', parasite_density_method_col='parasite_density_method')
        self.assertEqual([{'specimen_name': 'sample1',
                           'specimen_taxon_id': 5833,
                           'host_taxon_id': 9606,
                           'collection_date': '01/02/2023',
                           'collection_country': 'Mozambique',
                           'project_name': 'project1',
                           'parasite_density_info': [{
                               'parasite_density': 10,
                               'parasite_density_method': 'qPCR'
                           }]},
                          {'specimen_name': 'sample2',
                           'specimen_taxon_id': 5833,
                           'host_taxon_id': 9606,
                           'collection_date': '01/02/2023',
                           'collection_country': 'Ghana',
                           'project_name': 'project2',
                           'parasite_density_info': [{
                               'parasite_density': 100,
                               'parasite_density_method': 'microscopy'
                           }]}], result)

    def test_specimen_info_table_to_pmo_with_additional_columns(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'specimen_taxon_id': [5833, 5833],
            'host_taxon_id': [9606, 9606],
            'collection_date': ['01/02/2023', '01/02/2023'],
            'collection_country': ['Mozambique', 'Ghana'],
            'project_name': ['project1', 'project2'],
            'special_field_1': ['something', 'something_else'],
            'special_field_2': ['this', 'that']
        })

        result = specimen_info_table_to_pmo(
            df, additional_specimen_cols=['special_field_1', 'special_field_2'])
        self.assertEqual([{'specimen_name': 'sample1',
                           'specimen_taxon_id': 5833,
                           'host_taxon_id': 9606,
                           'collection_date': '01/02/2023',
                           'collection_country': 'Mozambique',
                           'project_name': 'project1',
                           'special_field_1': 'something',
                           'special_field_2': 'this'},
                          {'specimen_name': 'sample2',
                           'specimen_taxon_id': 5833,
                           'host_taxon_id': 9606,
                           'collection_date': '01/02/2023',
                           'collection_country': 'Ghana',
                           'project_name': 'project2',
                           'special_field_1': 'something_else',
                           'special_field_2': 'that'
                           }], result)

    def test_specimen_info_table_to_pmo_fails_with_col_duplicate(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'specimen_taxon_id': [5833, 5833],
            'host_taxon_id': [9606, 9606],
            'collection_date': ['01/02/2023', '01/02/2023'],
            'collection_country': ['Mozambique', 'Ghana'],
            'project_name': ['project1', 'project2']
        })
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_pmo(
                df, drug_usage_col='specimen_name')
        self.assertEqual(
            "Selected columns must be unique.", str(context.exception))

    def test_specimen_info_table_to_pmo_fails_with_missing_col(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'specimen_taxon_id': [5833, 5833],
            'host_taxon_id': [9606, 9606],
            'collection_date': ['01/02/2023', '01/02/2023'],
            'collection_country': ['Mozambique', 'Ghana'],
            'project_name': ['project1', 'project2']
        })
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_pmo(
                df, drug_usage_col='drug_usage')
        self.assertEqual(
            "The following columns are not in the DataFrame: ['drug_usage']", str(context.exception))

    def test_specimen_info_table_to_pmo_fails_without_df(self):
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_pmo('test')
        self.assertEqual(
            "contents must be a pandas DataFrame.", str(context.exception))

    def test_library_sample_info_table_to_pmo_default(self):
        df = pd.DataFrame({
            'library_sample_name': ['sample1_MH_run1', 'sample2_MH_run1'],
            'sequencing_info_name': ['run1', 'run1'],
            'specimen_name': ['sample1', 'sample2'],
            'panel_name': ['MH', 'MH'],
        })

        result = library_sample_info_table_to_pmo(df)
        self.assertEqual([{'library_sample_name': 'sample1_MH_run1',
                           'sequencing_info_name': 'run1',
                           'specimen_name': 'sample1',
                           'panel_name': 'MH'},
                          {'library_sample_name': 'sample2_MH_run1',
                           'sequencing_info_name': 'run1',
                           'specimen_name': 'sample2',
                           'panel_name': 'MH'}], result)

    def test_library_sample_info_table_to_pmo_with_plate(self):
        df = pd.DataFrame({
            'library_sample_name': ['sample1_MH_run1', 'sample2_MH_run1'],
            'sequencing_info_name': ['run1', 'run1'],
            'specimen_name': ['sample1', 'sample2'],
            'panel_name': ['MH', 'MH'],
            'library_prep_plate_col': [1, 2],
            'library_prep_plate_name': ['plate1', 'plate1'],
            'library_prep_plate_row': ['A', 'B']
        })

        result = library_sample_info_table_to_pmo(df,
                                                  library_prep_plate_name_col='library_prep_plate_name', library_prep_plate_col_col='library_prep_plate_col', library_prep_plate_row_col='library_prep_plate_row')
        self.assertEqual([{'library_sample_name': 'sample1_MH_run1',
                           'sequencing_info_name': 'run1',
                           'specimen_name': 'sample1',
                           'panel_name': 'MH',
                           'library_prep_plate_info': {
                               'plate_name': 'plate1',
                               'plate_row': 'A',
                               'plate_col': 1
                           }
                           },
                          {
                         'library_sample_name': 'sample2_MH_run1',
                         'sequencing_info_name': 'run1',
                         'specimen_name': 'sample2',
                         'panel_name': 'MH',
                         'library_prep_plate_info': {
                             'plate_name': 'plate1',
                             'plate_row': 'B',
                             'plate_col': 2
                         }}], result)

    def test_library_sample_info_table_to_pmo_with_additional_columns(self):
        df = pd.DataFrame({
            'library_sample_name': ['sample1_MH_run1', 'sample2_MH_run1'],
            'sequencing_info_name': ['run1', 'run1'],
            'specimen_name': ['sample1', 'sample2'],
            'panel_name': ['MH', 'MH'],
            'new_col1': ['test', 'this'],
            'new_col2': ['add', 'one'],
        })

        result = library_sample_info_table_to_pmo(
            df, additional_library_sample_info_cols=['new_col1', 'new_col2'])
        self.assertEqual([{'library_sample_name': 'sample1_MH_run1',
                           'sequencing_info_name': 'run1',
                           'specimen_name': 'sample1',
                           'panel_name': 'MH',
                           'new_col1': 'test',
                           'new_col2': 'add'},
                          {'library_sample_name': 'sample2_MH_run1',
                           'sequencing_info_name': 'run1',
                           'specimen_name': 'sample2',
                           'panel_name': 'MH',
                           'new_col1': 'this',
                           'new_col2': 'one'}], result)

    def test_library_sample_info_table_to_pmo_fails_with_duplicate_cols(self):
        df = pd.DataFrame({
            'library_sample_name': ['sample1_MH_run1', 'sample2_MH_run1'],
            'sequencing_info_name': ['run1', 'run1'],
            'specimen_name': ['sample1', 'sample2'],
            'panel_name': ['MH', 'MH'],
        })

        with self.assertRaises(ValueError) as context:
            library_sample_info_table_to_pmo(
                df, specimen_name_col='panel_name')
        self.assertEqual(
            "Selected columns must be unique.", str(context.exception))

    def test_library_sample_info_table_to_pmo_fails_with_missing_cols(self):
        df = pd.DataFrame({
            'library_sample_name': ['sample1_MH_run1', 'sample2_MH_run1'],
            'sequencing_info_name': ['run1', 'run1'],
        })

        with self.assertRaises(ValueError) as context:
            library_sample_info_table_to_pmo(df)
        self.assertEqual(
            "The following columns are not in the DataFrame: ['specimen_name', 'panel_name']", str(context.exception))

    def test_library_sample_info_table_to_pmo_fails_without_df(self):
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_pmo('test')
        self.assertEqual(
            "contents must be a pandas DataFrame.", str(context.exception))


if __name__ == '__main__':
    unittest.main()
