import unittest
import pandas as pd

from pmotools.pmo_builder.metatable_to_json import *


class TestMetatableToJson(unittest.TestCase):

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
            'density': [10, 100],
            'method': ['qPCR', 'microscopy']
        })
        result = add_parasite_density_info(
            'density', 'method', self.small_json_example, df, 'specimen_name')
        self.assertEqual(result[0]['parasite_density_info'][0]['density'], 10)
        self.assertEqual(
            result[1]['parasite_density_info'][0]['density'], 100)
        self.assertEqual(
            result[0]['parasite_density_info'][0]['method'], 'qPCR')
        self.assertEqual(result[1]['parasite_density_info'][0]
                         ['method'], 'microscopy')

    def test_add_parasite_density_info_from_list(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'density1': [15, 107],
            'method1': ['qPCR', 'qPCR'],
            'density2': [10, 100],
            'method2': ['microscopy', 'microscopy']
        })
        result = add_parasite_density_info(
            ['density1', 'density2'], ['method1', 'method2'], self.small_json_example, df, 'specimen_name')
        self.assertEqual(
            result[0]['parasite_density_info'][0]['density'], 15)
        self.assertEqual(result[1]['parasite_density_info'][0]['density'], 107)
        self.assertEqual(
            result[0]['parasite_density_info'][1]['density'], 10)
        self.assertEqual(
            result[1]['parasite_density_info'][1]['density'], 100)
        self.assertEqual(
            result[0]['parasite_density_info'][0]['method'], 'qPCR')
        self.assertEqual(result[1]['parasite_density_info'][0]
                         ['method'], 'qPCR')
        self.assertEqual(
            result[0]['parasite_density_info'][1]['method'], 'microscopy')
        self.assertEqual(result[1]['parasite_density_info'][1]
                         ['method'], 'microscopy')

    def test_add_parasite_density_adds_nothing(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
        })
        result = add_parasite_density_info(
            None, None, self.small_json_example, df, 'specimen_name')
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
                ['density1'], ['method1', 'method2'], self.small_json_example, df, 'specimen_name')

        self.assertEqual(
            "If both parasite_density_col and parasite_density_method_col are lists, they must be the same length.", str(context.exception))

    def test_add_parasite_density_fails_with_no_density(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'method1': ['qPCR', 'qPCR'],
        })

        with self.assertRaises(ValueError) as context:
            add_parasite_density_info(
                None, 'method1', self.small_json_example, df, 'specimen_name')

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
                'density1', ['method1', 'method2'], self.small_json_example, df, 'specimen_name')

        self.assertEqual(
            "If parasite_density_col is a string, parasite_density_method_col must be a string or None.", str(context.exception))

    def test_specimen_info_table_to_json_default(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'specimen_taxon_id': [5833, 5833],
            'host_taxon_id': [9606, 9606],
            'collection_date': ['01/02/2023', '01/02/2023'],
            'collection_country': ['Mozambique', 'Ghana'],
            'project_name': ['project1', 'project2']
        })

        result = specimen_info_table_to_json(df)
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

    def test_specimen_info_table_to_json_with_plate_info(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'specimen_taxon_id': [5833, 5833],
            'host_taxon_id': [9606, 9606],
            'collection_date': ['01/02/2023', '01/02/2023'],
            'collection_country': ['Mozambique', 'Ghana'],
            'project_name': ['project1', 'project2'],
            'plate_col': [1, 2],
            'plate_name': ['plate1', 'plate1'],
            'plate_row': ['A', 'B']
        })

        result = specimen_info_table_to_json(
            df, plate_col_col='plate_col', plate_name_col='plate_name', plate_row_col='plate_row')
        self.assertEqual([{'specimen_name': 'sample1',
                           'specimen_taxon_id': 5833,
                           'host_taxon_id': 9606,
                           'collection_date': '01/02/2023',
                           'collection_country': 'Mozambique',
                           'project_name': 'project1',
                           'plate_info': {
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
                           'plate_info': {
                               'plate_name': 'plate1',
                               'plate_row': 'B',
                               'plate_col': 2
                           }}], result)

    def test_specimen_info_table_to_json_with_parasitemia(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'specimen_taxon_id': [5833, 5833],
            'host_taxon_id': [9606, 9606],
            'collection_date': ['01/02/2023', '01/02/2023'],
            'collection_country': ['Mozambique', 'Ghana'],
            'project_name': ['project1', 'project2'],
            'density': [10, 100],
            'method': ['qPCR', 'microscopy']
        })

        result = specimen_info_table_to_json(
            df, parasite_density_col='density', parasite_density_method_col='method')
        self.assertEqual([{'specimen_name': 'sample1',
                           'specimen_taxon_id': 5833,
                           'host_taxon_id': 9606,
                           'collection_date': '01/02/2023',
                           'collection_country': 'Mozambique',
                           'project_name': 'project1',
                           'parasite_density_info': [{
                               'density': 10,
                               'method': 'qPCR'
                           }]},
                          {'specimen_name': 'sample2',
                           'specimen_taxon_id': 5833,
                           'host_taxon_id': 9606,
                           'collection_date': '01/02/2023',
                           'collection_country': 'Ghana',
                           'project_name': 'project2',
                           'parasite_density_info': [{
                               'density': 100,
                               'method': 'microscopy'
                           }]}], result)

    def test_specimen_info_table_to_json_with_additional_columns(self):
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

        result = specimen_info_table_to_json(
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

    def test_specimen_info_table_to_json_fails_with_col_duplicate(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'specimen_taxon_id': [5833, 5833],
            'host_taxon_id': [9606, 9606],
            'collection_date': ['01/02/2023', '01/02/2023'],
            'collection_country': ['Mozambique', 'Ghana'],
            'project_name': ['project1', 'project2']
        })
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_json(
                df, drug_usage_col='specimen_name')
        self.assertEqual(
            "Selected columns must be unique.", str(context.exception))

    def test_specimen_info_table_to_json_fails_with_missing_col(self):
        df = pd.DataFrame({
            'specimen_name': ['sample1', 'sample2'],
            'specimen_taxon_id': [5833, 5833],
            'host_taxon_id': [9606, 9606],
            'collection_date': ['01/02/2023', '01/02/2023'],
            'collection_country': ['Mozambique', 'Ghana'],
            'project_name': ['project1', 'project2']
        })
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_json(
                df, drug_usage_col='drug_usage')
        self.assertEqual(
            "The following columns are not in the DataFrame: ['drug_usage']", str(context.exception))

    def test_specimen_info_table_to_json_fails_without_df(self):
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_json('test')
        self.assertEqual(
            "contents must be a pandas DataFrame.", str(context.exception))

    def test_experiment_info_table_to_json_default(self):
        df = pd.DataFrame({
            'experiment_sample_name': ['sample1_MH_run1', 'sample2_MH_run1'],
            'sequencing_info_name': ['run1', 'run1'],
            'specimen_name': ['sample1', 'sample2'],
            'panel_name': ['MH', 'MH'],
        })

        result = experiment_info_table_to_json(df)
        self.assertEqual([{'experiment_sample_name': 'sample1_MH_run1',
                           'sequencing_info_name': 'run1',
                           'specimen_name': 'sample1',
                           'panel_name': 'MH'},
                          {'experiment_sample_name': 'sample2_MH_run1',
                           'sequencing_info_name': 'run1',
                           'specimen_name': 'sample2',
                           'panel_name': 'MH'}], result)

    def test_experiment_info_table_to_json_with_plate(self):
        df = pd.DataFrame({
            'experiment_sample_name': ['sample1_MH_run1', 'sample2_MH_run1'],
            'sequencing_info_name': ['run1', 'run1'],
            'specimen_name': ['sample1', 'sample2'],
            'panel_name': ['MH', 'MH'],
            'plate_col': [1, 2],
            'plate_name': ['plate1', 'plate1'],
            'plate_row': ['A', 'B']
        })

        result = experiment_info_table_to_json(df, extraction_plate_name_col='plate_name', extraction_plate_col_col='plate_col', extraction_plate_row_col='plate_row',
                                               sequencing_plate_name_col='plate_name', sequencing_plate_col_col='plate_col', sequencing_plate_row_col='plate_row')
        self.assertEqual([{'experiment_sample_name': 'sample1_MH_run1',
                           'sequencing_info_name': 'run1',
                           'specimen_name': 'sample1',
                           'panel_name': 'MH',
                           'extraction_plate_info': {
                               'plate_name': 'plate1',
                               'plate_row': 'A',
                               'plate_col': 1
                           },
                           'sequencing_prep_plate_info': {
                               'plate_name': 'plate1',
                               'plate_row': 'A',
                               'plate_col': 1
                           }
                           },
                          {
                         'experiment_sample_name': 'sample2_MH_run1',
                         'sequencing_info_name': 'run1',
                         'specimen_name': 'sample2',
                         'panel_name': 'MH',
                         'extraction_plate_info': {
                             'plate_name': 'plate1',
                             'plate_row': 'B',
                             'plate_col': 2
                         },
                         'sequencing_prep_plate_info': {
                             'plate_name': 'plate1',
                             'plate_row': 'B',
                             'plate_col': 2
                         }}], result)

    def test_experiment_info_table_to_json_with_additional_columns(self):
        df = pd.DataFrame({
            'experiment_sample_name': ['sample1_MH_run1', 'sample2_MH_run1'],
            'sequencing_info_name': ['run1', 'run1'],
            'specimen_name': ['sample1', 'sample2'],
            'panel_name': ['MH', 'MH'],
            'new_col1': ['test', 'this'],
            'new_col2': ['add', 'one'],
        })

        result = experiment_info_table_to_json(
            df, additional_experiment_info_cols=['new_col1', 'new_col2'])
        self.assertEqual([{'experiment_sample_name': 'sample1_MH_run1',
                           'sequencing_info_name': 'run1',
                           'specimen_name': 'sample1',
                           'panel_name': 'MH',
                           'new_col1': 'test',
                           'new_col2': 'add'},
                          {'experiment_sample_name': 'sample2_MH_run1',
                           'sequencing_info_name': 'run1',
                           'specimen_name': 'sample2',
                           'panel_name': 'MH',
                           'new_col1': 'this',
                           'new_col2': 'one'}], result)

    def test_experiment_info_table_to_json_fails_with_duplicate_cols(self):
        df = pd.DataFrame({
            'experiment_sample_name': ['sample1_MH_run1', 'sample2_MH_run1'],
            'sequencing_info_name': ['run1', 'run1'],
            'specimen_name': ['sample1', 'sample2'],
            'panel_name': ['MH', 'MH'],
        })

        with self.assertRaises(ValueError) as context:
            experiment_info_table_to_json(
                df, specimen_name_col='panel_name')
        self.assertEqual(
            "Selected columns must be unique.", str(context.exception))

    def test_experiment_info_table_to_json_fails_with_missing_cols(self):
        df = pd.DataFrame({
            'experiment_sample_name': ['sample1_MH_run1', 'sample2_MH_run1'],
            'sequencing_info_name': ['run1', 'run1'],
        })

        with self.assertRaises(ValueError) as context:
            experiment_info_table_to_json(df)
        self.assertEqual(
            "The following columns are not in the DataFrame: ['specimen_name', 'panel_name']", str(context.exception))

    def test_experiment_info_table_to_json_fails_without_df(self):
        with self.assertRaises(ValueError) as context:
            specimen_info_table_to_json('test')
        self.assertEqual(
            "contents must be a pandas DataFrame.", str(context.exception))


if __name__ == '__main__':
    unittest.main()
