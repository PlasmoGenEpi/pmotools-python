import unittest
import pandas as pd
import copy
from unittest.mock import patch

from pmotools.pmo_builder.mhap_table_to_pmo_json import *


class TestMhapTableToJson(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None
        self.small_representative_dict = {"targets": [{"target_name": 'target1', "microhaplotypes": [{"seq": "ACTG"}, {"seq": "ATTG"}, {"seq": "ATTA"}]}, {
            "target_name": 'target2', "microhaplotypes": [{"seq": "TTTT"}, {"seq": "TTTA"},]}, {"target_name": 'target3', "microhaplotypes": [{"seq": "GGG"}, {"seq": "AAG"}]}]}
        self.small_detected_dict = {"bioinformatics_run_name": 'run1', "library_samples": [{"library_sample_name": "sample1", "target_results": [{"mhaps_target_id": 0, "mhaps": [{"mhap_id": 0, "reads": 10}]}, {"mhaps_target_id": 1, "mhaps": [{"mhap_id": 0, "reads": 100}]}, {"mhaps_target_id": 2, "mhaps": [{"mhap_id": 0, "reads": 10}]}]},
                                                                                           {"library_sample_name": "sample2", "target_results": [{"mhaps_target_id": 0, "mhaps": [{"mhap_id": 0, "reads": 1}]}, {"mhaps_target_id": 1, "mhaps": [
                                                                                               {"mhap_id": 0, "reads": 12}, {"mhap_id": 1, "reads": 53}]}, {"mhaps_target_id": 2, "mhaps": [{"mhap_id": 1, "reads": 73}]}]},
                                                                                           {"library_sample_name": "sample3", "target_results": [{"mhaps_target_id": 0, "mhaps": [
                                                                                               {"mhap_id": 1, "reads": 76}]}, {"mhaps_target_id": 1, "mhaps": [{"mhap_id": 1, "reads": 6}]},]},
                                                                                           {"library_sample_name": "sample4", "target_results": [
                                                                                               {"mhaps_target_id": 0, "mhaps": [{"mhap_id": 0, "reads": 297}]}]},
                                                                                           {"library_sample_name": "sample5", "target_results": [{"mhaps_target_id": 0, "mhaps": [{"mhap_id": 1, "reads": 123}, {
                                                                                               "mhap_id": 2, "reads": 1}]}, {"mhaps_target_id": 2, "mhaps": [{"mhap_id": 0, "reads": 17}]}]},
                                                                                           ]}
        small_mhap_data = {'library_sample_name': ['sample1', 'sample1', 'sample1', 'sample2', 'sample2', 'sample2', 'sample2', 'sample3', 'sample3', 'sample4', 'sample5', 'sample5', 'sample5'],
                           'target_name':            ['target1', 'target2', 'target3', 'target1', 'target2', 'target2', 'target3', 'target1', 'target2', 'target1', 'target1', 'target1', 'target3'],
                           'seq': ['ACTG', 'TTTT', 'GGG', 'ACTG', 'TTTT', 'TTTA', 'AAG', 'ATTG', 'TTTA', 'ACTG', 'ATTG', 'ATTA', 'GGG'],
                           'reads': [10, 100, 10, 1, 12, 53, 73, 76, 6, 297, 123, 1, 17], }
        self.small_mhap_table = pd.DataFrame(data=small_mhap_data)
        self.small_df_mhaps_target_id_values = [
            0, 1, 2, 0, 1, 1, 2, 0, 1, 0, 0, 0, 2]
        self.small_df_mhap_id_values = [0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 2, 0]
        self.df_with_mhap_target_id = self.small_mhap_table.copy()
        self.df_with_mhap_target_id['mhaps_target_id'] = self.small_df_mhaps_target_id_values
        self.df_with_mhap_id = self.df_with_mhap_target_id.copy()
        self.df_with_mhap_id['mhap_id'] = self.small_df_mhap_id_values

        self.bioinformatics_run_name = 'run1'

    def test_get_target_id_in_representative_mhaps(self):
        expected_cols = self.small_mhap_table.columns.to_list() + \
            ['mhaps_target_id']
        actual = get_target_id_in_representative_mhaps(
            self.small_mhap_table.copy(), self.small_representative_dict)
        self.assertListEqual(actual.columns.to_list(), expected_cols)
        self.assertListEqual(
            actual.mhaps_target_id.to_list(), self.small_df_mhaps_target_id_values)

    def test_get_target_id_in_representative_mhaps_missing_target(self):
        df_with_missing = self.small_mhap_table.replace(
            'target1', 'missing_target').copy()
        with self.assertRaises(ValueError) as context:
            get_target_id_in_representative_mhaps(
                df_with_missing, self.small_representative_dict)
        self.assertEqual(
            "Missing target_name(s) in representative microhaplotype table: ['missing_target']", str(context.exception))

    def test_get_mhap_index_in_representative_mhaps(self):
        actual = get_mhap_index_in_representative_mhaps(
            self.df_with_mhap_target_id, self.small_representative_dict)
        expected_cols = self.small_mhap_table.columns.to_list() + \
            ['mhaps_target_id', 'mhap_id']
        self.assertListEqual(actual.columns.to_list(), expected_cols)
        self.assertListEqual(
            actual.mhap_id.to_list(), self.small_df_mhap_id_values)

    def test_get_mhap_index_in_representative_mhaps_missing_mhap(self):
        df_with_missing = self.small_mhap_table.replace(
            'ATTG', 'missing_seq1').copy()
        df_with_missing.replace('TTTT', 'missing_seq2', inplace=True)
        df_with_missing['mhaps_target_id'] = self.small_df_mhaps_target_id_values
        with self.assertRaises(ValueError) as context:
            get_mhap_index_in_representative_mhaps(
                df_with_missing, self.small_representative_dict)
        expected_missing = pd.DataFrame({'target_name': ['target2', 'target1'], 'seq': [
                                        'missing_seq2', 'missing_seq1']}, index=[1, 7])
        self.assertIn("Some seq values not found", str(context.exception))
        self.assertIn("missing_seq1", str(context.exception))
        self.assertIn("missing_seq2", str(context.exception))

    def test_build_detected_mhap_dict_minimum(self):
        mhap_cols = ['mhap_id', 'reads']
        actual = build_detected_mhap_dict(
            self.df_with_mhap_id, self.bioinformatics_run_name, mhap_cols)
        self.assertDictEqual(actual, self.small_detected_dict)

    def test_build_detected_mhap_dict_with_umis(self):
        mhap_cols = ['mhap_id', 'reads', 'umis']
        df_with_umis = self.df_with_mhap_id.copy()
        df_with_umis['umis'] = [10, 100, 10,
                                1, 12, 53, 73, 76, 6, 297, 123, 1, 17]
        # sample index, target index, mhap index
        small_detected_dict_with_umis = copy.deepcopy(
            self.small_detected_dict)
        small_detected_dict_with_umis["library_samples"][0]["target_results"][0]["mhaps"][0]["umis"] = 10
        small_detected_dict_with_umis["library_samples"][0]["target_results"][1]["mhaps"][0]["umis"] = 100
        small_detected_dict_with_umis["library_samples"][0]["target_results"][2]["mhaps"][0]["umis"] = 10
        small_detected_dict_with_umis["library_samples"][1]["target_results"][0]["mhaps"][0]["umis"] = 1
        small_detected_dict_with_umis["library_samples"][1]["target_results"][1]["mhaps"][0]["umis"] = 12
        small_detected_dict_with_umis["library_samples"][1]["target_results"][1]["mhaps"][1]["umis"] = 53
        small_detected_dict_with_umis["library_samples"][1]["target_results"][2]["mhaps"][0]["umis"] = 73
        small_detected_dict_with_umis["library_samples"][2]["target_results"][0]["mhaps"][0]["umis"] = 76
        small_detected_dict_with_umis["library_samples"][2]["target_results"][1]["mhaps"][0]["umis"] = 6
        small_detected_dict_with_umis["library_samples"][3]["target_results"][0]["mhaps"][0]["umis"] = 297
        small_detected_dict_with_umis["library_samples"][4]["target_results"][0]["mhaps"][0]["umis"] = 123
        small_detected_dict_with_umis["library_samples"][4]["target_results"][0]["mhaps"][1]["umis"] = 1
        small_detected_dict_with_umis["library_samples"][4]["target_results"][1]["mhaps"][0]["umis"] = 17
        # small_detected_dict_with_umis
        actual = build_detected_mhap_dict(
            df_with_umis, self.bioinformatics_run_name, mhap_cols)
        self.assertDictEqual(actual, small_detected_dict_with_umis)

    def test_build_detected_mhap_dict_with_additional_cols(self):
        mhap_cols = ['mhap_id', 'reads', 'ad_col1', 'ad_col2']

        small_detected_dict_with_ad_cols = copy.deepcopy(
            self.small_detected_dict)
        small_detected_dict_with_ad_cols["library_samples"][0]["target_results"][0]["mhaps"][0]["ad_col1"] = 'this1'
        small_detected_dict_with_ad_cols["library_samples"][0]["target_results"][1]["mhaps"][0]["ad_col1"] = 'this2'
        small_detected_dict_with_ad_cols["library_samples"][0]["target_results"][2]["mhaps"][0]["ad_col1"] = 'this3'
        small_detected_dict_with_ad_cols["library_samples"][1]["target_results"][0]["mhaps"][0]["ad_col1"] = 'this4'
        small_detected_dict_with_ad_cols["library_samples"][1]["target_results"][1]["mhaps"][0]["ad_col1"] = 'this1'
        small_detected_dict_with_ad_cols["library_samples"][1]["target_results"][1]["mhaps"][1]["ad_col1"] = 'this2'
        small_detected_dict_with_ad_cols["library_samples"][1]["target_results"][2]["mhaps"][0]["ad_col1"] = 'this3'
        small_detected_dict_with_ad_cols["library_samples"][2]["target_results"][0]["mhaps"][0]["ad_col1"] = 'this4'
        small_detected_dict_with_ad_cols["library_samples"][2]["target_results"][1]["mhaps"][0]["ad_col1"] = 'this1'
        small_detected_dict_with_ad_cols["library_samples"][3]["target_results"][0]["mhaps"][0]["ad_col1"] = 'this2'
        small_detected_dict_with_ad_cols["library_samples"][4]["target_results"][0]["mhaps"][0]["ad_col1"] = 'this3'
        small_detected_dict_with_ad_cols["library_samples"][4]["target_results"][0]["mhaps"][1]["ad_col1"] = 'this4'
        small_detected_dict_with_ad_cols["library_samples"][4]["target_results"][1]["mhaps"][0]["ad_col1"] = 'this4'
        small_detected_dict_with_ad_cols["library_samples"][0]["target_results"][1]["mhaps"][0]["ad_col2"] = 'that1'
        small_detected_dict_with_ad_cols["library_samples"][0]["target_results"][2]["mhaps"][0]["ad_col2"] = 'that2'
        small_detected_dict_with_ad_cols["library_samples"][1]["target_results"][1]["mhaps"][0]["ad_col2"] = 'that1'
        small_detected_dict_with_ad_cols["library_samples"][1]["target_results"][1]["mhaps"][1]["ad_col2"] = 'that2'
        small_detected_dict_with_ad_cols["library_samples"][2]["target_results"][0]["mhaps"][0]["ad_col2"] = 'that1'
        small_detected_dict_with_ad_cols["library_samples"][2]["target_results"][1]["mhaps"][0]["ad_col2"] = 'that2'
        small_detected_dict_with_ad_cols["library_samples"][4]["target_results"][0]["mhaps"][0]["ad_col2"] = 'that1'
        small_detected_dict_with_ad_cols["library_samples"][4]["target_results"][0]["mhaps"][1]["ad_col2"] = 'that2'

        df_with_ad_cols = self.df_with_mhap_id.copy()
        df_with_ad_cols['ad_col1'] = ['this1', 'this2', 'this3', 'this4', 'this1',
                                      'this2', 'this3', 'this4', 'this1', 'this2', 'this3', 'this4', 'this4']
        df_with_ad_cols['ad_col2'] = [None, 'that1', 'that2', None, 'that1',
                                      'that2', None, 'that1', 'that2', None, 'that1', 'that2', None]
        actual = build_detected_mhap_dict(
            df_with_ad_cols, self.bioinformatics_run_name, mhap_cols)
        self.assertDictEqual(actual, small_detected_dict_with_ad_cols)

    @patch("pmotools.pmo_builder.mhap_table_to_pmo_json.check_additional_columns_exist")
    @patch("pmotools.pmo_builder.mhap_table_to_pmo_json.get_target_id_in_representative_mhaps")
    @patch("pmotools.pmo_builder.mhap_table_to_pmo_json.get_mhap_index_in_representative_mhaps")
    @patch("pmotools.pmo_builder.mhap_table_to_pmo_json.build_detected_mhap_dict")
    def test_create_detected_microhaplotype_dict(self, mock_build_detected_mhap_dict, mock_get_mhap_index_in_representative_mhaps, mock_get_target_id_in_representative_mhaps, mock_check_additional_columns_exist):
        mock_get_target_id_in_representative_mhaps.return_value = self.small_df_mhaps_target_id_values
        mock_get_mhap_index_in_representative_mhaps.return_value = self.small_df_mhap_id_values
        mock_build_detected_mhap_dict.return_value = self.small_detected_dict
        actual = create_detected_microhaplotype_dict(
            self.small_mhap_table, self.bioinformatics_run_name, self.small_representative_dict)
        self.assertEqual(self.small_detected_dict, actual)

    def test_create_representative_microhaplotype_dict(self):
        actual = create_representative_microhaplotype_dict(
            self.small_mhap_table)
        self.assertDictEqual(self.small_representative_dict, actual)

    def test_create_representative_microhaplotype_dict_with_loc(self):
        # add locations to the df
        loc_df = pd.DataFrame({'target_name': ['target1', 'target2', 'target3'], 'chrom': [
            'chrom1', 'chrom2', 'chrom3'], 'start': [1, 2, 3], 'end': [4, 5, 6], 'ref': ['ACTG', 'TTTT', 'GGG']})
        mhap_table_with_loc = self.small_mhap_table.merge(loc_df)
        rep_dict_with_loc = copy.deepcopy(self.small_representative_dict)
        rep_dict_with_loc["targets"][0]['mhap_location'] = {
            'genome_id': 0, 'chrom': 'chrom1', 'start': 1, 'end': 4, 'ref_seq': 'ACTG'}
        rep_dict_with_loc["targets"][1]['mhap_location'] = {
            'genome_id': 0, 'chrom': 'chrom2', 'start': 2, 'end': 5, 'ref_seq': 'TTTT'}
        rep_dict_with_loc["targets"][2]['mhap_location'] = {
            'genome_id': 0, 'chrom': 'chrom3', 'start': 3, 'end': 6, 'ref_seq': 'GGG'}

        actual = create_representative_microhaplotype_dict(
            mhap_table_with_loc, chrom_col='chrom', start_col='start', end_col='end', ref_seq_col='ref')
        self.assertDictEqual(rep_dict_with_loc, actual)

    def test_create_representative_microhaplotype_dict_masking(self):
        masking_df = pd.DataFrame({'target_name': ['target1', 'target2'], 'masking_seq_start': ['5,92', '4'], 'masking_seq_segment_size': [
                                  '30,12', '6'], 'masking_replacement_size': ['30,13', '6']})
        mhap_table_with_masking = self.small_mhap_table.merge(
            masking_df, how='left')
        rep_dict_with_masking = copy.deepcopy(self.small_representative_dict)
        rep_dict_with_masking["targets"][0]['microhaplotypes'][0]["masking"] = [{'seq_start': 5, 'seq_segment_size': 30, 'replacement_size': 30}, {
            'seq_start': 92, 'seq_segment_size': 12, 'replacement_size': 13}]
        rep_dict_with_masking["targets"][0]['microhaplotypes'][1]["masking"] = [{'seq_start': 5, 'seq_segment_size': 30, 'replacement_size': 30}, {
            'seq_start': 92, 'seq_segment_size': 12, 'replacement_size': 13}]
        rep_dict_with_masking["targets"][0]['microhaplotypes'][2]["masking"] = [{'seq_start': 5, 'seq_segment_size': 30, 'replacement_size': 30}, {
            'seq_start': 92, 'seq_segment_size': 12, 'replacement_size': 13}]
        rep_dict_with_masking["targets"][1]['microhaplotypes'][0]["masking"] = [
            {'seq_start': 4, 'seq_segment_size': 6, 'replacement_size': 6}]
        rep_dict_with_masking["targets"][1]['microhaplotypes'][1]["masking"] = [
            {'seq_start': 4, 'seq_segment_size': 6, 'replacement_size': 6}]
        actual = create_representative_microhaplotype_dict(mhap_table_with_masking, masking_seq_start_col='masking_seq_start',
                                                           masking_seq_segment_size_col='masking_seq_segment_size', masking_replacement_size_col='masking_replacement_size')
        self.assertDictEqual(rep_dict_with_masking, actual)

    def test_create_representative_microhaplotype_dict_with_ad_cols(self):
        mhap_table_ad_cols = self.small_mhap_table.copy()
        mhap_table_ad_cols['adcol1'] = 'this'
        mhap_table_ad_cols['adcol2'] = 'that'
        mhap_table_ad_cols['cig'] = 'cigs'

        rep_dict_with_ad_cols = copy.deepcopy(self.small_representative_dict)
        rep_dict_with_ad_cols["targets"][0]['microhaplotypes'][0]["adcol1"] = 'this'
        rep_dict_with_ad_cols["targets"][0]['microhaplotypes'][0]["adcol2"] = 'that'
        rep_dict_with_ad_cols["targets"][0]['microhaplotypes'][0]["pseudo_cigar"] = 'cigs'
        rep_dict_with_ad_cols["targets"][0]['microhaplotypes'][1]["adcol1"] = 'this'
        rep_dict_with_ad_cols["targets"][0]['microhaplotypes'][1]["adcol2"] = 'that'
        rep_dict_with_ad_cols["targets"][0]['microhaplotypes'][1]["pseudo_cigar"] = 'cigs'
        rep_dict_with_ad_cols["targets"][0]['microhaplotypes'][2]["adcol1"] = 'this'
        rep_dict_with_ad_cols["targets"][0]['microhaplotypes'][2]["adcol2"] = 'that'
        rep_dict_with_ad_cols["targets"][0]['microhaplotypes'][2]["pseudo_cigar"] = 'cigs'
        rep_dict_with_ad_cols["targets"][1]['microhaplotypes'][0]["adcol1"] = 'this'
        rep_dict_with_ad_cols["targets"][1]['microhaplotypes'][0]["adcol2"] = 'that'
        rep_dict_with_ad_cols["targets"][1]['microhaplotypes'][0]["pseudo_cigar"] = 'cigs'
        rep_dict_with_ad_cols["targets"][1]['microhaplotypes'][1]["adcol1"] = 'this'
        rep_dict_with_ad_cols["targets"][1]['microhaplotypes'][1]["adcol2"] = 'that'
        rep_dict_with_ad_cols["targets"][1]['microhaplotypes'][1]["pseudo_cigar"] = 'cigs'
        rep_dict_with_ad_cols["targets"][2]['microhaplotypes'][0]["adcol1"] = 'this'
        rep_dict_with_ad_cols["targets"][2]['microhaplotypes'][0]["adcol2"] = 'that'
        rep_dict_with_ad_cols["targets"][2]['microhaplotypes'][0]["pseudo_cigar"] = 'cigs'
        rep_dict_with_ad_cols["targets"][2]['microhaplotypes'][1]["adcol1"] = 'this'
        rep_dict_with_ad_cols["targets"][2]['microhaplotypes'][1]["adcol2"] = 'that'
        rep_dict_with_ad_cols["targets"][2]['microhaplotypes'][1]["pseudo_cigar"] = 'cigs'

        actual = create_representative_microhaplotype_dict(
            mhap_table_ad_cols, additional_representative_mhap_cols=['adcol1', 'adcol2'], pseudocigar_col='cig')
        self.assertEqual(actual, rep_dict_with_ad_cols)

    @patch("pmotools.pmo_builder.mhap_table_to_pmo_json.create_representative_microhaplotype_dict")
    @patch("pmotools.pmo_builder.mhap_table_to_pmo_json.create_detected_microhaplotype_dict")
    def test_mhap_table_to_pmo_json(self, mock_create_detected_microhaplotype_dict, mock_create_representative_microhaplotype_dict):
        mock_create_detected_microhaplotype_dict.return_value = self.small_detected_dict
        mock_create_representative_microhaplotype_dict.return_value = self.small_representative_dict
        actual = mhap_table_to_pmo_json(
            self.small_mhap_table, self.bioinformatics_run_name)
        self.assertEqual(actual, {"representative_microhaplotypes": {"targets": [{"target_name": "target1", "microhaplotypes": [{"seq": "ACTG"}, {"seq": "ATTG"}, {"seq": "ATTA"}]}, {"target_name": "target2", "microhaplotypes": [{"seq": "TTTT"}, {"seq": "TTTA"}]}, {"target_name": "target3", "microhaplotypes": [{"seq": "GGG"}, {"seq": "AAG"}]}]}, "detected_microhaplotypes": [{"bioinformatics_run_name": "run1", "library_samples": [{"library_sample_name": "sample1", "target_results": [{"mhaps_target_id": 0, "mhaps": [{"mhap_id": 0, "reads": 10}]}, {"mhaps_target_id": 1, "mhaps": [{"mhap_id": 0, "reads": 100}]}, {"mhaps_target_id": 2, "mhaps": [{"mhap_id": 0, "reads": 10}]}]}, {"library_sample_name": "sample2", "target_results": [
                         {"mhaps_target_id": 0, "mhaps": [{"mhap_id": 0, "reads": 1}]}, {"mhaps_target_id": 1, "mhaps": [{"mhap_id": 0, "reads": 12}, {"mhap_id": 1, "reads": 53}]}, {"mhaps_target_id": 2, "mhaps": [{"mhap_id": 1, "reads": 73}]}]}, {"library_sample_name": "sample3", "target_results": [{"mhaps_target_id": 0, "mhaps": [{"mhap_id": 1, "reads": 76}]}, {"mhaps_target_id": 1, "mhaps": [{"mhap_id": 1, "reads": 6}]}]}, {"library_sample_name": "sample4", "target_results": [{"mhaps_target_id": 0, "mhaps": [{"mhap_id": 0, "reads": 297}]}]}, {"library_sample_name": "sample5", "target_results": [{"mhaps_target_id": 0, "mhaps": [{"mhap_id": 1, "reads": 123}, {"mhap_id": 2, "reads": 1}]}, {"mhaps_target_id": 2, "mhaps": [{"mhap_id": 0, "reads": 17}]}]}]}]})


if __name__ == '__main__':
    unittest.main()
