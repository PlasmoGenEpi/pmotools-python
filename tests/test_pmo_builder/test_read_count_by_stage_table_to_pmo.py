import unittest
import pandas as pd
import numpy as np

from pmotools.pmo_builder.read_count_by_stage_table_to_pmo import (
    read_count_by_stage_table_to_pmo,
    _process_total_raw_count_table,
    _process_reads_by_stage_table,
    _build_read_counts_by_stage_output,
)


class TestReadCountByStageTableToPMO(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

        # Sample total raw count table
        self.total_raw_count_table = pd.DataFrame(
            {
                "library_sample_name": ["sample1", "sample2", "sample3"],
                "total_raw_count": [1000, 2000, 1500],
                "additional_info": ["info1", "info2", "info3"],
            }
        )

        # Sample reads by stage table (long format)
        self.reads_by_stage_table_long = pd.DataFrame(
            {
                "library_sample_name": [
                    "sample1",
                    "sample1",
                    "sample1",
                    "sample2",
                    "sample2",
                    "sample2",
                ],
                "target_name": [
                    "target1",
                    "target1",
                    "target2",
                    "target1",
                    "target1",
                    "target2",
                ],
                "stage": [
                    "demultiplexed",
                    "denoised",
                    "demultiplexed",
                    "demultiplexed",
                    "denoised",
                    "demultiplexed",
                ],
                "reads": [100, 80, 50, 200, 150, 75],
            }
        )

        # Sample reads by stage table (wide format)
        self.reads_by_stage_table_wide = pd.DataFrame(
            {
                "library_sample_name": ["sample1", "sample2"],
                "target_name": ["target1", "target1"],
                "demultiplexed": [100, 200],
                "denoised": [80, 150],
                "filtered": [60, 120],
            }
        )

        # Expected output structure
        self.expected_output_structure = {
            "bioinformatics_run_name": "test_run",
            "read_counts_by_library_sample_by_stage": [
                {
                    "library_sample_name": "sample1",
                    "total_raw_count": 1000,
                    "read_counts_for_targets": [
                        {
                            "target_name": "target1",
                            "stages": [
                                {"stage": "demultiplexed", "reads": 100},
                                {"stage": "denoised", "reads": 80},
                            ],
                        }
                    ],
                }
            ],
        }

    def test_basic_functionality_long_format(self):
        """Test basic functionality with long format reads by stage table."""
        result = read_count_by_stage_table_to_pmo(
            total_raw_count_table=self.total_raw_count_table,
            bioinformatics_run_name="test_run",
            read_count_col="reads",
            reads_by_stage_table=self.reads_by_stage_table_long,
        )

        # Should return a list with one entry
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)

        run_data = result[0]
        self.assertIn("bioinformatics_run_name", run_data)
        self.assertIn("read_counts_by_library_sample_by_stage", run_data)
        self.assertEqual(run_data["bioinformatics_run_name"], "test_run")
        self.assertEqual(len(run_data["read_counts_by_library_sample_by_stage"]), 3)

    def test_basic_functionality_wide_format(self):
        """Test basic functionality with wide format reads by stage table."""
        result = read_count_by_stage_table_to_pmo(
            total_raw_count_table=self.total_raw_count_table,
            bioinformatics_run_name="test_run",
            reads_by_stage_table=self.reads_by_stage_table_wide,
            stage_col=["demultiplexed", "denoised", "filtered"],
        )

        # Should return a list with one entry
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)

        run_data = result[0]
        self.assertIn("bioinformatics_run_name", run_data)
        self.assertIn("read_counts_by_library_sample_by_stage", run_data)
        self.assertEqual(run_data["bioinformatics_run_name"], "test_run")

    def test_only_total_raw_count_table(self):
        """Test functionality with only total raw count table (no reads by stage)."""
        result = read_count_by_stage_table_to_pmo(
            total_raw_count_table=self.total_raw_count_table,
            bioinformatics_run_name="test_run",
        )

        # Should return a list with one entry
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)

        run_data = result[0]
        self.assertEqual(run_data["bioinformatics_run_name"], "test_run")
        self.assertEqual(len(run_data["read_counts_by_library_sample_by_stage"]), 3)

        # Check that each sample has total_raw_count but no read_counts_for_targets
        for sample in run_data["read_counts_by_library_sample_by_stage"]:
            self.assertIn("total_raw_count", sample)
            self.assertNotIn("read_counts_for_targets", sample)

    def test_custom_column_names(self):
        """Test with custom column names."""
        custom_total_table = pd.DataFrame(
            {"sample_id": ["sample1", "sample2"], "raw_count": [1000, 2000]}
        )

        custom_reads_table = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample1"],
                "locus": ["target1", "target1"],
                "pipeline_stage": ["demultiplexed", "denoised"],
                "counts": [100, 80],
            }
        )

        result = read_count_by_stage_table_to_pmo(
            total_raw_count_table=custom_total_table,
            bioinformatics_run_name="test_run",
            reads_by_stage_table=custom_reads_table,
            library_sample_name_col="sample_id",
            target_name_col="locus",
            total_raw_count_col="raw_count",
            stage_col="pipeline_stage",
            read_count_col="counts",
        )

        # Should return a list with one entry
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["bioinformatics_run_name"], "test_run")

    def test_additional_columns(self):
        """Test with additional columns."""
        result = read_count_by_stage_table_to_pmo(
            total_raw_count_table=self.total_raw_count_table,
            bioinformatics_run_name="test_run",
            additional_library_sample_cols=["additional_info"],
        )

        # Should return a list with one entry
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)

        run_data = result[0]
        # Check that additional columns are included
        for sample in run_data["read_counts_by_library_sample_by_stage"]:
            if sample["library_sample_name"] == "sample1":
                self.assertIn("additional_info", sample)
                self.assertEqual(sample["additional_info"], "info1")

    def test_additional_target_columns(self):
        """Test with additional target columns in reads by stage table."""
        # Create reads table with additional columns
        reads_table_with_additional = pd.DataFrame(
            {
                "library_sample_name": ["sample1", "sample1", "sample2"],
                "target_name": ["target1", "target1", "target1"],
                "stage": ["demultiplexed", "denoised", "demultiplexed"],
                "read_count": [100, 80, 200],
                "quality_score": [0.95, 0.98, 0.92],
                "coverage_depth": [50, 40, 100],
            }
        )

        result = read_count_by_stage_table_to_pmo(
            total_raw_count_table=self.total_raw_count_table,
            bioinformatics_run_name="test_run",
            reads_by_stage_table=reads_table_with_additional,
            additional_target_cols=["quality_score", "coverage_depth"],
        )

        # Should return a list with one entry
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)

        run_data = result[0]

        # Check that additional columns are included in stages
        sample1_data = next(
            s
            for s in run_data["read_counts_by_library_sample_by_stage"]
            if s["library_sample_name"] == "sample1"
        )

        self.assertIn("read_counts_for_targets", sample1_data)
        target1_data = sample1_data["read_counts_for_targets"][0]
        self.assertEqual(target1_data["target_name"], "target1")

        # Check that stages contain additional columns
        stages = target1_data["stages"]
        demultiplexed_stage = next(s for s in stages if s["stage"] == "demultiplexed")
        denoised_stage = next(s for s in stages if s["stage"] == "denoised")

        # Check demultiplexed stage has additional columns
        self.assertIn("quality_score", demultiplexed_stage)
        self.assertIn("coverage_depth", demultiplexed_stage)
        self.assertEqual(demultiplexed_stage["quality_score"], 0.95)
        self.assertEqual(demultiplexed_stage["coverage_depth"], 50)

        # Check denoised stage has additional columns
        self.assertIn("quality_score", denoised_stage)
        self.assertIn("coverage_depth", denoised_stage)
        self.assertEqual(denoised_stage["quality_score"], 0.98)
        self.assertEqual(denoised_stage["coverage_depth"], 40)

    def test_wide_format_melt_conversion(self):
        """Test that wide format is properly converted using pd.melt."""
        result = read_count_by_stage_table_to_pmo(
            total_raw_count_table=self.total_raw_count_table,
            bioinformatics_run_name="test_run",
            reads_by_stage_table=self.reads_by_stage_table_wide,
            stage_col=["demultiplexed", "denoised", "filtered"],
        )

        # Should return a list with one entry
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)

        run_data = result[0]
        # Check that the wide format was converted properly
        sample1_data = next(
            s
            for s in run_data["read_counts_by_library_sample_by_stage"]
            if s["library_sample_name"] == "sample1"
        )

        self.assertIn("read_counts_for_targets", sample1_data)
        target1_data = sample1_data["read_counts_for_targets"][0]
        self.assertEqual(target1_data["target_name"], "target1")

        # Check that all stages are present
        stage_names = [stage["stage"] for stage in target1_data["stages"]]
        self.assertIn("demultiplexed", stage_names)
        self.assertIn("denoised", stage_names)
        self.assertIn("filtered", stage_names)

    def test_input_validation(self):
        """Test input validation."""
        # Test non-DataFrame input for total_raw_count_table
        with self.assertRaises(ValueError):
            read_count_by_stage_table_to_pmo(
                total_raw_count_table="not_a_dataframe",
                bioinformatics_run_name="test_run",
            )

        # Test non-DataFrame input for reads_by_stage_table
        with self.assertRaises(ValueError):
            read_count_by_stage_table_to_pmo(
                total_raw_count_table=self.total_raw_count_table,
                bioinformatics_run_name="test_run",
                reads_by_stage_table="not_a_dataframe",
            )

    def test_missing_columns(self):
        """Test error handling for missing columns."""
        # Test missing required columns in total_raw_count_table
        bad_table = pd.DataFrame(
            {"wrong_col": ["sample1", "sample2"], "total_raw_count": [1000, 2000]}
        )

        with self.assertRaises(ValueError):
            read_count_by_stage_table_to_pmo(
                total_raw_count_table=bad_table, bioinformatics_run_name="test_run"
            )

        # Test missing required columns in reads_by_stage_table
        bad_reads_table = pd.DataFrame(
            {
                "library_sample_name": ["sample1"],
                "wrong_col": ["target1"],
                "stage": ["demultiplexed"],
                "read_count": [100],
            }
        )

        with self.assertRaises(ValueError):
            read_count_by_stage_table_to_pmo(
                total_raw_count_table=self.total_raw_count_table,
                bioinformatics_run_name="test_run",
                reads_by_stage_table=bad_reads_table,
            )

    def test_duplicate_library_samples(self):
        """Test error handling for duplicate library sample names."""
        duplicate_table = pd.DataFrame(
            {
                "library_sample_name": ["sample1", "sample1"],
                "total_raw_count": [1000, 2000],
            }
        )

        with self.assertRaises(ValueError):
            read_count_by_stage_table_to_pmo(
                total_raw_count_table=duplicate_table,
                bioinformatics_run_name="test_run",
            )

    def test_additional_columns_validation(self):
        """Test validation of additional columns."""
        with self.assertRaises(ValueError):
            read_count_by_stage_table_to_pmo(
                total_raw_count_table=self.total_raw_count_table,
                bioinformatics_run_name="test_run",
                additional_library_sample_cols=["nonexistent_column"],
            )

    def test_process_total_raw_count_table(self):
        """Test the _process_total_raw_count_table helper function."""
        result = _process_total_raw_count_table(
            self.total_raw_count_table,
            "library_sample_name",
            "total_raw_count",
            ["additional_info"],
        )

        self.assertIn("sample1", result)
        self.assertEqual(result["sample1"]["total_raw_count"], 1000)
        self.assertEqual(result["sample1"]["additional_info"], "info1")

    def test_process_reads_by_stage_table_long(self):
        """Test the _process_reads_by_stage_table helper function with long format."""
        result = _process_reads_by_stage_table(
            self.reads_by_stage_table_long,
            "library_sample_name",
            "target_name",
            "stage",
            "reads",
        )

        self.assertIn("sample1", result)
        self.assertIn("target1", result["sample1"])
        self.assertIn("demultiplexed", result["sample1"]["target1"])
        # Now stages are dictionaries with stage and read_count
        demultiplexed_data = result["sample1"]["target1"]["demultiplexed"]
        self.assertEqual(demultiplexed_data["stage"], "demultiplexed")
        self.assertEqual(demultiplexed_data["reads"], 100)

    def test_process_reads_by_stage_table_wide(self):
        """Test the _process_reads_by_stage_table helper function with wide format."""
        result = _process_reads_by_stage_table(
            self.reads_by_stage_table_wide,
            "library_sample_name",
            "target_name",
            ["demultiplexed", "denoised", "filtered"],
            "reads",
        )

        self.assertIn("sample1", result)
        self.assertIn("target1", result["sample1"])
        self.assertIn("demultiplexed", result["sample1"]["target1"])
        # Now stages are dictionaries with stage and read_count
        demultiplexed_data = result["sample1"]["target1"]["demultiplexed"]
        self.assertEqual(demultiplexed_data["stage"], "demultiplexed")
        self.assertEqual(demultiplexed_data["reads"], 100)

    def test_build_read_counts_by_stage_output(self):
        """Test the _build_read_counts_by_stage_output helper function."""
        library_sample_data = {
            "sample1": {"total_raw_count": 1000},
            "sample2": {"total_raw_count": 2000},
        }

        reads_by_stage_data = {
            "sample1": {"target1": {"demultiplexed": 100, "denoised": 80}}
        }

        result = _build_read_counts_by_stage_output(
            library_sample_data, reads_by_stage_data, "test_run"
        )

        self.assertEqual(result["bioinformatics_run_name"], "test_run")
        self.assertEqual(len(result["read_counts_by_library_sample_by_stage"]), 2)

        # Check sample1 has read_counts_for_targets
        sample1 = next(
            s
            for s in result["read_counts_by_library_sample_by_stage"]
            if s["library_sample_name"] == "sample1"
        )
        self.assertIn("read_counts_for_targets", sample1)

        # Check sample2 does not have read_counts_for_targets
        sample2 = next(
            s
            for s in result["read_counts_by_library_sample_by_stage"]
            if s["library_sample_name"] == "sample2"
        )
        self.assertNotIn("read_counts_for_targets", sample2)

    def test_empty_reads_by_stage_data(self):
        """Test behavior when reads_by_stage_data is None."""
        library_sample_data = {"sample1": {"total_raw_count": 1000}}

        result = _build_read_counts_by_stage_output(
            library_sample_data, None, "test_run"
        )

        sample1 = result["read_counts_by_library_sample_by_stage"][0]
        self.assertNotIn("read_counts_for_targets", sample1)

    def test_numeric_read_counts(self):
        """Test that read counts are properly converted to integers."""
        table_with_float_counts = pd.DataFrame(
            {
                "library_sample_name": ["sample1"],
                "target_name": ["target1"],
                "stage": ["demultiplexed"],
                "read_count": [100.5],  # Float value
            }
        )

        result = _process_reads_by_stage_table(
            table_with_float_counts,
            "library_sample_name",
            "target_name",
            "stage",
            "read_count",
        )

        # Should be converted to int
        demultiplexed_data = result["sample1"]["target1"]["demultiplexed"]
        self.assertEqual(demultiplexed_data["reads"], 100)

    def test_missing_values_handling(self):
        """Test handling of missing values in additional columns."""
        table_with_nan = pd.DataFrame(
            {
                "library_sample_name": ["sample1", "sample2"],
                "total_raw_count": [1000, 2000],
                "optional_info": ["info1", np.nan],
            }
        )

        result = _process_total_raw_count_table(
            table_with_nan, "library_sample_name", "total_raw_count", ["optional_info"]
        )

        # sample1 should have optional_info, sample2 should not
        self.assertIn("optional_info", result["sample1"])
        self.assertNotIn("optional_info", result["sample2"])

    def test_multiple_runs_from_column(self):
        """Test functionality when bioinformatics_run_name is a column in total_raw_count_table."""
        # Create table with multiple runs
        multi_run_table = pd.DataFrame(
            {
                "library_sample_name": ["sample1", "sample2", "sample1", "sample2"],
                "total_raw_count": [1000, 2000, 1500, 2500],
                "bioinformatics_run_name": ["run1", "run1", "run2", "run2"],
            }
        )

        # Create reads table with multiple runs
        multi_run_reads = pd.DataFrame(
            {
                "library_sample_name": ["sample1", "sample1", "sample2", "sample2"],
                "target_name": ["target1", "target1", "target1", "target1"],
                "stage": ["demultiplexed", "denoised", "demultiplexed", "denoised"],
                "read_count": [100, 80, 200, 150],
                "bioinformatics_run_name": ["run1", "run1", "run1", "run2"],
            }
        )

        result = read_count_by_stage_table_to_pmo(
            total_raw_count_table=multi_run_table,
            bioinformatics_run_name="bioinformatics_run_name",
            reads_by_stage_table=multi_run_reads,
        )

        # Should return a list of dictionaries
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)  # Two runs

        # Check run1
        run1_data = next(r for r in result if r["bioinformatics_run_name"] == "run1")
        self.assertEqual(
            len(run1_data["read_counts_by_library_sample_by_stage"]), 2
        )  # sample1 and sample2

        # Check run2
        run2_data = next(r for r in result if r["bioinformatics_run_name"] == "run2")
        self.assertEqual(
            len(run2_data["read_counts_by_library_sample_by_stage"]), 2
        )  # sample1 and sample2

        # Verify that samples are filtered correctly for each run
        run1_samples = [
            s["library_sample_name"]
            for s in run1_data["read_counts_by_library_sample_by_stage"]
        ]
        run2_samples = [
            s["library_sample_name"]
            for s in run2_data["read_counts_by_library_sample_by_stage"]
        ]

        self.assertIn("sample1", run1_samples)
        self.assertIn("sample2", run1_samples)
        self.assertIn("sample1", run2_samples)
        self.assertIn("sample2", run2_samples)

    def test_multiple_runs_without_reads_table(self):
        """Test multiple runs functionality without reads by stage table."""
        multi_run_table = pd.DataFrame(
            {
                "library_sample_name": ["sample1", "sample2", "sample1", "sample2"],
                "total_raw_count": [1000, 2000, 1500, 2500],
                "bioinformatics_run_name": ["run1", "run1", "run2", "run2"],
            }
        )

        result = read_count_by_stage_table_to_pmo(
            total_raw_count_table=multi_run_table,
            bioinformatics_run_name="bioinformatics_run_name",
        )

        # Should return a list of dictionaries
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)  # Two runs

        # Check that each run has the correct samples
        run1_data = next(r for r in result if r["bioinformatics_run_name"] == "run1")
        run2_data = next(r for r in result if r["bioinformatics_run_name"] == "run2")

        self.assertEqual(len(run1_data["read_counts_by_library_sample_by_stage"]), 2)
        self.assertEqual(len(run2_data["read_counts_by_library_sample_by_stage"]), 2)

    def test_multiple_runs_reads_table_without_run_column(self):
        """Test multiple runs when reads table doesn't have bioinformatics_run_name column."""
        multi_run_table = pd.DataFrame(
            {
                "library_sample_name": ["sample1", "sample2", "sample1", "sample2"],
                "total_raw_count": [1000, 2000, 1500, 2500],
                "bioinformatics_run_name": ["run1", "run1", "run2", "run2"],
            }
        )

        # Reads table without bioinformatics_run_name column
        reads_table = pd.DataFrame(
            {
                "library_sample_name": ["sample1", "sample2"],
                "target_name": ["target1", "target1"],
                "stage": ["demultiplexed", "demultiplexed"],
                "read_count": [100, 200],
            }
        )

        result = read_count_by_stage_table_to_pmo(
            total_raw_count_table=multi_run_table,
            bioinformatics_run_name="bioinformatics_run_name",
            reads_by_stage_table=reads_table,
        )

        # Should return a list of dictionaries
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)  # Two runs

        # Both runs should have the same reads data (since reads table doesn't have run column)
        run1_data = next(r for r in result if r["bioinformatics_run_name"] == "run1")
        run2_data = next(r for r in result if r["bioinformatics_run_name"] == "run2")

        # Both should have read_counts_for_targets since reads table applies to all runs
        for run_data in [run1_data, run2_data]:
            for sample in run_data["read_counts_by_library_sample_by_stage"]:
                if sample["library_sample_name"] in ["sample1", "sample2"]:
                    self.assertIn("read_counts_for_targets", sample)

    def test_single_run_behavior_unchanged(self):
        """Test that single run behavior returns a list with one entry."""
        result = read_count_by_stage_table_to_pmo(
            total_raw_count_table=self.total_raw_count_table,
            bioinformatics_run_name="test_run",
            reads_by_stage_table=self.reads_by_stage_table_long,
            read_count_col="reads",
        )

        # Should return a list with one entry
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["bioinformatics_run_name"], "test_run")


if __name__ == "__main__":
    unittest.main()
