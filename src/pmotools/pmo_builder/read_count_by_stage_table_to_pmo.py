#!/usr/bin/env python3
import pandas as pd

from ..pmo_builder.json_convert_utils import check_additional_columns_exist


def read_count_by_stage_table_to_pmo(
    total_raw_count_table: pd.DataFrame,
    bioinformatics_run_name: str | None = None,
    reads_by_stage_table: pd.DataFrame | None = None,
    library_sample_name_col: str = "library_sample_name",
    target_name_col: str = "target_name",
    total_raw_count_col: str = "total_raw_count",
    stage_col: str | list = "stage",
    read_count_col: str = "read_count",
    additional_library_sample_cols: list | None = None,
    additional_target_cols: list | None = None,
) -> list[dict]:
    """
    Convert tables of read counts by stage into PMO read_counts_by_stage format.


    :param total_raw_count_table: table with total raw counts per sample
    :type total_raw_count_table: pd.DataFrame
    :param bioinformatics_run_name: name for the bioinformatics run (column name or individual run name)
    :type bioinformatics_run_name: str, optional
    :param reads_by_stage_table: table of reads per sample, per locus, per stage. Can be long
        format (single stage column) or wide format (multiple stage columns)
    :type reads_by_stage_table: pd.DataFrame, optional
    :param library_sample_name_col: column name for library sample names. Default: library_sample_name
    :type library_sample_name_col: str
    :param target_name_col: column name for target names. Default: target_name
    :type target_name_col: str
    :param total_raw_count_col: column name for total raw counts. Default: total_raw_count
    :type total_raw_count_col: str
    :param stage_col: column name for pipeline stages (long format), or list of stage column
        names for wide format. Default: stage
    :type stage_col: str or list
    :param read_count_col: column name for read counts. Default: read_count
    :type read_count_col: str
    :param additional_library_sample_cols: additional columns to include for library samples
    :type additional_library_sample_cols: list of str, optional
    :param additional_target_cols: additional columns to include for targets
    :type additional_target_cols: list of str, optional
    :return: list of dicts formatted for the PMO read_counts_by_stage section. Always returns
        a list, with one entry for single runs or multiple entries when bioinformatics_run_name
        is a column in total_raw_count_table.
    :rtype: list of dict

    .. note::
        For wide-format data, provide ``stage_col`` as a list of column names. The function
        uses ``pd.melt()`` to convert wide format to long format before processing.
    """

    # Validate input
    if not isinstance(total_raw_count_table, pd.DataFrame):
        raise ValueError("total_raw_count_table must be a pandas DataFrame.")

    if reads_by_stage_table is not None and not isinstance(
        reads_by_stage_table, pd.DataFrame
    ):
        raise ValueError("reads_by_stage_table must be a pandas DataFrame or None.")

    # Check additional columns if provided
    if additional_library_sample_cols:
        check_additional_columns_exist(
            total_raw_count_table, additional_library_sample_cols
        )

    if reads_by_stage_table is not None and additional_target_cols:
        check_additional_columns_exist(reads_by_stage_table, additional_target_cols)

    # Check if bioinformatics_run_name is a column in total_raw_count_table
    if (
        bioinformatics_run_name is not None
        and bioinformatics_run_name in total_raw_count_table.columns
    ):
        # Create separate entries for each unique run
        output_data_list = []
        unique_runs = total_raw_count_table[bioinformatics_run_name].unique()

        for run_name in unique_runs:
            # Filter data for this specific run
            run_total_table = total_raw_count_table[
                total_raw_count_table[bioinformatics_run_name] == run_name
            ].drop(columns=[bioinformatics_run_name])

            run_reads_table = None
            if reads_by_stage_table is not None:
                if bioinformatics_run_name in reads_by_stage_table.columns:
                    run_reads_table = reads_by_stage_table[
                        reads_by_stage_table[bioinformatics_run_name] == run_name
                    ].drop(columns=[bioinformatics_run_name])
                else:
                    # If reads_by_stage_table doesn't have bioinformatics_run_name column,
                    # use all data for all runs
                    run_reads_table = reads_by_stage_table

            # Process data for this run
            library_sample_data = _process_total_raw_count_table(
                run_total_table,
                library_sample_name_col,
                total_raw_count_col,
                additional_library_sample_cols,
            )

            reads_by_stage_data = None
            if run_reads_table is not None:
                reads_by_stage_data = _process_reads_by_stage_table(
                    run_reads_table,
                    library_sample_name_col,
                    target_name_col,
                    stage_col,
                    read_count_col,
                    additional_target_cols,
                )

            # Build output for this run
            run_output = _build_read_counts_by_stage_output(
                library_sample_data,
                reads_by_stage_data,
                run_name,
            )
            output_data_list.append(run_output)

        return output_data_list
    else:
        # Single run - still return as list for consistency
        library_sample_data = _process_total_raw_count_table(
            total_raw_count_table,
            library_sample_name_col,
            total_raw_count_col,
            additional_library_sample_cols,
        )

        reads_by_stage_data = None
        if reads_by_stage_table is not None:
            reads_by_stage_data = _process_reads_by_stage_table(
                reads_by_stage_table,
                library_sample_name_col,
                target_name_col,
                stage_col,
                read_count_col,
                additional_target_cols,
            )

        output_data = _build_read_counts_by_stage_output(
            library_sample_data,
            reads_by_stage_data,
            bioinformatics_run_name,
        )

        return [output_data]  # Return as list for consistency


def _process_total_raw_count_table(
    total_raw_count_table: pd.DataFrame,
    library_sample_name_col: str,
    total_raw_count_col: str,
    additional_cols: list | None = None,
) -> dict:
    """
    Process the total raw count table into a dictionary mapping library samples to their data.

    :param total_raw_count_table: DataFrame with total raw counts
    :param library_sample_name_col: Column name for library sample names
    :param total_raw_count_col: Column name for total raw counts
    :param additional_cols: Additional columns to include

    :return: Dictionary mapping library_sample_name to sample data
    """
    # Validate required columns exist
    required_cols = [library_sample_name_col, total_raw_count_col]
    missing_cols = [
        col for col in required_cols if col not in total_raw_count_table.columns
    ]
    if missing_cols:
        raise ValueError(
            f"Missing required columns in total_raw_count_table: {missing_cols}"
        )

    # Check for duplicates in library_sample_name_col
    if total_raw_count_table[library_sample_name_col].duplicated().any():
        duplicates = total_raw_count_table[
            total_raw_count_table[library_sample_name_col].duplicated(keep=False)
        ]
        raise ValueError(
            f"Duplicate library sample names found in total_raw_count_table:\n{duplicates}"
        )

    # Build the output dictionary
    sample_data = {}
    for _, row in total_raw_count_table.iterrows():
        sample_name = row[library_sample_name_col]
        total_raw_count = int(row[total_raw_count_col])

        sample_info = {
            "total_raw_count": total_raw_count,
        }

        # Add additional columns if specified
        if additional_cols:
            for col in additional_cols:
                if col in row and pd.notna(row[col]):
                    sample_info[col] = row[col]

        sample_data[sample_name] = sample_info

    return sample_data


def _process_reads_by_stage_table(
    reads_by_stage_table: pd.DataFrame,
    library_sample_name_col: str,
    target_name_col: str,
    stage_col: str | list,
    read_count_col: str,
    additional_cols: list | None = None,
) -> dict:
    """
    Process the reads by stage table into a nested dictionary structure.

    :param reads_by_stage_table: DataFrame with reads by stage information
    :param library_sample_name_col: Column name for library sample names
    :param target_name_col: Column name for target names
    :param stage_col: Column name for pipeline stages, or list of stage column names for wide format
    :param read_count_col: Column name for read counts
    :param additional_cols: Additional columns to include

    :return: Nested dictionary: {library_sample_name: {target_name: {stage: read_count}}}
    """
    # Handle wide format conversion if stage_col is a list
    if isinstance(stage_col, list):
        # Include additional columns in id_vars for wide format conversion
        id_vars = [library_sample_name_col, target_name_col]
        if additional_cols:
            id_vars.extend(additional_cols)

        # Convert wide format to long format using pd.melt
        processed_table = pd.melt(
            reads_by_stage_table,
            id_vars=id_vars,
            value_vars=stage_col,
            var_name="stage",
            value_name=read_count_col,
        )
        # Update stage_col to the new column name after melt
        stage_col = "stage"
    else:
        # Use the table as-is for long format
        processed_table = reads_by_stage_table.copy()

    # Validate required columns exist
    required_cols = [
        library_sample_name_col,
        target_name_col,
        stage_col,
        read_count_col,
    ]
    missing_cols = [col for col in required_cols if col not in processed_table.columns]
    if missing_cols:
        raise ValueError(
            f"Missing required columns in reads_by_stage_table: {missing_cols}"
        )

    # Build the nested dictionary structure
    reads_data = {}

    for _, row in processed_table.iterrows():
        sample_name = row[library_sample_name_col]
        target_name = row[target_name_col]
        stage = row[stage_col]
        read_count = int(row[read_count_col])

        # Initialize nested structure if needed
        if sample_name not in reads_data:
            reads_data[sample_name] = {}
        if target_name not in reads_data[sample_name]:
            reads_data[sample_name][target_name] = {}

        # Create stage data with read count and additional columns
        stage_data = {"stage": stage, "reads": read_count}

        # Add additional columns if present and not null
        if additional_cols:
            for col in additional_cols:
                if col in row and pd.notna(row[col]):
                    stage_data[col] = row[col]

        # Store the stage data
        reads_data[sample_name][target_name][stage] = stage_data

    return reads_data


def _build_read_counts_by_stage_output(
    library_sample_data: dict,
    reads_by_stage_data: dict | None,
    bioinformatics_run_name: str | None = None,
) -> dict:
    """
    Build the final output structure for read_counts_by_stage.

    :param library_sample_data: Dictionary with library sample data
    :param reads_by_stage_data: Optional dictionary with reads by stage data
    :param bioinformatics_run_name: Optional name for the bioinformatics run

    :return: Dictionary formatted for PMO read_counts_by_stage section
    """
    read_counts_by_library_sample_by_stage = []

    for sample_name, sample_info in library_sample_data.items():
        # Start building the library sample entry
        library_sample_entry = {
            "library_sample_name": sample_name,
            "total_raw_count": sample_info["total_raw_count"],
        }

        # Add additional fields from sample_info
        for key, value in sample_info.items():
            if key != "total_raw_count":
                library_sample_entry[key] = value

        # Process reads by stage data if available
        read_counts_for_targets = []
        if reads_by_stage_data and sample_name in reads_by_stage_data:
            target_data = reads_by_stage_data[sample_name]

            for target_name, stages_data in target_data.items():
                # Build stages list
                stages = []
                for stage_name, stage_info in stages_data.items():
                    # stage_info is now a dictionary with stage, read_count, and additional columns
                    stages.append(stage_info)

                # Create target entry
                target_entry = {
                    "target_name": target_name,
                    "stages": stages,
                }
                read_counts_for_targets.append(target_entry)

        # Add read_counts_for_targets if we have data
        if read_counts_for_targets:
            library_sample_entry["read_counts_for_targets"] = read_counts_for_targets

        read_counts_by_library_sample_by_stage.append(library_sample_entry)

    # Build the final output structure
    output_data = {
        "read_counts_by_library_sample_by_stage": read_counts_by_library_sample_by_stage,
    }
    if bioinformatics_run_name is not None:
        output_data["bioinformatics_run_name"] = bioinformatics_run_name
    return output_data
