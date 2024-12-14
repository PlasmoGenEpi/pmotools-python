#!/usr/bin/env python3
import pandas as pd
import json
from .json_convert_utils import check_additional_columns_exist


def demultiplexed_targets_to_pmo_dict(
    contents: pd.DataFrame,
    bioinfo_id: str,
    sampleID_col: str = 'sampleID',
    target_id_col: str = 'target_id',
    read_count_col: str = 'reads',
    additional_hap_detected_cols: dict | None = None
):
    """
    Convert a dataframe of microhaplotype calls into a dictionary for detected haplotypes 
    and representative haplotype sequences.

    :param contents: DataFrame containing demultiplexed sample information
    :param bioinfo_id: Bioinformatics ID of the demultiplexed targets
    :param sampleID_col: Name of the column containing sample IDs
    :param target_id_col: Name of the column containing locus IDs
    :param read_count_col: Name of the column containing read counts
    :param additional_hap_detected_cols: Optional columns to include in the output,
                                         with keys as column names and values as their output names
    :return: JSON string containing the processed data
    """
    # Validate additional columns
    check_additional_columns_exist(contents, additional_hap_detected_cols)

    # Initialize the main dictionary
    demultiplexed_targets_dict = {}

    # Group by sampleID for efficiency
    grouped = contents.groupby(sampleID_col)
    for sample, group in grouped:
        sample_dict = {"demultiplexed_targets": {
            "experiment_sample_id": sample}}

        # Process each target within the sample group
        for _, row in group.iterrows():
            target = row[target_id_col]
            reads = int(row[read_count_col])  # Ensure integer

            target_dict = {
                "raw_read_count": reads,
                "target_id": target
            }

            # Add optional additional columns
            if additional_hap_detected_cols:
                target_dict.update({
                    output_name: row[col]
                    for col, output_name in additional_hap_detected_cols.items()
                })

            # Add target data to the sample dictionary
            sample_dict["demultiplexed_targets"][target] = target_dict

        # Update the main dictionary with the sample data
        demultiplexed_targets_dict[sample] = sample_dict

    # Final data structure
    output_data = {
        "target_demultiplexed_experiment_samples": {
            bioinfo_id: {
                "demultiplexed_experiment_samples": demultiplexed_targets_dict,
                "tar_amp_bioinformatics_info_id": bioinfo_id
            }
        }
    }

    # Serialize and return the result
    return output_data
