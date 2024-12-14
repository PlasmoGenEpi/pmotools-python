#!/usr/bin/env python3
import pandas as pd
import json
from .json_convert_utils import check_additional_columns_exist


def microhaplotype_table_to_pmo_dict(
    contents: pd.DataFrame,
    bioinfo_id: str,
    sampleID_col: str = 'sampleID',
    locus_col: str = 'locus',
    mhap_col: str = 'asv',
    reads_col: str = 'reads',
    additional_hap_detected_cols: list | None = None
):
    """
    Convert a dataframe of a microhaplotype calls into a dictionary containing a dictionary for the haplotypes_detected and a dictionary for the representative_haplotype_sequences.

    :param contents: The dataframe containing microhaplotype calls
    :param bioinfo_id: the bioinformatics ID of the microhaplotype table
    :param sampleID_col: the name of the column containing the sample IDs
    :param locus_col: the name of the column containing the locus IDs
    :param mhap_col: the name of the column containing the microhaplotype sequence
    :param reads_col: the name of the column containing the reads counts
    :param additional_hap_detected_cols: optional additional columns to add to the microhaplotype detected dictionary, the key is the pandas column and the value is what to name it in the output
    :return: a dict of both the haplotypes_detected and representative_haplotype_sequences
    """

    representative_microhaplotype_dict = create_representative_microhaplotype_dict(
        contents, locus_col, mhap_col)

    detected_mhap_dict = create_detected_microhaplotype_dict(contents, sampleID_col, locus_col,
                                                             mhap_col, reads_col, representative_microhaplotype_dict,
                                                             additional_hap_detected_cols)

    output_data = {"microhaplotypes_detected": {bioinfo_id: {'experiment_samples': detected_mhap_dict}},
                   "representative_microhaplotype_sequences": {bioinfo_id: {"representative_microhaplotype_id": bioinfo_id, 'targets': representative_microhaplotype_dict}}
                   }
    return output_data


def create_representative_microhaplotype_dict(
        microhaplotype_table: pd.DataFrame,
        locus_col: str,
        mhap_col: str
):
    """
    Convert the read-in microhaplotype calls table into a representative microhaplotype JSON-like dictionary.

    :param microhaplotype_table: The parsed microhaplotype calls table.
    :param locus_col: The name of the column containing the locus IDs.
    :param mhap_col: The name of the column containing the microhaplotype sequence.
    :return: A dictionary formatted for JSON output with representative microhaplotype sequences.
    """
    # Drop duplicates and reset index
    unique_table = microhaplotype_table[[
        locus_col, mhap_col]].drop_duplicates().reset_index(drop=True)

    json_data = {
        locus: {
            "seqs": {
                f"{locus}.{idx}": {
                    "microhaplotype_id": f"{locus}.{idx}",
                    "seq": seq
                }
                for idx, seq in enumerate(group[mhap_col])
            }
        }
        for locus, group in unique_table.groupby(locus_col)
    }

    return json_data


def create_detected_microhaplotype_dict(
    microhaplotype_table: pd.DataFrame,
    sampleID_col: str,
    locus_col: str,
    mhap_col: str,
    reads_col: str,
    representative_microhaplotype_dict: dict,
    additional_hap_detected_cols: list | None = None
):
    """
    Convert the read-in microhaplotype calls table into the detected microhaplotype dictionary.

    :param microhaplotype_table: Parsed microhaplotype calls table.
    :param sampleID_col: Column containing the sample IDs.
    :param locus_col: Column containing the locus IDs.
    :param mhap_col: Column containing the microhaplotype sequences.
    :param reads_col: Column containing the read counts.
    :param representative_microhaplotype_dict: Dictionary of representative microhaplotypes.
    :param additional_hap_detected_cols: Optional additional columns to add to the microhaplotypes detected, the key is the pandas column and the value is what to name it in the output.
    :return: A dictionary of detected microhaplotype results.
    """
    # Validate additional columns if provided
    if additional_hap_detected_cols:
        check_additional_columns_exist(
            microhaplotype_table, additional_hap_detected_cols)

    # Map sequences to representative haplotype IDs for fast lookup
    rep_hap_map = {
        (locus, seq["seq"]): seq["microhaplotype_id"]
        for locus, reps in representative_microhaplotype_dict.items()
        for seq in reps["seqs"].values()
    }

    def build_haplotype_info(row):
        """Helper to construct haplotype info for each row."""
        locus, seq = row[locus_col], row[mhap_col]
        matching_id = rep_hap_map.get((locus, seq))
        if not matching_id:
            raise ValueError(
                f"No representative haplotype ID found for {seq} at locus {locus}")

        haplotype_info = {
            "haplotype_id": matching_id,
            "read_count": row[reads_col],
        }

        if additional_hap_detected_cols:
            haplotype_info.update({input_col: row[input_col]
                                   for input_col in additional_hap_detected_cols})

        return matching_id, haplotype_info

    # Build the JSON-like structure
    json_data = (
        microhaplotype_table
        .groupby(sampleID_col)
        .apply(lambda sample_group: {
            "sample_id": sample_group[sampleID_col].iloc[0],
            "target_results": {
                locus: {
                    "microhaplotypes": {
                        hap_id: hap_info
                        for _, row in locus_group.iterrows()
                        for hap_id, hap_info in [build_haplotype_info(row)]
                    }
                }
                for locus, locus_group in sample_group.groupby(locus_col)
            }
        })
        .to_dict()
    )
    return json_data
