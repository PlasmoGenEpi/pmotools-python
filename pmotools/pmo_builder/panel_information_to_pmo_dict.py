#!/usr/bin/env python3
import pandas as pd
import json
import numpy as np
from .json_convert_utils import check_additional_columns_exist


def panel_info_table_to_pmo_dict(contents: pd.DataFrame,
                                 panel_name: str,
                                 genome_info: dict,
                                 target_name_col: str = 'target_name',
                                 reaction_name_col: str = 'reaction_name',
                                 forward_primers_seq_col: str = 'fwd_primer',
                                 reverse_primers_seq_col: str = 'rev_primer',
                                 forward_primers_start_col: str | None = None,
                                 forward_primers_end_col: str | None = None,
                                 reverse_primers_start_col: str | None = None,
                                 reverse_primers_end_col: str | None = None,
                                 insert_start_col: str | None = None,
                                 insert_end_col: str | None = None,
                                 insert_ref_seq_col: str | None = None,
                                 chrom_col: str | None = None,
                                 strand_col: str | None = None,
                                 gene_name_col: str | None = None,
                                 target_attributes_col: str | None = None,
                                 additional_target_info_cols: list | None = None,
                                 ):
    """
    Convert a dataframe containing panel information into dictionary of targets and reference information

    :param contents (pd.DataFrame): The dataframe containing the target information
    :param panel_name (str): The panel ID assigned to the panel
    :param genome_info (str): A dictionary containing the genome information. Must include name, genome_version, taxon_id, and url. Can also include chromosomes, gff_url, and additional fields.
    :param target_name_col (str): the name of the column containing the target IDs
    :param reaction_name_col (str): the name of column containing the reaction that the target was in
    :param forward_primers_seq_col (str): the name of the column containing the sequence of the forward primer
    :param reverse_primers_seq_col (str): the name of the column containing the sequence of the reverse primer
    :param forward_primers_start_col (Optional[str]): the name of the column containing the 0-based start coordinate of the forward primer
    :param forward_primers_end_col (Optional[str]): the name of the column containing the 0-based end coordinate of the forward primer
    :param reverse_primers_start_col (Optional[str]): the name of the column containing the 0-based start coordinate of the reverse primer
    :param reverse_primers_end_col (Optional[str]): the name of the column containing the 0-based end coordinate of the reverse primer
    :param insert_start_col (Optional[str]): the name of the column containing the 0-based start coordinate of the insert
    :param insert_end_col (Optional[str]): the name of the column containing the 0-based end coordinate of the insert
    :param chrom_col (Optional[str]): the name of the column containing the chromosome for the target
    :param strand_col (Optional[str]): the name of the column containing the strand for the target
    :param gene_name_col (Optional[str]): the name of the column containing the gene name being targeted. 
    :param target_attributes_col (Optional[str]): A classification type for the target, or other tags for target
    :param additional_target_info_cols (Optional[str]): dictionary of optional additional columns to add to the target information dictionary. Keys are column names and values are the type.

    :return: a dict of the panel information
    """

    if not isinstance(contents, pd.DataFrame):
        raise ValueError("contents must be a pandas DataFrame.")
    if not isinstance(genome_info, dict):
        raise ValueError("genome_info must be a dictionary.")

    # Check additional columns if any are added
    check_additional_columns_exist(contents, additional_target_info_cols)

    # If one location column set, check all location columns are set
    location_cols = check_and_create_location_columns(forward_primers_start_col, forward_primers_end_col, reverse_primers_start_col,
                                                      reverse_primers_end_col, insert_start_col, insert_end_col, chrom_col, strand_col, insert_ref_seq_col)

    # Create dictionary of targets
    targets_dict = create_target_info_dict(contents, target_name_col, forward_primers_seq_col, reverse_primers_seq_col,
                                           location_info_cols=location_cols, gene_name_col=gene_name_col, target_attributes_col=target_attributes_col,
                                           additional_target_info_cols=additional_target_info_cols)
    # Put together components
    panel_info_dict = {"panel_info": {panel_name: {"panel_name": panel_name,
                       "target_genome": genome_info, "targets": targets_dict}}}

    return panel_info_dict


def create_target_info_dict(
    contents: pd.DataFrame,
    target_name_col: str,
    forward_primers_seq_col: str,
    reverse_primers_seq_col: str,
    location_info_cols: list | None = None,
    gene_name_col: str = None,
    target_attributes_col: str = None,
    additional_target_info_cols: list = None
):
    """
    Convert the read in target information into the panel information dictionary


    :param contents: The dataframe containing the target information
    :param target_name_col: the name of the column containing the target IDs
    :param forward_primers_seq_col: the name of the column containing the sequence of the forward primer
    :param reverse_primers_seq_col: the name of the column containing the sequence of the reverse primer
    :param location_info_cols: list of column names for location information
    :param gene_name_col: the name of the column containing the gene name
    :param chrom_col: the name of the column containing the chromosome the set of primers target
    :param strand_col: the name of the column containing the strand for the target
    :return: a dictionary of the target information
    """
    # Check targets before putting into JSON
    columns_to_check = []
    if location_info_cols:
        forward_primers_start_col, forward_primers_end_col, reverse_primers_start_col, reverse_primers_end_col, insert_start_col, insert_end_col, chrom_col, strand_col = location_info_cols
        columns_to_check = columns_to_check + \
            [chrom_col, strand_col, insert_start_col, insert_end_col]
    if gene_name_col:
        columns_to_check.append(gene_name_col)
    if target_attributes_col:
        columns_to_check.append(target_attributes_col)
    if additional_target_info_cols:
        columns_to_check = columns_to_check + additional_target_info_cols

    # Check these columns are unique for each target
    check_columns_unique_for_target(
        contents, target_name_col, columns_to_check)

    # Put targets together in dictionary
    targets_dict = {}
    for target_name in contents[target_name_col].unique():
        target_df = contents.loc[contents[target_name_col] == target_name]
        target_dict = {
            "target_name": target_name,
            "forward_primers": [],
            "reverse_primers": [],
        }
        if gene_name_col:
            target_dict["gene_name"] = target_df[gene_name_col].iloc[0]
        if target_attributes_col:
            target_dict["target_attributes"] = target_df[target_attributes_col].iloc[0]
        if additional_target_info_cols:
            for col in additional_target_info_cols:
                value = target_df[col].iloc[0]
                # Convert numpy types to native Python types
                if isinstance(value, (np.integer, np.int64)):
                    value = int(value)
                elif isinstance(value, (np.floating, np.float64)):
                    value = float(value)
                elif pd.isna(value):
                    value = None
                target_dict[col] = value
        # Add insert location info if location_info_cols are provided
        if location_info_cols:
            target_dict["insert_location"] = {
                "chrom": target_df[chrom_col].iloc[0],
                "start": int(target_df[insert_start_col].iloc[0]),
                "end": int(target_df[insert_end_col].iloc[0]),
                "strand": target_df[strand_col].iloc[0]
            }
        # Extract primer information for each row
        for _, row in target_df.iterrows():
            fwd_primer_dict = {"seq": row[forward_primers_seq_col]}
            rev_primer_dict = {"seq": row[reverse_primers_seq_col]}
            if location_info_cols:
                fwd_primer_dict["location"] = {
                    "chrom": row[chrom_col],
                    "end": int(row[forward_primers_end_col]),
                    "start": int(row[forward_primers_start_col]),
                    "strand": row[strand_col]
                }
                rev_primer_dict["location"] = {
                    "chrom": row[chrom_col],
                    "end": int(row[reverse_primers_end_col]),
                    "start": int(row[reverse_primers_start_col]),
                    "strand": row[strand_col]
                }
            target_dict["forward_primers"].append(fwd_primer_dict)
            target_dict["reverse_primers"].append(rev_primer_dict)

        targets_dict[target_name] = target_dict
    return targets_dict


def check_and_create_location_columns(
    forward_primers_start_col,
    forward_primers_end_col,
    reverse_primers_start_col,
    reverse_primers_end_col,
    insert_start_col,
    insert_end_col,
    chrom_col,
    strand_col,
    insert_ref_seq_col
):
    location_cols = [
        forward_primers_start_col,
        forward_primers_end_col,
        reverse_primers_start_col,
        reverse_primers_end_col,
        insert_start_col,
        insert_end_col,
        insert_ref_seq_col,
        chrom_col,
        strand_col,
    ]
    if any(location_cols):
        if any(location_cols[0:3, 8]):
            if not all(location_cols):
                raise ValueError(
                    "If one of the location params (forward_primers_start_col, forward_primers_end_col, reverse_primers_start_col, reverse_primers_end_col, chrom_col) is set then all must be set."
                )
        if (insert_start_col == None) != (insert_end_col == None):
            raise ValueError(
                'If either insert_start_col or insert_end_col set, then both need to be.'
            )
        if insert_ref_seq_col and not insert_start_col:
            raise ValueError(
                'If insert_ref_seq_col is set then insert_start_col')
        return location_cols
    return None


def check_columns_unique_for_target(df, target_name_col, columns_to_check):
    for col in columns_to_check:
        duplicates = df.groupby(target_name_col)[col].nunique()
        duplicates = duplicates[duplicates > 1]
        if not duplicates.empty:
            duplicate_targets = duplicates.index.tolist()
            raise ValueError(
                f"The following target_names have multiple unique {col}: {duplicate_targets}")
