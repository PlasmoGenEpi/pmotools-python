#!/usr/bin/env python3
import pandas as pd
import numpy as np
from ..pmo_builder.json_convert_utils import check_additional_columns_exist
import warnings


def panel_info_table_to_pmo_dict(target_table: pd.DataFrame,
                                 panel_name: str,
                                 genome_info: dict,
                                 target_name_col: str = 'target_name',
                                 forward_primers_seq_col: str = 'fwd_primer',
                                 reverse_primers_seq_col: str = 'rev_primer',
                                 reaction_name_col: str | None = None,
                                 forward_primers_start_col: int | None = None,
                                 forward_primers_end_col: int | None = None,
                                 reverse_primers_start_col: int | None = None,
                                 reverse_primers_end_col: int | None = None,
                                 insert_start_col: int | None = None,
                                 insert_end_col: int | None = None,
                                 chrom_col: str | None = None,
                                 strand_col: str | None = None,
                                 ref_seq_col: str | None = None,
                                 gene_name_col: str | None = None,
                                 target_attributes_col: str | None = None,
                                 additional_target_info_cols: list | None = None,
                                 ):
    """
    Convert a dataframe containing panel information into dictionary of targets and reference information


    :param target_table: The dataframe containing the target information
    :param panel_name: the panel ID assigned to the panel
    :param genome_info: A dictionary containing the genome information
    :param target_name_col: the name of the column containing the target IDs
    :param forward_primers_seq_col: the name of the column containing the sequence of the forward primer
    :param reverse_primers_seq_col: the name of the column containing the sequence of the reverse primer
    :param reaction_name_col: the name of the column containing which reaction the target was part of. By default they will all be put in one reaction.
    :param forward_primers_start_col (Optional): the name of the column containing the 0-based start coordinate of the forward primer
    :param forward_primers_end_col (Optional): the name of the column containing the 0-based end coordinate of the forward primer
    :param reverse_primers_start_col (Optional): the name of the column containing the 0-based start coordinate of the reverse primer
    :param reverse_primers_end_col (Optional): the name of the column containing the 0-based end coordinate of the reverse primer
    :param insert_start_col (Optional): the name of the column containing the 0-based start coordinate of the insert
    :param insert_end_col (Optional): the name of the column containing the 0-based end coordinate of the insert
    :param chrom_col (Optional): the name of the column containing the chromosome for the target
    :param gene_name_col (Optional): the name of the column containing the gene id
    :param strand_col (Optional): the name of the column containing the strand for the target
    :param target_attributes_col (Optional): a list of classification type for the primer target
    :param additional_target_info_cols (Optional): dictionary of optional additional columns to add to the target information dictionary. Keys are column names and values are the type.
    :return: a dict of the panel information
    """

    if not isinstance(target_table, pd.DataFrame):
        raise ValueError("target_table must be a pandas DataFrame.")
    # check_genome_info(genome_info)

    # Check additional columns if any are added
    check_additional_columns_exist(target_table, additional_target_info_cols)

    # Check start, end and chrom includied for each if location is added to any of forward primer, reverse primer or insert
    location_cols = check_location_columns(forward_primers_start_col, forward_primers_end_col, reverse_primers_start_col,
                                           reverse_primers_end_col, insert_start_col, insert_end_col, chrom_col, strand_col, ref_seq_col)

    # Create dictionary of targets
    targets_dict = create_targets_dict(target_table, target_name_col, forward_primers_seq_col, reverse_primers_seq_col,
                                       location_info_cols=location_cols, gene_name_col=gene_name_col, target_attributes_col=target_attributes_col,
                                       additional_target_info_cols=additional_target_info_cols)
    # Put together components
    panel_info_dict = {"panel_info": [], "targeted_genomes": [
        genome_info], "target_info": targets_dict}
    return panel_info_dict


def build_panel_info(panel_name, target_table, reaction_name_col):
    panel_dict = {'panel_name': panel_name, "reactions": []}
    if reaction_name_col:
        for reaction in target_table[reaction_name_col].unique():
            reaction_target_table = target_table[target_table[reaction_name_col] == reaction]
            reaction_dict = {'reaction_name': reaction, 'panel_targets': []}
            # TODO: USE get_index_of_specimen_names from recent PR


def check_genome_info(genome_info):
    if isinstance(genome_info, dict):
        required_keys = {"name", "genome_version", "taxon_id", "url"}
        missing_keys = required_keys - genome_info.keys()
        if missing_keys:
            raise ValueError(
                f"genome_info missing required keys: {', '.join(missing_keys)}")
    else:
        raise TypeError(
            f"genome_info must be a dict, but got {type(genome_info).__name__}")


def create_targets_dict(
    target_table: pd.DataFrame,
    target_name_col: str,
    forward_primers_seq_col: str,
    reverse_primers_seq_col: str,
    genome_id: int = 0,
    location_info_cols: list | None = None,
    gene_name_col: str = None,
    target_attributes_col: str = None,
    additional_target_info_cols: list = None
):
    """
    Convert the read in target information into the panel information dictionary


    :param target_table: The dataframe containing the target information
    :param target_name_col: the name of the column containing the target IDs
    :param forward_primers_seq_col: the name of the column containing the sequence of the forward primer
    :param reverse_primers_seq_col: the name of the column containing the sequence of the reverse primer
    :param location_info_cols: list of column names for location information
    :param gene_name_col: the name of the column containing the gene id
    :param target_attributes_col: a list of classification type for the primer target
    :return: a dictionary of the target information
    """
    # Check targets before putting into JSON
    (
        forward_primers_start_col,
        forward_primers_end_col,
        reverse_primers_start_col,
        reverse_primers_end_col,
        insert_start_col,
        insert_end_col,
        chrom_col,
        strand_col,
        ref_seq_col,
    ) = location_info_cols if location_info_cols else [None] * 9

    # Check target information in the dataframe
    check_targets_are_unique(target_table, target_name_col)
    columns_to_check = [forward_primers_seq_col, reverse_primers_seq_col]
    if location_info_cols:
        columns_to_check += [col for col in location_info_cols if col]
    check_unique_target_info(target_table, target_name_col, columns_to_check)
    missing_insert_loc, missing_fwd_primer_loc, missing_rev_primer_loc = summarise_targets_missing_optional_info(target_table, target_name_col, chrom_col, insert_start_col, insert_end_col,
                                                                                                                 forward_primers_start_col, forward_primers_end_col, reverse_primers_start_col, reverse_primers_end_col)

    # Put targets together in dictionary
    targets_dicts = []
    for _, row in target_table.iterrows():
        target_name = row[target_name_col]
        target_dict = {
            "target_name": target_name,
        }
        if gene_name_col:
            target_dict["gene_name"] = row[gene_name_col]
        if target_attributes_col:
            target_dict["target_attributes_col"] = row[target_attributes_col]
        if additional_target_info_cols:
            for col in additional_target_info_cols:
                value = row[col]
                # Convert numpy types to native Python types
                if isinstance(value, (np.integer, np.int64)):
                    value = int(value)
                elif isinstance(value, (np.floating, np.float64)):
                    value = float(value)
                elif pd.isna(value):
                    value = None
                target_dict[col] = value

        # Add insert location info if location_info_cols are provided
        if insert_start_col and target_name not in missing_insert_loc:
            target_dict["insert_location"] = {
                "genome_id": genome_id,
                "chrom": row[chrom_col],
                "start": int(row[insert_start_col]),
                "end": int(row[insert_end_col]),
            }
            if strand_col and pd.notna(row[strand_col]):
                target_dict["insert_location"]["strand"] = row[strand_col]
            if ref_seq_col and pd.notna(row[ref_seq_col]):
                target_dict["insert_location"]["ref_seq"] = row[ref_seq_col]

        # Extract primer information for each row
        fwd_primer_dict = {"seq": row[forward_primers_seq_col]}
        rev_primer_dict = {"seq": row[reverse_primers_seq_col]}
        if forward_primers_start_col and target_name not in missing_fwd_primer_loc:
            fwd_primer_dict["location"] = {
                "genome_id": genome_id,
                "chrom": row[chrom_col],
                "end": int(row[forward_primers_start_col]),
                "start": int(row[forward_primers_end_col]),
            }
            if strand_col and pd.notna(row[strand_col]):
                fwd_primer_dict["location"]["strand"] = row[strand_col]
        if reverse_primers_start_col and target_name not in missing_rev_primer_loc:
            rev_primer_dict["location"] = {
                "genome_id": genome_id,
                "chrom": row[chrom_col],
                "end": int(row[reverse_primers_end_col]),
                "start": int(row[reverse_primers_start_col]),
            }
            if strand_col and pd.notna(row[strand_col]):
                rev_primer_dict["location"]["strand"] = row[strand_col]
        target_dict["forward_primer"] = fwd_primer_dict
        target_dict["reverse_primer"] = rev_primer_dict

        targets_dicts.append(target_dict)

    return targets_dicts


def check_location_columns(
    forward_primers_start_col, forward_primers_end_col, reverse_primers_start_col,
    reverse_primers_end_col, insert_start_col, insert_end_col, chrom_col, strand_col, ref_seq_col
):
    location_cols = [
        forward_primers_start_col,
        forward_primers_end_col,
        reverse_primers_start_col,
        reverse_primers_end_col,
        insert_start_col,
        insert_end_col,
        chrom_col,
        strand_col,
        ref_seq_col
    ]
    if any(location_cols):
        warnings = []
        if not chrom_col:
            warnings.append(
                "If including location information (any of forward_primers_start_col, forward_primers_end_col, reverse_primers_start_col, reverse_primers_end_col, insert_start_col, insert_end_col) chrom_col must be set."
            )
        if (forward_primers_start_col is None) != (forward_primers_end_col is None):
            warnings.append(
                "If one of forward_primers_start_col or forward_primers_end_col is set, then both must be."
            )
        if (reverse_primers_start_col is None) != (reverse_primers_end_col is None):
            warnings.append(
                "If one of reverse_primers_start_col or reverse_primers_end_col is set, then both must be."
            )
        if (insert_start_col is None) != (insert_end_col is None):
            warnings.append(
                "If one of insert_start_col or insert_end_col is set, then both must be."
            )
        if warnings:
            raise ValueError(
                "Errors with location column configuration:\n- " + "\n- ".join(warnings))
        return location_cols
    return None


def check_targets_are_unique(df, target_name_col):
    duplications = df[df[target_name_col].duplicated(keep=False)]
    if not duplications.empty:
        raise ValueError(
            f"The following target_ids are duplicated: {duplications[target_name_col].unique()}")


def check_unique_target_info(df, target_name_col, columns_to_check):
    groups = (
        df.groupby(columns_to_check)[target_name_col]
        .apply(list)
        .reset_index(name=target_name_col)
    )

    # Keep only groups where more than one target shares the same primer pair
    duplicated_groups = groups[groups[target_name_col].str.len() > 1]

    if not duplicated_groups.empty:
        msg_lines = ["The following targets have duplicated information:"]
        for _, row in duplicated_groups.iterrows():
            cols_info = ", ".join(
                f"{col}={row[col]}" for col in columns_to_check)
            targets = ", ".join(map(str, row[target_name_col]))
            msg_lines.append(f"targets: {targets} → {cols_info}")

        raise ValueError("\n".join(msg_lines))


def summarise_targets_missing_optional_info(df, target_name_col, chrom_col, insert_start_col, insert_end_col, forward_primers_start_col, forward_primers_end_col, reverse_primers_start_col, reverse_primers_end_col):
    missing_insert_loc = None
    missing_fwd_primer_loc = None
    missing_rev_primer_loc = None

    def check_missing(name, cols):
        missing = df[df[cols].isnull().any(axis=1)][target_name_col].tolist()
        if len(missing) > 0:
            warnings.warn(
                f"{name} location information was not added for the following targets that had empty fields: {', '.join(missing)}")
        return missing

    missing_insert_loc = check_missing(
        'Insert', [chrom_col, insert_start_col, insert_end_col]) if insert_start_col else None
    missing_fwd_primer_loc = check_missing('Forward primer', [
        chrom_col, forward_primers_start_col, forward_primers_end_col]) if forward_primers_start_col else None
    missing_rev_primer_loc = check_missing('Reverse primer', [
        chrom_col, reverse_primers_start_col, reverse_primers_end_col]) if reverse_primers_start_col else None

    return missing_insert_loc, missing_fwd_primer_loc, missing_rev_primer_loc
