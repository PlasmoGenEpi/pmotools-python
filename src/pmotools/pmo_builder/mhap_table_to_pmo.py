#!/usr/bin/env python3
import pandas as pd
import warnings

from ..pmo_builder.json_convert_utils import check_additional_columns_exist


def mhap_table_to_pmo(
    microhaplotype_table: pd.DataFrame,
    bioinformatics_run_name: str,
    library_sample_name_col: str = "library_sample_name",
    target_name_col: str = "target_name",
    seq_col: str = "seq",
    reads_col: str = "reads",
    genome_id=0,
    umis_col: str | None = None,
    chrom_col: str | None = None,
    start_col: str | None = None,
    end_col: str | None = None,
    ref_seq_col: str | None = None,
    strand_col: str | None = None,
    alt_annotations_col: str | None = None,
    masking_seq_start_col: str | None = None,
    masking_seq_segment_size_col: str | None = None,
    masking_replacement_size_col: str | None = None,
    masking_delim: str = ",",
    microhaplotype_name_col: str | None = None,
    pseudocigar_col: str | None = None,
    quality_col: str | None = None,
    additional_representative_mhap_cols: str | None = None,
    additional_mhap_detected_cols: list | None = None,
):
    """
    Convert a dataframe of a microhaplotype calls into a dictionary containing a dictionary for the haplotypes_detected and a dictionary for the representative_haplotype_sequences.

    :param microhaplotype_table (pd.DataFrame): The dataframe containing microhaplotype calls
    :param bioinformatics_run_name (str) : Unique name for the bioinformatics run that generated the data (column name or individual run name).
    :param library_sample_name_col (str) : the name of the column containing the experiment sample names. Default: library_sample_name
    :param target_name_col (str) : the name of the column containing the targets. Default: target_name
    :param seq_col (str) : the name of the column containing the microhaplotype sequences. Default: seq
    :param reads_col (str) : the name of the column containing the reads counts. Default: reads
    :param genome_id (int) : the ID of the genome used as reference. Default: 0
    :param umis_col (Optional[str]) : the name of the column with unique molecular identifier count associated with this microhaplotype
    :param chrom_col (Optional[str]) : the name of the column containing the chromosome name of the microhaplotype
    :param start_col (Optional[str]) : the name of the column containing the start of the microhaplotype
    :param end_col (Optional[str]) : the name of the column containing the end of the microhaplotype
    :param ref_seq_col (Optional[str]) : the name of the column containing the reference sequence for the microhaplotype
    :param strand_col (Optional[str]) : the name containing the strand of the microhaplotype
    :param alt_annotations_col (Optional[str]) : the name of the column containing any alternative annotations
    :param masking_seq_start_col (Optional[str]) : the name ofthe column containing a list of start positions for masking
    :param masking_seq_segment_size_col (Optional[str]) : the name of the column containing a list of lengths of the segments in seq being masked
    :param masking_replacement_size_col (Optional[str]) : the name of the column containing a list of lengths of the masking replacements
    :param masking_delim (Optional[str]) : delim of the masking information. Default: ','
    :param microhaplotype_name_col (Optional[str]) : the name of the column containing an optional name for this microhaplotype
    :param pseudocigar_col (Optional[str]) : the name of the column containing a pseudocigar for the microhaplotype
    :param quality_col (Optional[str]) : the name of the column containing the ansi fastq per base quality score for this sequence
    :param additional_representative_mhap_cols (Optional[List[str], None]]): additional columns to add to the representative microhaplotypes table.
    :param additional_mhap_detected_cols (Optional[List[str], None]]): additional columns to add to the detected microhaplotypes table.

    :return: a dict of both the haplotypes_detected and representative_haplotype_sequences
    """

    representative_microhaplotype_dict = create_representative_microhaplotype_dict(
        microhaplotype_table,
        target_name_col,
        seq_col,
        genome_id,
        chrom_col,
        start_col,
        end_col,
        ref_seq_col,
        strand_col,
        alt_annotations_col,
        masking_seq_start_col,
        masking_seq_segment_size_col,
        masking_replacement_size_col,
        masking_delim,
        microhaplotype_name_col,
        pseudocigar_col,
        quality_col,
        additional_representative_mhap_cols,
    )

    detected_mhap_dict_list = []
    if bioinformatics_run_name in microhaplotype_table.columns:
        for bioinfo_run in microhaplotype_table[bioinformatics_run_name].unique():
            microhaplotype_table_per_run = microhaplotype_table[
                microhaplotype_table[bioinformatics_run_name] == bioinfo_run
            ]
            detected_mhap_dict = create_detected_microhaplotype_dict(
                microhaplotype_table_per_run,
                bioinfo_run,
                representative_microhaplotype_dict,
                library_sample_name_col,
                target_name_col,
                seq_col,
                reads_col,
                umis_col,
                additional_mhap_detected_cols,
            )
            detected_mhap_dict_list.append(detected_mhap_dict)
    else:
        detected_mhap_dict = create_detected_microhaplotype_dict(
            microhaplotype_table,
            bioinformatics_run_name,
            representative_microhaplotype_dict,
            library_sample_name_col,
            target_name_col,
            seq_col,
            reads_col,
            umis_col,
            additional_mhap_detected_cols,
        )
        detected_mhap_dict_list.append(detected_mhap_dict)

    output_data_dict = {
        "representative_microhaplotypes": representative_microhaplotype_dict,
        "detected_microhaplotypes": detected_mhap_dict_list,
    }
    return output_data_dict


def create_representative_microhaplotype_dict(
    microhaplotype_table: pd.DataFrame,
    target_name_col: str = "target_name",
    seq_col: str = "seq",
    genome_id: int = 0,
    chrom_col: str | None = None,
    start_col: str | None = None,
    end_col: str | None = None,
    ref_seq_col: str | None = None,
    strand_col: str | None = None,
    alt_annotations_col: str | None = None,
    masking_seq_start_col: str | None = None,
    masking_seq_segment_size_col: str | None = None,
    masking_replacement_size_col: str | None = None,
    masking_delim: str = ",",
    microhaplotype_name_col: str | None = None,
    pseudocigar_col: str | None = None,
    quality_col: str | None = None,
    additional_representative_mhap_cols: list[str] | None = None,
):
    """
    Convert the read-in microhaplotype calls table into a representative microhaplotype JSON-like dictionary.

    :param microhaplotype_table (pd.DataFrame): The dataframe containing microhaplotype calls
    :param target_name_col (str) : the name of the column containing the targets. Default: target_name
    :param seq_col (str) : the name of the column containing the microhaplotype sequences. Default: seq
    :param genome_id (int) : the genome ID
    :param chrom_col (Optional[str]) : the name of the column containing the chromosome name of the microhaplotype
    :param start_col (Optional[str]) : the name of the column containing the start of the microhaplotype
    :param end_col (Optional[str]) : the name of the column containing the end of the microhaplotype
    :param ref_seq_col (Optional[str]) : the name of the column containing the reference sequence for the microhaplotype
    :param strand_col (Optional[str]) : the name containing the strand of the microhaplotype
    :param alt_annotations_col (Optional[str]) : the name of the column containing any alternative annotations
    :param masking_seq_start_col (Optional[str]) : the name ofthe column containing a list of start positions for masking
    :param masking_seq_segment_size_col (Optional[str]) : the name of the column containing a list of lengths of the segments in seq being masked
    :param masking_replacement_size_col (Optional[str]) : the name of the column containing a list of lengths of the masking replacements
    :param masking_delim (Optional[str]) : delim of the masking information. Default: ','
    :param microhaplotype_name_col (Optional[str]) : the name of the column containing an optional name for this microhaplotype
    :param pseudocigar_col (Optional[str]) : the name of the column containing a pseudocigar for the microhaplotype
    :param quality_col (Optional[str]) : the name of the column containing the ansi fastq per base quality score for this sequence
    :param additional_representative_mhap_cols (Optional[List[str], None]]): additional columns to add to the representative microhaplotypes table.

    :return: A dictionary formatted for JSON output with representative microhaplotype sequences.
    """

    if additional_representative_mhap_cols:
        check_additional_columns_exist(
            microhaplotype_table, additional_representative_mhap_cols
        )

    def get_if_present(row, col):
        return row[col] if col and pd.notna(row[col]) else None

    def extract_masking(row):
        if not (
            masking_seq_start_col
            and masking_seq_segment_size_col
            and masking_replacement_size_col
        ):
            return []
        if all(
            [
                pd.notna(row[masking_seq_start_col]),
                pd.notna(row[masking_seq_segment_size_col]),
                pd.notna(row[masking_replacement_size_col]),
            ]
        ):
            starts = str(row[masking_seq_start_col]).split(masking_delim)
            sizes = str(row[masking_seq_segment_size_col]).split(masking_delim)
            replacements = str(row[masking_replacement_size_col]).split(masking_delim)
            return [
                {
                    "seq_start": int(s),
                    "seq_segment_size": int(sz),
                    "replacement_size": int(r),
                }
                for s, sz, r in zip(starts, sizes, replacements)
                if s and sz and r
            ]
        else:
            return []

    def warn_if_duplicated_seqs(df, target_col, seq_col):
        dup_counts = df.groupby([target_col, seq_col]).size()
        duplicate_combos = dup_counts[dup_counts > 1]

        if not duplicate_combos.empty:
            warnings.warn(
                f"Duplicate (target, asv) combinations found:\n{duplicate_combos}",
                UserWarning,
            )

    # Determine which columns to keep
    optional_cols = [
        chrom_col,
        start_col,
        end_col,
        ref_seq_col,
        strand_col,
        alt_annotations_col,
        microhaplotype_name_col,
        pseudocigar_col,
        quality_col,
    ]
    masking_cols = [
        masking_seq_start_col,
        masking_seq_segment_size_col,
        masking_replacement_size_col,
    ]

    # Check location cols are set correctly
    if any(masking_cols):
        if not all(masking_cols):
            raise ValueError(
                "If one of masking_seq_start_col, masking_seq_segment_size_col, masking_replacement_size_col is set, then all must be."
            )
    if any([chrom_col, start_col, end_col, ref_seq_col, strand_col]):
        if not all([chrom_col, start_col, end_col]):
            raise ValueError(
                "If any location columns set (chrom_col, start_col, end_col, ref_seq_col, strand_col), then all required ones must be (chrom_col, start_col, end_col)."
            )

    all_cols = [target_name_col, seq_col] + [
        c for c in optional_cols + masking_cols if c
    ]
    if additional_representative_mhap_cols:
        all_cols += additional_representative_mhap_cols
    all_cols = list(set(all_cols))

    unique_table = (
        microhaplotype_table[all_cols].drop_duplicates().reset_index(drop=True)
    )

    warn_if_duplicated_seqs(unique_table, target_name_col, seq_col)
    mhap_data = {"targets": []}
    for target, group in unique_table.groupby(target_name_col):
        target_dict = {"target_name": target, "microhaplotypes": []}
        first_row = group.iloc[0]
        if chrom_col and pd.notna(group[chrom_col].iloc[0]):
            loc = {
                "genome_id": genome_id,
                "chrom": first_row[chrom_col],
                "start": first_row[start_col],
                "end": first_row[end_col],
            }
            if ref_seq_col and pd.notna(first_row[ref_seq_col]):
                loc["ref_seq"] = first_row[ref_seq_col]
            if strand_col and pd.notna(first_row[strand_col]):
                loc["strand"] = first_row[strand_col]
            target_dict["mhap_location"] = loc

        for _, row in group.iterrows():
            mhap = {"seq": row[seq_col]}
            if val := get_if_present(row, alt_annotations_col):
                mhap["alt_annotations"] = val
            if val := get_if_present(row, microhaplotype_name_col):
                mhap["microhaplotype_name"] = val
            if val := get_if_present(row, pseudocigar_col):
                mhap["pseudo_cigar"] = val
            if val := get_if_present(row, quality_col):
                mhap["quality"] = val
            if additional_representative_mhap_cols:
                for col in additional_representative_mhap_cols:
                    if val := get_if_present(row, col):
                        mhap[col] = val

            # Add masking if present
            masking = extract_masking(row)
            if masking:
                mhap["masking"] = masking

            target_dict["microhaplotypes"].append(mhap)

        mhap_data["targets"].append(target_dict)
    return mhap_data


def create_detected_microhaplotype_dict(
    microhaplotype_table: pd.DataFrame,
    bioinformatics_run_name: str,
    representative_microhaplotype_dict: dict,
    library_sample_name_col: str = "library_sample_name",
    target_name_col: str = "target_name",
    seq_col: str = "seq",
    reads_col: str = "reads",
    umis_col: str | None = None,
    additional_mhap_detected_cols: list | None = None,
):
    """
    Convert the read-in microhaplotype calls table into the detected microhaplotype dictionary.

    :param microhaplotype_table: Parsed microhaplotype calls table.
    :param bioinformatics_run_name:  Unique name for the bioinformatics run that generated the data.
    :param representative_microhaplotype_dict: Dictionary of representative microhaplotypes.
    :param library_sample_name_col: Column containing the sample IDs.
    :param target_name_col: Column containing the locus IDs.
    :param seq_col: Column containing the microhaplotype sequences.
    :param reads_col: Column containing the read counts.
    :param umis_col: : Ccolumn with unique molecular identifier count associated with this microhaplotype
    :param additional_mhap_detected_cols: Optional additional columns to add to the microhaplotypes detected, the key is the pandas column and the value is what to name it in the output.
    :return: A dictionary of detected microhaplotype results.
    """
    # Rename columns in dataframe and gather columns
    column_mapping = {
        library_sample_name_col: "library_sample_name",
        target_name_col: "target_name",
        seq_col: "seq",
        reads_col: "reads",
    }
    mhap_cols = ["mhap_id", "reads"]
    if umis_col:
        column_mapping[umis_col] = "umis"
        mhap_cols.append("umis")
    df = microhaplotype_table.rename(columns=column_mapping).copy()

    # Validate additional columns if provided
    if additional_mhap_detected_cols:
        check_additional_columns_exist(df, additional_mhap_detected_cols)
        mhap_cols += additional_mhap_detected_cols

    # Find IDs for targets and mhaps in the representative table
    df = get_target_id_in_representative_mhaps(df, representative_microhaplotype_dict)
    df = get_mhap_index_in_representative_mhaps(df, representative_microhaplotype_dict)

    # Build detected mhap table
    mhap_detected = build_detected_mhap_dict(df, bioinformatics_run_name, mhap_cols)
    return mhap_detected


def build_detected_mhap_dict(
    df, bioinformatics_run_name, mhap_cols, always_include=None
):
    if always_include is None:
        always_include = ["mhap_id", "reads"]

    mhap_detected = {
        "bioinformatics_run_name": bioinformatics_run_name,
        "library_samples": [],
    }

    for sample, sample_df in df.groupby("library_sample_name"):
        target_results = []
        for target_id, target_df in sample_df.groupby("mhaps_target_id"):
            mhaps = target_df.apply(
                lambda row: {
                    col: row[col]
                    for col in mhap_cols
                    if col in always_include or pd.notna(row[col])
                },
                axis=1,
            ).to_list()
            target_results.append({"mhaps_target_id": target_id, "mhaps": mhaps})
        mhap_detected["library_samples"].append(
            {"library_sample_name": sample, "target_results": target_results}
        )

    return mhap_detected


def get_target_id_in_representative_mhaps(df, representative_dict):
    target_name_to_mhaps_target_id = {
        entry["target_name"]: i
        for i, entry in enumerate(representative_dict["targets"])
    }
    df["mhaps_target_id"] = df["target_name"].map(target_name_to_mhaps_target_id)
    if df["mhaps_target_id"].isnull().any():
        missing_targets = df[df.mhaps_target_id.isnull()]["target_name"].unique()
        raise ValueError(
            f"Missing target_name(s) in representative microhaplotype table: {missing_targets}"
        )
    return df


def get_mhap_index_in_representative_mhaps(df, representative_dict):
    target_seq_to_mhap_id = {
        (target_id, mhap["seq"]): i
        for target_id, target_entry in enumerate(representative_dict["targets"])
        for i, mhap in enumerate(target_entry["microhaplotypes"])
    }
    df["mhap_id"] = df.apply(
        lambda row: target_seq_to_mhap_id.get((row["mhaps_target_id"], row["seq"])),
        axis=1,
    )
    if df["mhap_id"].isnull().any():
        missing_seqs = df[df["mhap_id"].isnull()][
            ["target_name", "seq"]
        ].drop_duplicates()
        raise ValueError(
            f"Some seq values not found in representative microhaplotype table:\n{missing_seqs}"
        )
    return df
