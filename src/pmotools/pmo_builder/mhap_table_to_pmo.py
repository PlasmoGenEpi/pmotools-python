#!/usr/bin/env python3
import pandas as pd
import warnings

from ..pmo_builder.json_convert_utils import check_additional_columns_exist


def mhap_table_to_pmo(
    microhaplotype_table: pd.DataFrame,
    bioinformatics_run_name: str | None = None,
    library_sample_name_col: str = "library_sample_name",
    target_name_col: str = "target_name",
    seq_col: str = "seq",
    reads_col: str = "reads",
    genome_id: int = 0,
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
    pseudocigar_chrom_col: str | None = None,
    pseudocigar_start_col: str | None = None,
    pseudocigar_end_col: str | None = None,
    pseudocigar_ref_seq_col: str | None = None,
    pseudocigar_strand_col: str | None = None,
    pseudocigar_genome_id: int | None = None,
    pseudocigar_generation_description_col: str | None = None,
    quality_col: str | None = None,
    additional_representative_mhap_cols: list | None = None,
    additional_mhap_detected_cols: list | None = None,
):
    """
    Convert a dataframe of microhaplotype calls into a dictionary containing a dictionary for the haplotypes_detected and a dictionary for the representative_haplotype_sequences.

    :param microhaplotype_table: the dataframe containing microhaplotype calls
    :type microhaplotype_table: pd.DataFrame
    :param bioinformatics_run_name: unique name for the bioinformatics run that generated the data (column name or individual run name). Default: None
    :type bioinformatics_run_name: str, optional
    :param library_sample_name_col: the name of the column containing the library sample names. Default: library_sample_name
    :type library_sample_name_col: str
    :param target_name_col: the name of the column containing the targets. Default: target_name
    :type target_name_col: str
    :param seq_col: the name of the column containing the microhaplotype sequences. Default: seq
    :type seq_col: str
    :param reads_col: the name of the column containing the read counts. Default: reads
    :type reads_col: str
    :param genome_id: the ID of the genome used as reference. Default: None
    :type genome_id: int, optional
    :param umis_col: the name of the column with the unique molecular identifier count associated with this microhaplotype
    :type umis_col: str, optional
    :param chrom_col: the name of the column containing the chromosome name of the microhaplotype
    :type chrom_col: str, optional
    :param start_col: the name of the column containing the start of the microhaplotype
    :type start_col: str, optional
    :param end_col: the name of the column containing the end of the microhaplotype
    :type end_col: str, optional
    :param ref_seq_col: the name of the column containing the reference sequence for the microhaplotype
    :type ref_seq_col: str, optional
    :param strand_col: the name of the column containing the strand of the microhaplotype
    :type strand_col: str, optional
    :param alt_annotations_col: the name of the column containing any alternative annotations
    :type alt_annotations_col: str, optional
    :param masking_seq_start_col: the name of the column containing a list of start positions for masking
    :type masking_seq_start_col: str, optional
    :param masking_seq_segment_size_col: the name of the column containing a list of lengths of the segments in seq being masked
    :type masking_seq_segment_size_col: str, optional
    :param masking_replacement_size_col: the name of the column containing a list of lengths of the masking replacements
    :type masking_replacement_size_col: str, optional
    :param masking_delim: delimiter of the masking information. Default: ','
    :type masking_delim: str, optional
    :param microhaplotype_name_col: the name of the column containing an optional name for this microhaplotype
    :type microhaplotype_name_col: str, optional
    :param pseudocigar_col: the name of the column containing a pseudocigar for the microhaplotype
    :type pseudocigar_col: str, optional
    :param pseudocigar_chrom_col: column for the pseudocigar ref_loc chromosome; defaults to chrom_col
    :type pseudocigar_chrom_col: str, optional
    :param pseudocigar_start_col: column for the pseudocigar ref_loc start (required if pseudocigar_col is set)
    :type pseudocigar_start_col: str, optional
    :param pseudocigar_end_col: column for the pseudocigar ref_loc end (required if pseudocigar_col is set)
    :type pseudocigar_end_col: str, optional
    :param pseudocigar_ref_seq_col: optional column for the pseudocigar ref_loc reference sequence
    :type pseudocigar_ref_seq_col: str, optional
    :param pseudocigar_strand_col: optional column for the pseudocigar ref_loc strand
    :type pseudocigar_strand_col: str, optional
    :param pseudocigar_genome_id: genome id for the pseudocigar ref_loc; defaults to genome_id
    :type pseudocigar_genome_id: int, optional
    :param pseudocigar_generation_description_col: optional column describing how the pseudocigar was generated
    :type pseudocigar_generation_description_col: str, optional
    :param quality_col: the name of the column containing the ANSI FASTQ per-base quality score for this sequence
    :type quality_col: str, optional
    :param additional_representative_mhap_cols: additional columns to add to the representative microhaplotypes table
    :type additional_representative_mhap_cols: list of str, optional
    :param additional_mhap_detected_cols: additional columns to add to the detected microhaplotypes table
    :type additional_mhap_detected_cols: list of str, optional
    :return: a dict of both the haplotypes_detected and representative_haplotype_sequences
    :rtype: dict
    """

    representative_microhaplotype_dict = create_representative_microhaplotype_dict(
        microhaplotype_table=microhaplotype_table,
        target_name_col=target_name_col,
        seq_col=seq_col,
        genome_id=genome_id,
        chrom_col=chrom_col,
        start_col=start_col,
        end_col=end_col,
        ref_seq_col=ref_seq_col,
        strand_col=strand_col,
        alt_annotations_col=alt_annotations_col,
        masking_seq_start_col=masking_seq_start_col,
        masking_seq_segment_size_col=masking_seq_segment_size_col,
        masking_replacement_size_col=masking_replacement_size_col,
        masking_delim=masking_delim,
        microhaplotype_name_col=microhaplotype_name_col,
        pseudocigar_col=pseudocigar_col,
        pseudocigar_chrom_col=pseudocigar_chrom_col,
        pseudocigar_start_col=pseudocigar_start_col,
        pseudocigar_end_col=pseudocigar_end_col,
        pseudocigar_ref_seq_col=pseudocigar_ref_seq_col,
        pseudocigar_strand_col=pseudocigar_strand_col,
        pseudocigar_genome_id=pseudocigar_genome_id,
        pseudocigar_generation_description_col=pseudocigar_generation_description_col,
        quality_col=quality_col,
        additional_representative_mhap_cols=additional_representative_mhap_cols,
    )

    detected_mhap_dict_list = []
    if bioinformatics_run_name in microhaplotype_table.columns:
        for bioinfo_run in microhaplotype_table[bioinformatics_run_name].unique():
            microhaplotype_table_per_run = microhaplotype_table[
                microhaplotype_table[bioinformatics_run_name] == bioinfo_run
            ]
            detected_mhap_dict = create_detected_microhaplotype_dict(
                microhaplotype_table=microhaplotype_table_per_run,
                representative_microhaplotype_dict=representative_microhaplotype_dict,
                bioinformatics_run_name=bioinfo_run,
                library_sample_name_col=library_sample_name_col,
                target_name_col=target_name_col,
                seq_col=seq_col,
                reads_col=reads_col,
                umis_col=umis_col,
                additional_mhap_detected_cols=additional_mhap_detected_cols,
            )
            detected_mhap_dict_list.append(detected_mhap_dict)
    else:
        detected_mhap_dict = create_detected_microhaplotype_dict(
            microhaplotype_table=microhaplotype_table,
            representative_microhaplotype_dict=representative_microhaplotype_dict,
            bioinformatics_run_name=bioinformatics_run_name,
            library_sample_name_col=library_sample_name_col,
            target_name_col=target_name_col,
            seq_col=seq_col,
            reads_col=reads_col,
            umis_col=umis_col,
            additional_mhap_detected_cols=additional_mhap_detected_cols,
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
    pseudocigar_chrom_col: str | None = None,
    pseudocigar_start_col: str | None = None,
    pseudocigar_end_col: str | None = None,
    pseudocigar_ref_seq_col: str | None = None,
    pseudocigar_strand_col: str | None = None,
    pseudocigar_genome_id: int | None = None,
    pseudocigar_generation_description_col: str | None = None,
    quality_col: str | None = None,
    additional_representative_mhap_cols: list[str] | None = None,
):
    """
    Convert the read-in microhaplotype calls table into a representative microhaplotype JSON-like dictionary.

    :param microhaplotype_table: the dataframe containing microhaplotype calls
    :type microhaplotype_table: pd.DataFrame
    :param target_name_col: the name of the column containing the targets. Default: target_name
    :type target_name_col: str
    :param seq_col: the name of the column containing the microhaplotype sequences. Default: seq
    :type seq_col: str
    :param genome_id: the genome ID
    :type genome_id: int
    :param chrom_col: the name of the column containing the chromosome name of the microhaplotype
    :type chrom_col: str, optional
    :param start_col: the name of the column containing the start of the microhaplotype
    :type start_col: str, optional
    :param end_col: the name of the column containing the end of the microhaplotype
    :type end_col: str, optional
    :param ref_seq_col: the name of the column containing the reference sequence for the microhaplotype
    :type ref_seq_col: str, optional
    :param strand_col: the name of the column containing the strand of the microhaplotype
    :type strand_col: str, optional
    :param alt_annotations_col: the name of the column containing any alternative annotations
    :type alt_annotations_col: str, optional
    :param masking_seq_start_col: the name of the column containing a list of start positions for masking
    :type masking_seq_start_col: str, optional
    :param masking_seq_segment_size_col: the name of the column containing a list of lengths of the segments in seq being masked
    :type masking_seq_segment_size_col: str, optional
    :param masking_replacement_size_col: the name of the column containing a list of lengths of the masking replacements
    :type masking_replacement_size_col: str, optional
    :param masking_delim: delimiter of the masking information. Default: ','
    :type masking_delim: str, optional
    :param microhaplotype_name_col: the name of the column containing an optional name for this microhaplotype
    :type microhaplotype_name_col: str, optional
    :param pseudocigar_col: the name of the column containing a pseudocigar for the microhaplotype
    :type pseudocigar_col: str, optional
    :param pseudocigar_chrom_col: column for the pseudocigar ref_loc chromosome; defaults to chrom_col
    :type pseudocigar_chrom_col: str, optional
    :param pseudocigar_start_col: column for the pseudocigar ref_loc start (required if pseudocigar_col is set)
    :type pseudocigar_start_col: str, optional
    :param pseudocigar_end_col: column for the pseudocigar ref_loc end (required if pseudocigar_col is set)
    :type pseudocigar_end_col: str, optional
    :param pseudocigar_ref_seq_col: optional column for the pseudocigar ref_loc reference sequence
    :type pseudocigar_ref_seq_col: str, optional
    :param pseudocigar_strand_col: optional column for the pseudocigar ref_loc strand
    :type pseudocigar_strand_col: str, optional
    :param pseudocigar_genome_id: genome id for the pseudocigar ref_loc; defaults to genome_id
    :type pseudocigar_genome_id: int, optional
    :param pseudocigar_generation_description_col: optional column describing how the pseudocigar was generated
    :type pseudocigar_generation_description_col: str, optional
    :param quality_col: the name of the column containing the ANSI FASTQ per-base quality score for this sequence
    :type quality_col: str, optional
    :param additional_representative_mhap_cols: additional columns to add to the representative microhaplotypes table
    :type additional_representative_mhap_cols: list of str, optional
    :return: a dictionary formatted for JSON output with representative microhaplotype sequences
    :rtype: dict
    """

    if additional_representative_mhap_cols:
        check_additional_columns_exist(
            microhaplotype_table, additional_representative_mhap_cols
        )

    # the pseudocigar ref_loc shares the chromosome / genome with the
    # microhaplotype location by default, but each may use its own column
    pseudocigar_chrom_col = (
        pseudocigar_chrom_col if pseudocigar_chrom_col else chrom_col
    )
    ps_genome_id = (
        pseudocigar_genome_id if pseudocigar_genome_id is not None else genome_id
    )
    if pseudocigar_col and not (
        pseudocigar_chrom_col and pseudocigar_start_col and pseudocigar_end_col
    ):
        raise ValueError(
            "pseudocigar_col is set, so a Pseudocigar ref_loc must be "
            "constructable: set a chromosome (pseudocigar_chrom_col or "
            "chrom_col), pseudocigar_start_col, and pseudocigar_end_col."
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
        pseudocigar_chrom_col,
        pseudocigar_start_col,
        pseudocigar_end_col,
        pseudocigar_ref_seq_col,
        pseudocigar_strand_col,
        pseudocigar_generation_description_col,
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
                if not (
                    pd.notna(row[pseudocigar_chrom_col])
                    and pd.notna(row[pseudocigar_start_col])
                    and pd.notna(row[pseudocigar_end_col])
                ):
                    raise ValueError(
                        f"pseudocigar present for target {target}, seq "
                        f"{row[seq_col]} but its ref_loc chrom/start/end "
                        "is missing"
                    )
                ref_loc = {
                    "genome_id": ps_genome_id,
                    "chrom": row[pseudocigar_chrom_col],
                    "start": row[pseudocigar_start_col],
                    "end": row[pseudocigar_end_col],
                }
                if pseudocigar_ref_seq_col and pd.notna(row[pseudocigar_ref_seq_col]):
                    ref_loc["ref_seq"] = row[pseudocigar_ref_seq_col]
                if pseudocigar_strand_col and pd.notna(row[pseudocigar_strand_col]):
                    ref_loc["strand"] = row[pseudocigar_strand_col]
                pseudocigar = {"pseudocigar_seq": val, "ref_loc": ref_loc}
                if pseudocigar_generation_description_col and pd.notna(
                    row[pseudocigar_generation_description_col]
                ):
                    pseudocigar["pseudocigar_generation_description"] = row[
                        pseudocigar_generation_description_col
                    ]
                mhap["pseudocigar"] = pseudocigar
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
    representative_microhaplotype_dict: dict,
    bioinformatics_run_name: str | None = None,
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
    :param representative_microhaplotype_dict: Dictionary of representative microhaplotypes.
    :param bioinformatics_run_name: Optional Unique name for the bioinformatics run that generated the data.
    :param library_sample_name_col: Column containing the sample IDs.
    :param target_name_col: Column containing the locus IDs.
    :param seq_col: Column containing the microhaplotype sequences.
    :param reads_col: Column containing the read counts.
    :param umis_col: Optional Column with unique molecular identifier count associated with this microhaplotype
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
        "library_samples": [],
    }
    if bioinformatics_run_name is not None:
        mhap_detected["bioinformatics_run_name"] = bioinformatics_run_name

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


def create_minimum_library_specimen_dict_from_mhap_table(
    detected_microhaps: list[dict],
    panel_name: str,
    library_sample_field_name: str = "library_sample_name",
    library_sample_specimen_key: dict[str, str] | pd.DataFrame | None = None,
    library_sample_name_col: str = "library_sample_name",
    specimen_name_col: str = "specimen_name",
    missing_library_sample_becomes_specimen_name: bool = False,
):
    """
    Create a minimum library_sample_info and specimen_info dicts from the detected microhaps

    :param detected_microhaps: the detected microhaps object created by create_detected_microhaplotype_dict
    :param panel_name: the panel_name for the library_sample
    :param library_sample_field_name: the field name to use to extract the library_sample_name from the detected_michrohaplotypes
    :param library_sample_specimen_key: a dict mapping library_sample_name -> specimen_name,
                                        or a pandas DataFrame with two columns for renaming controlled by library_sample_name_col and specimen_name_col
                                        if None, specimen_name == library_sample_name
    :param library_sample_name_col: the column name in library_sample_specimen_key that contains the library_sample_name
    :param specimen_name_col: the column name in library_sample_specimen_key that contains the specimen_name
    :param missing_library_sample_becomes_specimen_name: if True and a library_sample_name is missing
                                                         from library_sample_specimen_key, fall back to
                                                         using the library_sample_name as the specimen_name;
                                                         if False, raise an error
    :return: dict with keys 'library_sample_info' and 'specimen_info'
    """
    # Collect all sample dicts across every entry in detected_microhaps
    all_samples: list[dict] = []
    for entry in detected_microhaps:
        all_samples.extend(entry.get("library_samples", []))

    # check that every sample has the expected key
    missing_key_indices = [
        i for i, s in enumerate(all_samples) if library_sample_field_name not in s
    ]
    if missing_key_indices:
        raise KeyError(
            f"The following sample indices are missing the field name '{library_sample_field_name}': "
            f"{missing_key_indices}"
        )

    # check that all library_sample_name values are unique
    raw_names: list[str] = [s[library_sample_field_name] for s in all_samples]
    seen: set[str] = set()
    duplicates: set[str] = set()
    for name in raw_names:
        if name in seen:
            duplicates.add(name)
        seen.add(name)
    if duplicates:
        raise ValueError(f"Duplicate library sample names found: {sorted(duplicates)}")
    actual_library_sample_specimen_key = None
    if library_sample_specimen_key is not None and isinstance(
        library_sample_specimen_key, dict
    ):
        actual_library_sample_specimen_key = library_sample_specimen_key
    elif library_sample_specimen_key is not None and isinstance(
        library_sample_specimen_key, pd.DataFrame
    ):
        actual_library_sample_specimen_key = library_sample_specimen_key.set_index(
            library_sample_name_col
        )[specimen_name_col].to_dict()
    # now construct library_sample_info
    library_sample_info: list[dict] = []
    for sample in all_samples:
        lib_name: str = sample[library_sample_field_name]
        # use look up table to get specimen_name if provided, otherwise use library_sample_name as specimen_name
        if actual_library_sample_specimen_key is not None:
            if lib_name in actual_library_sample_specimen_key:
                specimen_name = actual_library_sample_specimen_key[lib_name]
            elif missing_library_sample_becomes_specimen_name:
                # if not in key but allowing missing to become specimen_name, use library_sample_name as specimen_name
                specimen_name = lib_name
            else:
                raise KeyError(
                    f"library_sample_name '{lib_name}' not found in library_sample_specimen_key "
                    f"and missing_library_sample_becomes_specimen_name is False."
                )
        else:
            specimen_name = lib_name

        library_sample_info.append(
            {
                "library_sample_name": lib_name,
                "panel_name": panel_name,
                "specimen_name": specimen_name,
            }
        )

    # build specimen_info from unique specimen_names (preserving first-seen order)
    seen_specimens: set[str] = set()
    specimen_info: list[dict] = []
    for entry in library_sample_info:
        sp = entry["specimen_name"]
        if sp not in seen_specimens:
            seen_specimens.add(sp)
            specimen_info.append({"specimen_name": sp})

    return {
        "library_sample_info": library_sample_info,
        "specimen_info": specimen_info,
    }
