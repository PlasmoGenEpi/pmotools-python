#!/usr/bin/env python3
import pandas as pd


def microhaplotype_table_to_pmo_dict(file : str, bioinfo_id : str, sampleID_col : str, locus_col : str, mhap_col : str, reads_col : str, delim : str, additional_hap_detected_cols : dict | None = None):
    """
    Convert a delimited file of a microhaplotype calls into a dictionary containing a dictionary for the haplotypes_detected and a dictionary for the representative_haplotype_sequences


    :param file: The name of the microhaplotype table file
    :param bioinfo_id: the bioinformatics ID of the microhaplotype table
    :param sampleID_col: the name of the column containing the sample IDs
    :param locus_col: the name of the column containing the locus IDs
    :param mhap_col: the name of the column containing the microhaplotype sequence
    :param reads_col: the name of the column containing the reads counts
    :param delim: the delimiter character of the input file
    :param additional_hap_detected_cols: optional additional columns to add to the microhaplotype detected dictionary, the key is the pandas column and the value is what to name it in the output
    :return: a dict of both the haplotypes_detected and representative_haplotype_sequences
    """

    contents = pd.read_csv(file, sep=delim)

    representative_microhaplotype_dict = create_representative_microhaplotype_dict(
        contents, locus_col, mhap_col)

    detected_mhap_dict = create_detected_microhaplotype_dict(contents, sampleID_col, locus_col,
                                                             mhap_col, reads_col, representative_microhaplotype_dict,
                                                             additional_hap_detected_cols)

    output_data = {"haplotypes_detected": {'bioinformatics_id': bioinfo_id, 'samples': detected_mhap_dict},
                   "representative_haplotype_sequences": representative_microhaplotype_dict}
    return output_data


def create_representative_microhaplotype_dict(microhaplotype_table : pd.DataFrame, locus_col : str, mhap_col : str):
    """
    Convert the read in microhaplotype calls table into the representative microhaplotype dictionary

    :param microhaplotype_table: the parsed microhaplotype calls table
    :param locus_col: the name of the column containing the locus IDs
    :param mhap_col: the name of the column containing the microhaplotype sequence
    :return: a dictionary of the representative microhaplotype sequences
    """
    microhaplotype_table = microhaplotype_table[[locus_col,
                                                 mhap_col]].drop_duplicates()
    microhaplotype_table.reset_index(inplace=True, drop=True)

    # Group the dataframe by 'locus'
    grouped = microhaplotype_table.groupby('locus')
    json_data = {}
    # Populate the dictionary
    for locus, group in grouped:
        microhaplotypes = []
        microhaplotype_index = 0
        for _, row in group.iterrows():
            microhaplotypes.append({
                "microhaplotype_id": '.'.join([locus, str(microhaplotype_index)]),
                "seq": row[mhap_col]
            })
            microhaplotype_index += 1
        json_data[locus] = {
            "seqs": microhaplotypes}
    return json_data

# TODO: add on additional columns


def create_detected_microhaplotype_dict(microhaplotype_table : pd.DataFrame, sampleID_col : str,
                                        locus_col : str, mhap_col : str,
                                        reads_col : str, representative_microhaplotype_dict : dict,
                                        additional_hap_detected_cols : dict | None = None):
    """
    Convert the read in microhaplotype calls table into the detected microhaplotype dictionary


    :param microhaplotype_table: the parsed microhaplotype calls table
    :param sampleID_col: the name of the column containing the sample IDs
    :param locus_col: the name of the column containing the locus IDs
    :param mhap_col: the name of the column containing the microhaplotype sequence
    :param reads_col: the name of the column containing the reads counts
    :param representative_microhaplotype_dict: the representative microhaplotype dictionary created by the function create_representative_microhaplotype_dict
    :param additional_hap_detected_cols: optional additional columns to add to the microhaplotype detected dictionary, the key is the pandas column and the value is what to name it in the output
    :return: a dictionary of the detected microhaplotype results
    """
    json_data = {}
    sample_grouped = microhaplotype_table.groupby(sampleID_col)

    # check if adding additional columns and if you are then make sure the data table has them
    if additional_hap_detected_cols is not None:
        not_found_cols = []
        for col in additional_hap_detected_cols:
            if col not in microhaplotype_table.columns:
                not_found_cols.append(col)
        if len(not_found_cols) > 0:
            raise ValueError("Could not find the following columns: " + ",".join(not_found_cols), " when trying to add additional columns")

    for sample_id, sample_group in sample_grouped:
        target_results = {}
        locus_grouped = sample_group.groupby(locus_col)
        for locus, locus_group in locus_grouped:
            microhaplotypes = []
            representative_microhaplotype_for_locus = representative_microhaplotype_dict[
                locus]['seqs']
            for _, row in locus_group.iterrows():
                matching_ids = [item['microhaplotype_id']
                                for item in representative_microhaplotype_for_locus if item['seq'] == row[mhap_col]]
                if len(matching_ids) != 1:
                    raise Exception(
                        "Representative microhaplotype ids are not unique")
                else:
                    matching_id = matching_ids[0]
                haplotype_info = {"haplotype_id": matching_id, "read_count": row[reads_col]}

                if additional_hap_detected_cols is not None:
                    for col in additional_hap_detected_cols:
                        haplotype_info[additional_hap_detected_cols[col]] = row[col]
                microhaplotypes.append(haplotype_info)
            target_results[locus] = {
                "microhaplotypes": microhaplotypes,
            }
        json_data[sample_id] = {
            "sample_id": sample_id,
            "target_results": target_results
        }
    return json_data
