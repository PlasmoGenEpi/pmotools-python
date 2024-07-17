#!/usr/bin/env python3
import pandas as pd


def microhaplotype_table_to_pmo_dict(file, bioinfo_id, sampleID_col, locus_col, mhap_col, reads_col, delim):

    contents = pd.read_csv(file, sep=delim)
    representative_microhaplotype_dict = create_representative_microhaplotype_dict(
        contents, locus_col, mhap_col)

    detected_mhap_dict = create_detected_microhaplotype_dict(contents, sampleID_col, locus_col,
                                                             mhap_col, reads_col, representative_microhaplotype_dict)

    output_data = {"haplotypes_detected": {'bioinformatics_id': bioinfo_id, 'samples': detected_mhap_dict},
                   "representative_haplotype_sequences": representative_microhaplotype_dict}
    return output_data


def create_representative_microhaplotype_dict(microhaplotype_table, locus_col, mhap_col):
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


def create_detected_microhaplotype_dict(microhaplotype_table, sampleID_col, locus_col, mhap_col, reads_col, representative_microhaplotype_dict):
    json_data = {}
    sample_grouped = microhaplotype_table.groupby(sampleID_col)
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
                        "Representatice microhaplotype ids are not unique")
                else:
                    matching_id = matching_ids[0]
                microhaplotypes.append({
                    "haplotype_id": matching_id,
                    "read_count": row[reads_col]
                })
            target_results[locus] = {
                "microhaplotypes": microhaplotypes,
            }
        json_data[sample_id] = {
            "sample_id": sample_id,
            "target_results": target_results
        }
        return json_data
