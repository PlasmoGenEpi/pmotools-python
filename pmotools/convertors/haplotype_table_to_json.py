import os
import argparse
import json
import pandas as pd


def parse_args_excel_meta_to_json_meta():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, required=True,
                        help='Input excel file path')
    parser.add_argument('--bioinfo_id', type=str, required=True,
                        help='Identifier of bioinformatics processing run')
    parser.add_argument('--sampleID_col', type=str,
                        default='sampleID', help='Column name containing sampleIDs')
    parser.add_argument('--locus_col', type=str, default='locus',
                        help='Column name containing locus information')
    parser.add_argument('--mhap_col', type=str, default='asv',
                        help='Column name containing microhaplotypes')
    parser.add_argument('--reads_col', type=str, default='reads',
                        help='Column name containing reads per microhaplotype')
    parser.add_argument('--delim', type=str, default='\t',
                        help='Delimiter of input file')
    parser.add_argument('--output', type=str, required=True,
                        help='Output json file path')
    parser.add_argument('--overwrite', action='store_true',
                        help='If output file exists, overwrite it')
    return parser.parse_args()


def microhaplotype_table_to_json():
    args = parse_args_excel_meta_to_json_meta()
    if not args.output.endswith('.json'):
        args.output += '.json'
    # make sure file exists
    if not os.path.exists(args.file):
        raise FileNotFoundError(args.file)
    # only overwrite an existing file if --overwrite is on
    if os.path.exists(args.output) and not args.overwrite:
        raise Exception(
            "Output file already exists, use --overWrite to overwrite it")

    contents = pd.read_csv(args.file, sep=args.delim)
    representative_microhaplotype_dict = create_representative_microhaplotype_dict(
        contents, args.locus_col, args.mhap_col)

    detected_mhap_dict = create_detected_microhaplotype_dict(contents, args.sampleID_col, args.locus_col,
                                                             args.mhap_col, args.reads_col, representative_microhaplotype_dict)

    output_data = {"haplotypes_detected": {'bioinformatics_id': args.bioinfo_id, 'samples': detected_mhap_dict},
                   "representative_haplotype_sequences": representative_microhaplotype_dict}

    # Write output as json
    json_str = json.dumps(output_data, indent=4)
    with open(args.output, 'w') as json_file:
        json_file.write(json_str)


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


if __name__ == "__main__":
    microhaplotype_table_to_json()
