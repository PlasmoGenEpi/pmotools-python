#!/usr/bin/env python3
import argparse
import json
from collections import defaultdict

import pandas as pd
from pmotools.utils.small_utils import Utils


def parse_args_terra_amp_output_to_json():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True, help="Input excel file path")
    parser.add_argument(
        "--gt_sheet",
        type=str,
        default="gt",
        required=False,
        help="The gt sheet to convert, if none provided will default to gt",
    )
    parser.add_argument(
        "--asv_table_sheet",
        type=str,
        default="asv_table",
        required=False,
        help="The asv_table sheet to convert, if none provided will default to asv_table",
    )
    parser.add_argument(
        "--asv_seqs_sheet",
        type=str,
        default="asv_seqs",
        required=False,
        help="The asv_seqs sheet to convert, if none provided will default to asv_seqs",
    )

    # parser.add_argument('--index_col_name', type=str, required=False, help='by default output is a list, if an index column name is supplied it will be a dict with this column as index')

    parser.add_argument(
        "--output", type=str, required=True, help="Output json file path"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="If output file exists, overwrite it"
    )
    return parser.parse_args()


def terra_amp_output_to_json():
    args = parse_args_terra_amp_output_to_json()
    args.output = Utils.appendStrAsNeeded(args.output, ".json")
    # check if input file exists and if output file exists check if --overwrite flag is set
    Utils.inputOutputFileCheckFromArgParse(args)

    # print(asv_table_look_up)
    representative_microhaplotype_id = "terra"
    tar_amp_bioinformatics_info_id = "terra"
    asv_seqs = pd.read_excel(
        args.file, sheet_name=args.asv_seqs_sheet, index_col="asv_id"
    )
    asv_seq_table_look_up = dict()
    for index, row in asv_seqs.iterrows():
        asv_seq_table_look_up[index] = row["asv_seq"]

    representative_ref_seqs = {
        "representative_microhaplotype_id": representative_microhaplotype_id,
        "targets": {},
    }
    asv_table = pd.read_excel(
        args.file, sheet_name=args.asv_table_sheet, index_col="hapid"
    )
    asv_table_look_up = dict()
    asv_table_look_up_by_target_by_cigar = defaultdict(lambda: defaultdict(str))

    for index, row in asv_table.iterrows():
        asv_table_look_up.update(
            {
                index: {
                    "asv_id": index,
                    "CIGAR": row["CIGAR"],
                    "CIGAR_masked": row["CIGAR_masked"],
                    "Amplicon": row["Amplicon"],
                }
            }
        )
        asv_table_look_up_by_target_by_cigar[row["Amplicon"]][
            row["CIGAR_masked"]
        ] = index
        if row["Amplicon"] not in representative_ref_seqs["targets"]:
            representative_ref_seqs["targets"][row["Amplicon"]] = {
                "target_id": row["Amplicon"],
                "seqs": {},
            }
        CIGAR_masked = row["CIGAR_masked"]
        if isinstance(row["CIGAR_masked"], float) and pd.isna(row["CIGAR_masked"]):
            CIGAR_masked = "."

        representative_ref_seqs["targets"][row["Amplicon"]]["seqs"].update(
            {
                index: {
                    "microhaplotype_id": index,
                    "seq": asv_seq_table_look_up[index],
                    "CIGAR": row["CIGAR"],
                    "CIGAR_masked": CIGAR_masked,
                }
            }
        )

    microhaplotypes_detected = {
        "representative_microhaplotype_id": representative_microhaplotype_id,
        "tar_amp_bioinformatics_info_id": tar_amp_bioinformatics_info_id,
        "experiment_samples": {},
    }

    gt_contents = pd.read_excel(
        args.file, sheet_name=args.gt_sheet, index_col="Sample_id"
    )
    warnings = []

    for index, row in gt_contents.iterrows():
        sample_data = {"experiment_sample_id": index, "target_results": {}}

        for amplicon_id in row.keys():
            amplicon_data = []
            # print(row[amplicon_id])
            # print(type(row[amplicon_id]))
            # print(str(row[amplicon_id]))
            # # print(row[amplicon_id].isna())
            # print(isinstance(row[amplicon_id], float))
            # print(isinstance(row[amplicon_id], float) and pd.isna(row[amplicon_id]))
            if not (isinstance(row[amplicon_id], float) and pd.isna(row[amplicon_id])):
                current_haps = row[amplicon_id].split("_")
                for current_hap in current_haps:
                    current_hap_split = current_hap.split(":")
                    if len(current_hap_split) != 2:
                        warnings.append(
                            "sample: "
                            + str(index)
                            + " target: "
                            + amplicon_id
                            + "failed to process "
                            + current_hap
                            + " expected two values separated by :"
                        )
                    else:
                        read_cnt = float(current_hap_split[1])
                        cigar = current_hap_split[0]
                        if (
                            cigar
                            not in asv_table_look_up_by_target_by_cigar[amplicon_id]
                        ):
                            warnings.append(
                                "sample: "
                                + str(index)
                                + " target: "
                                + amplicon_id
                                + "failed to find "
                                + cigar
                                + " in cigar look up"
                            )
                        else:
                            amplicon_data.append(
                                {
                                    "microhaplotype_id": asv_table_look_up_by_target_by_cigar[
                                        amplicon_id
                                    ][cigar],
                                    "read_count": read_cnt,
                                }
                            )
            if len(amplicon_data) > 0:
                sample_data["target_results"].update(
                    {
                        amplicon_id: {
                            "target_id": amplicon_id,
                            "microhaplotypes": amplicon_data,
                        }
                    }
                )
        microhaplotypes_detected["experiment_samples"].update({index: sample_data})
    if len(warnings) > 0:
        raise Exception("\n".join(warnings))
    output_json = {
        "microhaplotypes_detected": {
            tar_amp_bioinformatics_info_id: microhaplotypes_detected
        },
        "representative_microhaplotype_sequences": {
            representative_microhaplotype_id: representative_ref_seqs
        },
    }

    json.dump(output_json, open(args.output, "w"), indent=4)


if __name__ == "__main__":
    terra_amp_output_to_json()
