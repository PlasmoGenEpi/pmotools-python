#!/usr/bin/env python3
import os, argparse, json
import sys
from collections import defaultdict

import pandas as pd

from pmotools.extract_from_pmo.PMOExtractor import PMOExtractor
from pmotools.extract_from_pmo.PMOReader import PMOReader
from pmotools.extract_from_pmo.PMOWriter import PMOWriter
from pmotools.utils.small_utils import Utils


def parse_args_extract_pmo_with_selected_meta():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, required=True, help='PMO file')
    parser.add_argument('--output', type=str, required=True, help='Output json file path')
    parser.add_argument('--overwrite', action = 'store_true', help='If output file exists, overwrite it')
    parser.add_argument('--verbose', action = 'store_true', help='write out various messages about extraction')
    parser.add_argument('--metaFieldsValues', type=str, required=True, help='Meta Fields to include, should either be a table with columns field, values (and optionally group) or supplied command line as field1:value1,value2,value3;field2:value1,value2')
    return parser.parse_args()

def extract_pmo_with_selected_meta():
    args = parse_args_extract_pmo_with_selected_meta()

    # check files
    Utils.inputOutputFileCheck(args.file, args.output, args.overwrite)

    selected_meta_groups = {}
    # parse meta values
    if os.path.exists(args.metaFieldsValues):
        selected_meta_groups = defaultdict(dict)
        meta_tab = pd.read_csv(args.metaFieldsValues, sep='\t')
        if "field" not in meta_tab or "values" not in meta_tab:
            raise Exception(args.metaFieldsValues + " doesn't have columns field and values, has " + ",".join(meta_tab.columns))
        if "group" in meta_tab:

            for index, row in meta_tab.iterrows():
                values_toks = row["values"].split(",")
                selected_meta_groups[row["group"]][row["field"]] = values_toks
        else:
            group_criteria = {}
            for index, row in meta_tab.iterrows():
                values_toks = row["values"].split(",")
                group_criteria[row["field"]] = values_toks
            selected_meta_groups[0] = group_criteria
    else:
        group_toks = args.metaFieldsValues.split(";")
        for idx, group_tok in enumerate(group_toks):
            group_criteria = {}
            field_with_values_toks = group_tok.split(":")
            for field_with_values_tok in field_with_values_toks:
                field_values_toks = field_with_values_tok.split("=")
                if len(field_values_toks) != 2:
                    raise Exception("error processing " + group_tok, " should be field and values separated by =")
                values_toks = field_values_toks[1].split(",")
                group_criteria[field_values_toks[0]] = values_toks
            selected_meta_groups[idx] = group_criteria

    # read in pmo
    pmo = PMOReader.read_in_pmo(args.file)

    # get count of fields
    fields_counts = PMOExtractor.count_specimen_meta_fields(pmo)

    # check to see if the fields supplied actually exit
    warnings = []
    fields_found = fields_counts["field"].tolist()
    for group in selected_meta_groups.values():
        for field in group.keys():
            if field not in fields_found:
                warnings.append("missing the field: " + field + " in pmo file: " + args.file)
    if len(warnings) > 0:
        raise Exception("\n".join(warnings))


    group_counts = defaultdict(int)
    all_specimen_ids = []
    for specimen in pmo["specimen_infos"].values():
        for group_name, meta in selected_meta_groups.items():
            passes_criteria = True
            for field, values in meta.items():
                if not (field in specimen and str(specimen[field]) in values):
                    passes_criteria = False
                    break
            if passes_criteria:
                group_counts[group_name] += 1
                specimen_id = specimen["specimen_id"]
                all_specimen_ids.append(specimen_id)


    pmo_out = PMOExtractor.extract_from_pmo_select_specimens(pmo, all_specimen_ids)

    # write out the extracted
    PMOWriter.write_out_pmo(pmo_out, args.output, args.overwrite)

    if args.verbose:
        sys.stdout.write("Extracted the following number of specimens per group:" + "\n")
        group_counts_df = pd.DataFrame(list(group_counts.items()), columns=["group", "counts"])
        group_counts_df.to_csv(sys.stdout, sep = "\t", index = False)

if __name__ == "__main__":
    extract_pmo_with_selected_meta()

