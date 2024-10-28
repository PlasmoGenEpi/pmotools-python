#!/usr/bin/env python3



import sys, os, argparse, json
import pandas as pd
sys.path.append(os.path.join(os.path.dirname(__file__), "convertors_to_pmo"))
sys.path.append(os.path.join(os.path.dirname(__file__), "extractors_from_pmo"))
sys.path.append(os.path.join(os.path.dirname(__file__), "pmo_utils"))
sys.path.append(os.path.join(os.path.dirname(__file__), "extract_info_from_pmo"))

from text_meta_to_json_meta import text_meta_to_json_meta
from excel_meta_to_json_meta import excel_meta_to_json_meta
from extract_pmo_with_selected_meta import extract_pmo_with_selected_meta
from extract_pmo_with_select_specimens import extract_pmo_with_select_specimens
from extract_pmo_with_select_targets import extract_pmo_with_select_targets
from microhaplotype_table_to_json_file import microhaplotype_table_to_json_file
from extract_allele_table import extract_for_allele_table
from combine_pmos import combine_pmos
from list_specimen_meta_fields import list_specimen_meta_fields
from list_tar_amp_bioinformatics_info_ids import list_tar_amp_bioinformatics_info_ids
from count_specimen_meta import count_specimen_meta
from count_targets_per_sample import count_targets_per_sample
from count_samples_per_target import count_samples_per_target
from terra_amp_output_to_json import terra_amp_output_to_json



from pmotools.utils.color_text import ColorText as CT

class pmofunction:
    def __init__(self, func, shortDescription):
        self.func = func
        self.shortDescription = shortDescription

class pmotools_runner   :
    def __init__(self) :
        self.functions = {
            "convertors_to_json" : {
                "text_meta_to_json_meta": pmofunction(text_meta_to_json_meta, "Convert text file meta to JSON Meta"),
                "excel_meta_to_json_meta": pmofunction(excel_meta_to_json_meta, "Convert excel file meta to JSON Meta"),
                "microhaplotype_table_to_json_file": pmofunction(microhaplotype_table_to_json_file, "Convert microhaplotype table to JSON Meta"),
                "terra_amp_output_to_json": pmofunction(terra_amp_output_to_json,
                                                                 "Convert terra output table to JSON seq table"),


            },
            "extractors_from_pmo" : {
                "extract_pmo_with_selected_meta": pmofunction(extract_pmo_with_selected_meta, "Extract from PMO samples and associated haplotypes with selected meta"),
                "extract_pmo_with_select_specimens": pmofunction(extract_pmo_with_select_specimens,"Extract from PMO specific samples"),
                "extract_pmo_with_select_targets" : pmofunction(extract_pmo_with_select_targets,"Extract from PMO specific targets"),
                "extract_allele_table": pmofunction(extract_for_allele_table,"Extract allele tables which can be as used as input to dcifer or moire"),
            },
            "working_with_multiple_pmos" : {
                "combine_pmos": pmofunction(combine_pmos,
                                                    "Combine multiple pmos of the same panel into a single pmo"),
            },
            "extract_basic_info_from_pmo" : {
                "list_specimen_meta_fields": pmofunction(list_specimen_meta_fields,
                                                    "List out the specimen meta fields in the specimen_info section"),
                "list_tar_amp_bioinformatics_info_ids": pmofunction(list_tar_amp_bioinformatics_info_ids,
                                                         "List out all the tar_amp_bioinformatics_info_ids in a PMO file"),
                "count_specimen_meta": pmofunction(count_specimen_meta,
                                                         "Count the values of specific specimen meta fields in the specimen_info section"),
                "count_targets_per_sample": pmofunction(count_targets_per_sample,
                                                   "Count the number of targets per sample"),
                "count_samples_per_target": pmofunction(count_samples_per_target,
                                                        "Count the number of samples per target"),

            }
        }
        self.version = "1.0.0"

    def printAvailableFunctions(self):
        sys.stdout.write("pmotools v" + self.version + " - A suite of tools for interacting with " + CT.boldGreen("Portable Microhaplotype Object (pmo)") + " file format" + "\n\n")

        sys.stdout.write("Available functions are" + "\n")
        for functionClass in self.functions:
            sys.stdout.write(CT.boldBlue(functionClass) + "\n")
            for functionName in self.functions[functionClass]:
                sys.stdout.write("\t" + functionName + " - " + self.functions[functionClass][functionName].shortDescription + "\n")

    def hasFunction(self, funcName):
        hasFunction = False
        for functionClass in self.functions:
            for functionName in self.functions[functionClass]:
                if functionName == funcName:
                    hasFunction = True
        return hasFunction

    def getFunction(self, funcName):
        for functionClass in self.functions:
            for functionName in self.functions[functionClass]:
                if functionName == funcName:
                    return self.functions[functionClass][functionName]





if __name__ == "__main__":
    runner = pmotools_runner()
    if len(sys.argv) == 1:
        runner.printAvailableFunctions()
        sys.exit(0)
    if not runner.hasFunction(sys.argv[1]):
        sys.stdout.write(CT.boldRed("Did not find Function ") + CT.boldWhite(sys.argv[1]) + "\n")
        runner.printAvailableFunctions()
        sys.exit(0)
    funcName = sys.argv[1]
    sys.argv.remove(funcName)
    sys.argv[0] = sys.argv[0] + " " + funcName
    runner.getFunction(funcName).func()

