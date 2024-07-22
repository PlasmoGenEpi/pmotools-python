#!/usr/bin/env python3



import sys, os, argparse, json
import pandas as pd
sys.path.append(os.path.join(os.path.dirname(__file__), "convertors"))
sys.path.append(os.path.join(os.path.dirname(__file__), "extractors"))

from text_meta_to_json_meta import text_meta_to_json_meta
from excel_meta_to_json_meta import excel_meta_to_json_meta
from extractor_microhapseq_with_selected_meta import extractor_microhapseq_with_selected_meta
from microhaplotype_table_to_json_file import microhaplotype_table_to_json_file


from pmotools.utils.color_text import ColorText as CT

class pmofunction:
    def __init__(self, func, shortDescription):
        self.func = func
        self.shortDescription = shortDescription

class pmotools_runner   :
    def __init__(self) :
        self.functions = {
            "convertors" : {
                "text_meta_to_json_meta": pmofunction(text_meta_to_json_meta, "Convert text file meta to JSON Meta"),
                "excel_meta_to_json_meta": pmofunction(excel_meta_to_json_meta, "Convert excel file meta to JSON Meta"),
                "microhaplotype_table_to_json_file": pmofunction(microhaplotype_table_to_json_file, "Convert microhaplotype table to JSON Meta"),
            },
            "extractors" : {
                "extractor_microhapseq_with_selected_meta": pmofunction(extractor_microhapseq_with_selected_meta, "Extract microhaplotype sequence with selected meta"),
            }
        }
        self.version = "1.0.0"

    def printAvailableFunctions(self):
        sys.stdout.write("pmotools v" + self.version + " - A suite of tools for interacting with " + CT.boldGreen("Portable Microhaplotype Object (pmo)") + " file format" + "\n\n")

        sys.stdout.write("Available functions are" + "\n")
        for functionClass in self.functions:
            sys.stdout.write(CT.boldBlue(functionClass) + "\n")
            for functionName in self.functions[functionClass]:
                sys.stdout.write("\t" + CT.boldWhite(functionName) + " - " + self.functions[functionClass][functionName].shortDescription + "\n")

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

