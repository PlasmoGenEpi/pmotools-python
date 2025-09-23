#!/usr/bin/env python3


import json
import gzip
import os
import sys

from pmotools.utils.small_utils import Utils


class PMOWriter:
    """
    A class for writing a PMO to file
    """

    @staticmethod
    def write_out_pmo(pmo, fnp: str | os.PathLike[str], overwrite: bool = False):
        """
        Write out a PMO, will write to zip file if the output fnp name ends with .gz
        :param pmo: the PMO to write
        :param fnp: the output filename path
        :param overwrite: whether to overwrite output file if it exists
        :return: nothing
        """
        Utils.outputfile_check(fnp, overwrite)
        if fnp == "STDOUT":
            json.dump(pmo, sys.stdout, indent=2)
        elif fnp.endswith(".gz"):
            with gzip.open(fnp, "wt", encoding="utf-8") as zipfile:
                json.dump(pmo, zipfile, indent=2)
        else:
            with open(fnp, "w", encoding="utf-8") as f:
                json.dump(pmo, f, indent=2)

    @staticmethod
    def add_pmo_extension_as_needed(output_fnp, gzip: bool = True):
        """
        Add on json or json.gz as needed to output pmo file

        :param output_fnp: the original output filename path
        :param gzip: whether to gzip the output file
        :return: the output filename path with the extension added if needed
        """

        # if piping to standard out then leave alone, else append as needed for gzipped output extensions
        if output_fnp == "STDOUT":
            return output_fnp
        elif gzip:
            return Utils.appendStrAsNeededDoubleEnding(output_fnp, ".json", ".gz")
        else:
            return Utils.appendStrAsNeeded(output_fnp, ".json")
