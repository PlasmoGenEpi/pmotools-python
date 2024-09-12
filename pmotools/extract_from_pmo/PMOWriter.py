#!/usr/bin/env python3


import json
import gzip
import os
from pmotools.utils.PMOChecker import PMOChecker
from pmotools.utils.small_utils import Utils


class PMOWriter:

    @staticmethod
    def write_out_pmo(pmo, fnp : str | os.PathLike[str], overwrite : bool = False):
        """
        Write out a PMO, will write to zip file if the output fnp name ends with .gz
        :param pmo: the PMO to write
        :param fnp: the output filename path
        :param overwrite: whether or not to overwrite output file if it exists
        :return: nothing
        """
        Utils.outputfile_check(fnp, overwrite)
        if fnp.endswith('.gz'):
            with gzip.open(fnp, 'wt', encoding="utf-8") as zipfile:
                json.dump(pmo, zipfile, indent=2)
        else:
            json.dump(pmo, open(fnp, 'w'), indent=2)

    @staticmethod
    def add_pmo_extension_as_needed(output_fnp, gzip : bool = True):
        """
        Add on json or json.gz as needed to output pmo file

        :param output_fnp: the original output filename path
        :param gzip: whether or not to gzip the output file
        :return: the output filename path with the extension added if needed
        """
        pmo_extension = ".json"
        if gzip:
            pmo_extension += ".gz"
        return Utils.appendStrAsNeeded(output_fnp, pmo_extension)
