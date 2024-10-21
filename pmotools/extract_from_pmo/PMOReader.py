#!/usr/bin/env python3


import json
import gzip
import os
from pmotools.utils.PMOChecker import PMOChecker


class PMOReader:
    """
    A class for reading in PMO from files
    """
    @staticmethod
    def read_in_pmo(fnp : str | os.PathLike[str]):
        """
        Read in a PMO file, can either be compressed(.gz) or uncompressed
        :param fnp: the file name path of the PMO file to read in
        :return: a PMO like object
        """
        # check input
        checker = PMOChecker()
        if fnp.endswith(".gz"):
            with gzip.open(fnp) as f:
                pmodata = json.load(f)
        else:
            with open(fnp) as f:
                pmodata = json.load(f)

        checker.check_for_required_base_fields(pmodata)
        return pmodata

    @staticmethod
    def read_in_pmos(fnps : list[str] | list[os.PathLike[str]]):
        """
        Read in a PMO file, can either be compressed(.gz) or uncompressed
        :param fnps: the file name path of the PMO file to read in
        :return: a list of PMO like object
        """
        ret = []
        for fnp in fnps:
            ret.append(PMOReader.read_in_pmo(fnp))
        return ret

