#!/usr/bin/env python3
import pandas as pd
import json


def pandas_table_to_json(contents: pd.DataFrame, return_indexed_dict: bool = False):
    """
    Convert a pandas dataframe table into a json dictionary, if there is an index column create a dictionary with the keys being the index

    :param contents: the dataframe to be converted
    :param index_col: whether to return an indexed dictionary
    :return: a dictionary of the input table data
    """

    # Custom object_hook to replace None with an empty string
    def custom_object_hook(d):
        return {k: ("" if v is None else v) for k, v in d.items()}
    if return_indexed_dict:
        contents_json = json.loads(contents.to_json(orient="index",
                                                    index=True,
                                                    date_format="iso"),
                                   object_hook=custom_object_hook)
    else:
        contents_json = json.loads(contents.to_json(
            orient="records", date_format="iso"), object_hook=custom_object_hook)
    return contents_json


def experiment_info_table_to_json(contents: pd.DataFrame, experiment_sample_id_col: str = 'experiment_sample_id'):
    copy_contents = contents.copy()
    copy_contents.set_index(experiment_sample_id_col, drop=False, inplace=True)
    meta_json = pandas_table_to_json(copy_contents, return_indexed_dict=True)
    return meta_json


def specimen_info_table_to_json(contents: pd.DataFrame, specimen_id_col: str = 'specimen_id'):
    copy_contents = contents.copy()
    copy_contents.set_index(specimen_id_col, drop=False, inplace=True)
    meta_json = pandas_table_to_json(copy_contents, return_indexed_dict=True)
    return meta_json
