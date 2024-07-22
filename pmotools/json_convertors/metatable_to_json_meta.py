#!/usr/bin/env python3
import pandas as pd
import json

def pandas_table_to_json(contents : pd.DataFrame, index_col : str | None = None ):
    """
    Convert a pandas dataframe table into a json dictionary, if there is an index column create a dictionary with the keys being the index

    :param contents: the dataframe to be converted
    :param index_col: an index column if any
    :return: a dictionary of the input table data
    """

    # Custom object_hook to replace None with an empty string
    def custom_object_hook(d):
        return {k: ("" if v is None else v) for k, v in d.items()}
    if index_col is not None:
        contents_json = json.loads(contents.to_json(orient="index",
                                                   index= True,
                                                   date_format="iso"),
                                   object_hook=custom_object_hook)
    else:
        contents_json = json.loads(contents.to_json(orient="records", date_format = "iso"), object_hook=custom_object_hook)
    return contents_json
