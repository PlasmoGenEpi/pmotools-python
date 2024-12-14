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


def experiment_info_table_to_json(
        contents: pd.DataFrame,
        experiment_sample_id_col: str = 'experiment_sample_id',
        sequencing_info_id: str = 'sequencing_info_id',
        specimen_id: str = 'specimen_id',
        panel_id: str = 'panel_id',
        accession: str = None,
        plate_col: int = None,
        plate_name: str = None,
        plate_row: str = None,
        additional_experiment_cols: list | None = None,
):
    """
    Converts a DataFrame containing experiment information into JSON.

    :param contents (pd.DataFrame): Input DataFrame containing experiment data.
    :param experiment_sample_id_col (str): Column name for experiment sample IDs. Default: experiment_sample_id
    :param sequencing_info_id (str): Column name for sequencing information IDs. Default: sequencing_info_id
    :param specimen_id (str): Column name for specimen IDs. Default: specimen_id
    :param panel_id (str): Column name for panel IDs. Default: panel_id
    :param accession (Optional[str]): Column name for accession information.
    :param plate_col (Optional[int]): Column index for plate information.
    :param plate_name (Optional[str]): Column name for plate names.
    :param plate_row (Optional[str]): Column name for plate rows.
    :param additional_experiment_cols (Optional[List[str], None]]): Additional column names to include.

    :return: JSON format where keys are `experiment_sample_id` and values are corresponding row data.
    """
    copy_contents = contents.copy()
    selected_columns = [
        experiment_sample_id_col,
        sequencing_info_id,
        specimen_id,
        panel_id
    ]

    # Add optional columns
    optional_columns = [accession, plate_col, plate_name, plate_row]
    selected_columns += [col for col in optional_columns if col]

    # Include additional user-defined columns if provided
    if additional_experiment_cols:
        selected_columns += additional_experiment_cols

    # Subset to columns
    copy_contents = copy_contents[selected_columns]

    # Convert to JSON
    copy_contents.set_index(experiment_sample_id_col, drop=False, inplace=True)
    meta_json = pandas_table_to_json(copy_contents, return_indexed_dict=True)
    return meta_json


def specimen_info_table_to_json(
        contents: pd.DataFrame,
        specimen_id_col: str = 'specimen_id',
        samp_taxon_id: int = 'samp_taxon_id',
        collection_date: str = 'collection_date',
        collection_country: str = 'collection_country',
        collector: str = 'collector',
        samp_store_loc: str = 'samp_store_loc',
        samp_collect_device: str = 'samp_collect_device',
        project_name: str = 'project_name',
        alternate_identifiers: str = None,
        geo_admin1: str = None,
        geo_admin2: str = None,
        geo_admin3: str = None,
        host_taxon_id: int = None,
        individual_id: str = None,
        lat_lon: str = None,
        parasite_density: float = None,
        plate_col: int = None,
        plate_name: str = None,
        plate_row: str = None,
        sample_comments: str = None,
        additional_specimen_cols: list | None = None,
):
    """
    Converts a DataFrame containing specimen information into JSON.

    :param contents (pd.DataFrame): The input DataFrame containing experiment data.
    :param specimen_id_col (str): The column name for specimen sample IDs. Default: specimen_id
    :param samp_taxon_id (int): NCBI taxonomy number of the organism. Default: samp_taxon_id
    :param collection_date (string): Date of the sample collection. Default: collection_date
    :param collection_country (string): Name of country collected in (admin level 0). Default : collection_country
    :param collector (string): Name of the primary person managing the specimen. Default: collector
    :param samp_store_loc (string): Sample storage site. Default: samp_store_loc
    :param samp_collect_device (string): The way the sample was collected. Default : samp_collect_device
    :param project_name (string): Name of the project. Default : project_name
    :param alternate_identifiers (Optional[str]): List of optional alternative names for the samples
    :param geo_admin1 (Optional[str]): Geographical admin level 1
    :param geo_admin2 (Optional[str]): Geographical admin level 2
    :param geo_admin3 (Optional[str]): Geographical admin level 3
    :param host_taxon_id (Optional[int]): NCBI taxonomy number of the host
    :param individual_id (Optional[str]): ID for the individual a specimen was collected from
    :param lat_lon (Optional[str]): Latitude and longitude of the collection site
    :param parasite_density (Optional[float]): The parasite density
    :param plate_col (Optional[int]): Column the specimen was in in the plate
    :param plate_name (Optional[str]): Name of plate the specimen was in
    :param plate_row (Optional[str]): Row the specimen was in in the plate
    :param sample_comments (Optional[str]): Additional comments about the sample
    :param additional_specimen_cols (Optional[List[str], None]]): Additional column names to include

    :return: JSON format where keys are `specimen_id` and values are corresponding row data.
    """
    copy_contents = contents.copy()
    selected_columns = [
        specimen_id_col,
        samp_taxon_id,
        collection_date,
        collection_country,
        collector,
        samp_store_loc,
        samp_collect_device,
        project_name
    ]

    # Add optional columns
    optional_columns = [alternate_identifiers, geo_admin1, geo_admin2, geo_admin3, host_taxon_id,
                        individual_id, lat_lon, parasite_density, plate_col, plate_name, plate_row, sample_comments]
    selected_columns += [col for col in optional_columns if col]
    # Include additional user-defined columns if provided
    if additional_specimen_cols:
        selected_columns += additional_specimen_cols

    # Subset to columns
    copy_contents = copy_contents[selected_columns]

    copy_contents.set_index(specimen_id_col, drop=False, inplace=True)
    meta_json = pandas_table_to_json(copy_contents, return_indexed_dict=True)
    return meta_json
