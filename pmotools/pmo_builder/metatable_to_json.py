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
        experiment_sample_name_col: str = 'experiment_sample_name',
        # TODO: This gets converted to ID on merge
        sequencing_info_name_col: str = 'sequencing_info_name',
        # TODO: This gets converted to ID on merge
        specimen_name_col: str = 'specimen_name',
        # TODO: This gets converted to ID on merge
        panel_name_col: str = 'panel_name',
        accession_col: str = None,
        extraction_plate_name_col: str = None,
        extraction_plate_col_col: int = None,
        extraction_plate_row_col: str = None,
        sequencing_plate_name_col: str = None,
        sequencing_plate_col_col: int = None,
        sequencing_plate_row_col: str = None,
        additional_experiment_info_cols: list | None = None,
):
    """
    Converts a DataFrame containing experiment information into JSON.

    :param contents (pd.DataFrame): Input DataFrame containing experiment data.
    :param experiment_sample_name_col (str): Column name for experiment sample names. Default: experiment_sample_name
    :param sequencing_info_name_col (str): Column name for sequencing information IDs. Default: sequencing_info_name
    :param specimen_name_col (str): Column name for specimen IDs. Default: specimen_name
    :param panel_name_col (str): Column name for panel IDs. Default: panel_name
    :param accession_col (Optional[str]): Column name for accession information.
    :param extraction_plate_name_col (Optional[str]): Column name containing plate name for extraction.
    :param extraction_plate_col_col (Optional[int]): Column name for col of sample on extraction plate.
    :param extraction_plate_row_col (Optional[str]): Column name for row of sample on extraction plate.
    :param sequencing_plate_name_col (Optional[str]): Column name containing plate name for sequencing.
    :param sequencing_plate_col_col (Optional[int]): Column name for col of sample on sequencing plate.
    :param sequencing_plate_row_col (Optional[str]): Column name for row of sample on sequencing plate.
    :param additional_experiment_info_cols (Optional[List[str], None]]): Additional column names to include.

    :return: JSON format where keys are `experiment_sample_id` and values are corresponding row data.
    """
    copy_contents = contents.copy()
    selected_columns = [
        experiment_sample_name_col,
        sequencing_info_name_col,
        specimen_name_col,
        panel_name_col
    ]

    # Add optional columns
    optional_columns = [accession_col, extraction_plate_name_col, extraction_plate_col_col,
                        extraction_plate_row_col, sequencing_plate_name_col, sequencing_plate_col_col, sequencing_plate_row_col]
    selected_columns += [col for col in optional_columns if col]

    # Include additional user-defined columns if provided
    if additional_experiment_info_cols:
        selected_columns += additional_experiment_info_cols
    # TODO: Think about what to do with the plate information is stored. The current way would need to make a table for each sample and for each type of plate. Should be flexible
    # Subset to columns
    copy_contents = copy_contents[selected_columns]

    # Convert to JSON
    copy_contents.set_index(experiment_sample_name_col,
                            drop=False, inplace=True)
    meta_json = pandas_table_to_json(copy_contents, return_indexed_dict=True)
    return meta_json


def specimen_info_table_to_json(
        contents: pd.DataFrame,
        specimen_name_col: str = 'specimen_name',
        specimen_taxon_id_col: int = 'specimen_taxon_id',
        host_taxon_id_col: int = 'host_taxon_id',
        collection_date_col: str = 'collection_date',
        collection_country_col: str = 'collection_country',
        project_name_col: str = 'project_name',
        alternate_identifiers_col: str = None,
        collector_chief_scientist_col: str = None,
        drug_usage_col: str = None,
        env_broad_scale_col: str = None,
        env_local_scale_col: str = None,
        env_medium_col: str = None,
        geo_admin1_col: str = None,
        geo_admin2_col: str = None,
        geo_admin3_col: str = None,
        host_age_col: float = None,
        host_sex_col: float = None,
        host_subject_id: str = None,
        lat_lon_col: str = None,
        parasite_density_col: float = None,
        parasite_density_method_col: str = None,
        plate_col_col: int = None,
        plate_name_col: str = None,
        plate_row_col: str = None,
        specimen_collect_device_col: str = None,
        specimen_comments_col: str = None,
        specimen_store_loc_col: str = None,
        additional_specimen_cols: list | None = None,
):
    """
    Converts a DataFrame containing specimen information into JSON.

    :param contents (pd.DataFrame): The input DataFrame containing experiment data.
    :param specimen_name_col (str): The column name for specimen sample IDs. Default: specimen_id
    :param specimen_taxon_id_col (int): NCBI taxonomy number of the organism. Default: samp_taxon_id
    :param host_taxon_id_col (int): NCBI taxonomy number of the host. Default: host_taxon_id
    :param collection_date_col (string): Date of the sample collection. Default: collection_date
    :param collection_country_col (string): Name of country collected in (admin level 0). Default : collection_country
    :param project_name_col (string): Name of the project. Default : project_name
    :param alternate_identifiers_col (Optional[str]): List of optional alternative names for the samples
    :param collector_chief_scientist_col (Optional[str]): Name of the primary person managing the specimen. Default: collector
    :param drug_usage_col (Optional[str]): Any drug used by subject and the frequency of usage; can include multiple drugs used
    :param env_broad_scale_col (Optional[str]): The broad environment from which the specimen was collected
    :param env_local_scale_col (Optional[str]): The local environment from which the specimen was collected
    :param env_medium_col (Optional[str]): The environment medium from which the specimen was collected from
    :param geo_admin1_col (Optional[str]): Geographical admin level 1
    :param geo_admin2_col (Optional[str]): Geographical admin level 2
    :param geo_admin3_col (Optional[str]): Geographical admin level 3
    :param host_age_col (Optional[float]): The age in years of the person
    :param host_sex_col (Optional[str]): If specimen is from a person, the sex of that person
    :param host_subject_id (Optional[str]): ID for the individual a specimen was collected from
    :param lat_lon_col (Optional[str]): Latitude and longitude of the collection site
    :param parasite_density_col (Optional[float]): The parasite density in parasites per microliters
    :param parasite_density_method_col (Optional[str]): The method of how this density was obtained
    :param plate_col_col (Optional[int]): Column the specimen was in in the plate
    :param plate_name_col (Optional[str]): Name of plate the specimen was in
    :param plate_row_col (Optional[str]): Row the specimen was in in the plate
    :param specimen_collect_device_col (Optional[str]): The way the specimen was collected
    :param specimen_comments_col (Optional[str]): Additional comments about the specimen
    :param specimen_store_loc_col (Optional[str]): Specimen storage site
    :param additional_specimen_cols (Optional[List[str], None]]): Additional column names to include

    :return: JSON format where keys are `specimen_name_col` and values are corresponding row data.
    """
    copy_contents = contents.copy()
    selected_columns = [
        specimen_name_col,
        specimen_taxon_id_col,
        host_taxon_id_col,
        collection_date_col,
        collection_country_col,
        project_name_col,
    ]
    # Add optional columns
    optional_columns = [alternate_identifiers_col,
                        collector_chief_scientist_col,
                        drug_usage_col,
                        env_broad_scale_col,
                        env_local_scale_col,
                        env_medium_col,
                        geo_admin1_col,
                        geo_admin2_col,
                        geo_admin3_col,
                        host_age_col,
                        host_sex_col,
                        host_subject_id,
                        lat_lon_col,
                        parasite_density_col,
                        parasite_density_method_col,
                        plate_col_col,
                        plate_name_col,
                        plate_row_col,
                        specimen_collect_device_col,
                        specimen_comments_col,
                        specimen_store_loc_col]

    selected_columns += [col for col in optional_columns if col]
    # Include additional user-defined columns if provided
    if additional_specimen_cols:
        selected_columns += additional_specimen_cols

    # Subset to columns
    copy_contents = copy_contents[selected_columns]

    meta_json = pandas_table_to_json(copy_contents, return_indexed_dict=True)

    # TODO: sort out the parasitemia and plate sections
    return meta_json
