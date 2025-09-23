#!/usr/bin/env python3
import pandas as pd
import json


def pandas_table_to_json(contents: pd.DataFrame, return_indexed_dict: bool = False):
    """
    Convert a pandas dataframe table into a json dictionary, if there is an index column create a dictionary with the keys being the index

    :param contents: the dataframe to be converted
    :param return_indexed_dict: whether to return an indexed dictionary
    :return: a dictionary of the input table data
    """

    # Custom object_hook to replace None with an empty string
    def custom_object_hook(d):
        return {k: ("" if v is None else v) for k, v in d.items()}

    if return_indexed_dict:
        contents_json = json.loads(
            contents.to_json(orient="index", index=True, date_format="iso"),
            object_hook=custom_object_hook,
        )
    else:
        contents_json = json.loads(
            contents.to_json(orient="records", date_format="iso"),
            object_hook=custom_object_hook,
        )
    return contents_json


def library_sample_info_table_to_pmo(
    contents: pd.DataFrame,
    library_sample_name_col: str = "library_sample_name",
    sequencing_info_name_col: str = "sequencing_info_name",
    specimen_name_col: str = "specimen_name",
    panel_name_col: str = "panel_name",
    accession_col: str = None,
    library_prep_plate_name_col: str = None,
    library_prep_plate_col_col: str = None,
    library_prep_plate_row_col: str = None,
    library_prep_plate_position_col: str = None,
    additional_library_sample_info_cols: list | None = None,
):
    """
    Converts a DataFrame containing library information into JSON.

    :param contents (pd.DataFrame): Input DataFrame containing library data.
    :param library_sample_name_col (str): Column name for library sample names. Default: library_sample_name
    :param sequencing_info_name_col (str): Column name for sequencing information names. Default: sequencing_info_name
    :param specimen_name_col (str): Column name for specimen IDs. Default: specimen_name
    :param panel_name_col (str): Column name for panel IDs. Default: panel_name
    :param accession_col (Optional[str]): Column name for accession information.
    :param library_prep_plate_name_col (Optional[str]): Column name containing plate name for sequencing.
    :param library_prep_plate_col_col (Optional[str]): Column name for col of sample on sequencing plate.
    :param library_prep_plate_row_col (Optional[str]): Column name for row of sample on sequencing plate.
    :param library_prep_plate_position_col (Optional[str]): Column name for position on sequencing plate (e.g. A01). Can't be set if library_prep_plate_col_col and library_prep_plate_row_col are specified.
    :param additional_library_sample_info_cols (Optional[List[str], None]]): Additional column names to include.

    :return: JSON format where keys are `library_sample_id` and values are corresponding row data.
    """
    # Check contents is a dataframe
    if not isinstance(contents, pd.DataFrame):
        raise ValueError("contents must be a pandas DataFrame.")

    copy_contents = contents.copy()
    column_mapping = {
        library_sample_name_col: "library_sample_name",
        sequencing_info_name_col: "sequencing_info_name",
        specimen_name_col: "specimen_name",
        panel_name_col: "panel_name",
    }

    # Add optional columns
    optional_column_mapping = {accession_col: "accession"}
    column_mapping.update(
        {k: v for k, v in optional_column_mapping.items() if k is not None}
    )

    # Include additional user-defined columns if provided
    if additional_library_sample_info_cols:
        for col in additional_library_sample_info_cols:
            column_mapping[col] = col

    # Checks on columns selected
    check_unique_columns(
        [
            library_sample_name_col,
            sequencing_info_name_col,
            specimen_name_col,
            panel_name_col,
            accession_col,
        ]
    )
    check_columns_exist(copy_contents, list(column_mapping.keys()))

    # Rename and subset columns
    selected_pmo_fields = list(column_mapping.values())
    copy_contents = copy_contents.rename(columns=column_mapping)
    subset_contents = copy_contents[selected_pmo_fields]

    # Convert to format
    meta_json = pandas_table_to_json(subset_contents)
    meta_json = add_plate_info(
        library_prep_plate_col_col,
        library_prep_plate_name_col,
        library_prep_plate_row_col,
        library_prep_plate_position_col,
        meta_json,
        copy_contents,
        "specimen_name",
        "library_prep_plate_info",
    )

    return meta_json


def specimen_info_table_to_pmo(
    contents: pd.DataFrame,
    specimen_name_col: str = "specimen_name",
    specimen_taxon_id_col: int = "specimen_taxon_id",
    host_taxon_id_col: str = "host_taxon_id",
    collection_date_col: str = "collection_date",
    collection_country_col: str = "collection_country",
    project_name_col: str = "project_name",
    alternate_identifiers_col: str = None,
    # collector_chief_scientist_col: str = None,
    drug_usage_col: str = None,
    env_broad_scale_col: str = None,
    env_local_scale_col: str = None,
    env_medium_col: str = None,
    geo_admin1_col: str = None,
    geo_admin2_col: str = None,
    geo_admin3_col: str = None,
    host_age_col: str = None,
    host_sex_col: str = None,
    host_subject_id: str = None,
    lat_lon_col: str = None,
    parasite_density_col: str = None,
    parasite_density_method_col: str = None,
    storage_plate_col_col: str = None,
    storage_plate_name_col: str = None,
    storage_plate_row_col: str = None,
    storage_plate_position_col: str = None,
    specimen_collect_device_col: str = None,
    specimen_comments_col: str = None,
    specimen_store_loc_col: str = None,
    additional_specimen_cols: list | None = None,
    list_values_specimen_columns: list | None = ["alternate_identifiers_col"],
    list_values_specimen_columns_delimiter: str = ",",
):
    """
    Converts a DataFrame containing specimen information into JSON.

    :param contents (pd.DataFrame): The input DataFrame containing library data.
    :param specimen_name_col (string): The column name for specimen sample IDs. Default: specimen_id
    :param specimen_taxon_id_col (string): NCBI taxonomy number of the organism. Default: samp_taxon_id
    :param host_taxon_id_col (string): NCBI taxonomy number of the host. Default: host_taxon_id
    :param collection_date_col (string): Date of the sample collection. Default: collection_date
    :param collection_country_col (string): Name of country collected in (admin level 0). Default : collection_country
    :param project_name_col (string): Name of the project. Default : project_name
    :param alternate_identifiers_col (Optional[str]): List of optional alternative names for the samples
    :param drug_usage_col (Optional[str]): Any drug used by subject and the frequency of usage; can include multiple drugs used
    :param env_broad_scale_col (Optional[str]): The broad environment from which the specimen was collected
    :param env_local_scale_col (Optional[str]): The local environment from which the specimen was collected
    :param env_medium_col (Optional[str]): The environment medium from which the specimen was collected from
    :param geo_admin1_col (Optional[str]): Geographical admin level 1
    :param geo_admin2_col (Optional[str]): Geographical admin level 2
    :param geo_admin3_col (Optional[str]): Geographical admin level 3
    :param host_age_col (Optional[str]): The age in years of the person
    :param host_sex_col (Optional[str]): If specimen is from a person, the sex of that person
    :param host_subject_id (Optional[str]): ID for the individual a specimen was collected from
    :param lat_lon_col (Optional[str]): Latitude and longitude of the collection site
    :param parasite_density_col (Optional[str, list[str]]): The parasite density in parasites per microliters
    :param parasite_density_method_col (Optional[str or list[str]]): The method of how the density was obtained. If set parasite_density_col must also be specified.
    :param storage_plate_col_col (Optional[str]): Column the specimen was in in the plate. If set storage_plate_row_col must also be specified.
    :param storage_plate_name_col (Optional[str]): Name of plate the specimen was in
    :param storage_plate_row_col (Optional[str]): Row the specimen was in in the plate. If set storage_plate_col_col must also be specified.
    :param storage_plate_position_col (Optional[str]): Position of the specimen on the plate (e.g. A01). Can't be set if storage_plate_col_col and storage_plate_row_col are specified.
    :param specimen_collect_device_col (Optional[str]): The way the specimen was collected
    :param specimen_comments_col (Optional[str]): Additional comments about the specimen
    :param specimen_store_loc_col (Optional[str]): Specimen storage site
    :param additional_specimen_cols (Optional[List[str], None]]): Additional column names to include
    :param list_values_specimen_columns (Optional[List[str], None]): columns that contain values that could be list, are delimited by the argument list_values_specimen_columns_delimiter
    :param list_values_specimen_columns_delimiter (','): delimiter between list_values_specimen_columns

    :return: JSON format where keys are `specimen_name_col` and values are corresponding row data.
    """
    # Check contents is a dataframe
    if not isinstance(contents, pd.DataFrame):
        raise ValueError("contents must be a pandas DataFrame.")

    copy_contents = contents.copy()

    column_mapping = {
        specimen_name_col: "specimen_name",
        specimen_taxon_id_col: "specimen_taxon_id",
        host_taxon_id_col: "host_taxon_id",
        collection_date_col: "collection_date",
        collection_country_col: "collection_country",
        project_name_col: "project_name",
    }

    optional_column_mapping = {
        alternate_identifiers_col: "alternate_identifiers",
        drug_usage_col: "drug_usage",
        env_broad_scale_col: "env_broad_scale",
        env_local_scale_col: "env_local_scale",
        env_medium_col: "env_medium",
        geo_admin1_col: "geo_admin1",
        geo_admin2_col: "geo_admin2",
        geo_admin3_col: "geo_admin3",
        host_age_col: "host_age",
        host_sex_col: "host_sex",
        host_subject_id: "host_subject_id",
        lat_lon_col: "lat_lon",
        specimen_collect_device_col: "specimen_collect_device",
        specimen_comments_col: "specimen_comments",
        specimen_store_loc_col: "specimen_store_loc",
    }

    column_mapping.update(
        {k: v for k, v in optional_column_mapping.items() if k is not None}
    )

    # Include additional user-defined columns if provided
    if additional_specimen_cols:
        # selected_columns += additional_specimen_cols
        for col in additional_specimen_cols:
            column_mapping[col] = col

    # Check column selection
    check_unique_columns(
        [
            specimen_name_col,
            specimen_taxon_id_col,
            host_taxon_id_col,
            collection_date_col,
            collection_country_col,
            project_name_col,
            alternate_identifiers_col,
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
            storage_plate_col_col,
            storage_plate_name_col,
            storage_plate_row_col,
            storage_plate_position_col,
            specimen_collect_device_col,
            specimen_comments_col,
            specimen_store_loc_col,
        ]
    )
    check_columns_exist(copy_contents, list(column_mapping.keys()))

    # Rename and subset columns
    selected_pmo_fields = list(column_mapping.values())
    copy_contents = copy_contents.rename(columns=column_mapping)
    subset_contents = copy_contents[selected_pmo_fields]
    meta_json = pandas_table_to_json(subset_contents)
    meta_json = add_parasite_density_info(
        parasite_density_col,
        parasite_density_method_col,
        meta_json,
        copy_contents,
        "specimen_name",
        entry_name="parasite_density_info",
    )

    meta_json = add_plate_info(
        storage_plate_col_col,
        storage_plate_name_col,
        storage_plate_row_col,
        storage_plate_position_col,
        meta_json,
        copy_contents,
        "specimen_name",
        entry_name="storage_plate_info",
    )

    for col in list_values_specimen_columns:
        if col in meta_json:
            meta_json[col] = meta_json[col].split(
                list_values_specimen_columns_delimiter
            )
    return meta_json


def check_unique_columns(columns):
    cols_to_check = [col for col in columns if col is not None]
    if len(cols_to_check) != len(set(cols_to_check)):
        raise ValueError("Selected columns must be unique.")


def check_columns_exist(df, columns):
    missing_cols = []
    df_columns = df.columns
    for col in columns:
        if col not in df_columns:
            missing_cols.append(col)
    if missing_cols:
        raise ValueError(
            f"The following columns are not in the DataFrame: {missing_cols}"
        )


def add_plate_info(
    plate_col_col,
    plate_name_col,
    plate_row_col,
    plate_position_col,
    meta_json,
    df,
    specimen_name_col,
    entry_name="plate_info",
):
    if all(
        col is None
        for col in [plate_col_col, plate_name_col, plate_row_col, plate_position_col]
    ):
        return meta_json

    # If one of col or row are set both must be
    if (plate_row_col is None) != (plate_col_col is None):
        raise ValueError("If either plate row or column is set, then both must be.")
    # Check position isn't specified in multiple ways
    if plate_position_col:
        if plate_col_col:
            raise ValueError(
                "Plate position can be specified using either row and col, or position, but not both."
            )
        else:
            plate_row_col = "plate_row"
            plate_col_col = "plate_col"

            try:
                df[plate_row_col] = (
                    df[plate_position_col].str.extract(r"(?i)^([A-H])")[0].str.upper()
                )
                df[plate_col_col] = (
                    df[plate_position_col]
                    .str.extract(r"(?i)^[A-H]0*([1-9]|1[0-2])$")[0]
                    .astype(int)
                )
            except (AttributeError, ValueError, IndexError, KeyError) as e:
                raise ValueError(
                    f"Values in '{plate_position_col}' must start with a single letter A-H/a-h followed by number 1-12."
                ) from e

    for row in meta_json:
        content_row = df[df[specimen_name_col] == row[specimen_name_col]]
        plate_name_val = content_row[plate_name_col].iloc[0] if plate_name_col else None
        plate_row_val = (
            content_row[plate_row_col].iloc[0].upper() if plate_row_col else None
        )
        plate_col_val = content_row[plate_col_col].iloc[0] if plate_col_col else None
        plate_info = {}
        if plate_name_val:
            plate_info["plate_name"] = plate_name_val
        if plate_row_val:
            plate_info["plate_row"] = plate_row_val
        if plate_col_val:
            plate_info["plate_col"] = plate_col_val

        if plate_info:
            row[entry_name] = plate_info
    return meta_json


def add_parasite_density_info(
    parasite_density_col,
    parasite_density_method_col,
    meta_json,
    df,
    specimen_name_col,
    entry_name,
):
    density_method_pairs = []
    if parasite_density_col is None and parasite_density_method_col is None:
        pass

    elif isinstance(parasite_density_col, list):
        if parasite_density_method_col is None:
            density_method_pairs = [(d_col, None) for d_col in parasite_density_col]
        elif isinstance(parasite_density_method_col, list):
            if len(parasite_density_col) != len(parasite_density_method_col):
                raise ValueError(
                    "If both parasite_density_col and parasite_density_method_col are lists, they must be the same length."
                )
            density_method_pairs = list(
                zip(parasite_density_col, parasite_density_method_col)
            )
        else:
            raise TypeError(
                "If parasite_density_col is a list, parasite_density_method_col must be a list or None."
            )

    elif isinstance(parasite_density_col, str):
        if parasite_density_method_col is None:
            density_method_pairs = [(parasite_density_col, None)]
        elif isinstance(parasite_density_method_col, str):
            density_method_pairs = [(parasite_density_col, parasite_density_method_col)]
        else:
            raise TypeError(
                "If parasite_density_col is a string, parasite_density_method_col must be a string or None."
            )

    elif parasite_density_col is None:
        if isinstance(parasite_density_method_col, list) or isinstance(
            parasite_density_method_col, str
        ):
            raise ValueError(
                "parasite_density_method_col is set but parasite_density_col is None. Cannot proceed."
            )

    else:
        raise TypeError(
            "Invalid types for parasite_density_col and parasite_density_method_col."
        )

    # Add parasite density info to meta_json
    for row in meta_json:
        content_row = df[df[specimen_name_col] == row[specimen_name_col]]
        density_infos = []
        for density_col, method_col in density_method_pairs:
            density_val = content_row[density_col].iloc[0] if density_col else None
            method_val = content_row[method_col].iloc[0] if method_col else None
            if density_val is not None:
                info = {"parasite_density": density_val}
                if method_val is not None:
                    info["parasite_density_method"] = method_val
                density_infos.append(info)
        if density_infos:
            row[entry_name] = density_infos
    return meta_json
