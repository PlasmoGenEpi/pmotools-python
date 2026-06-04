#!/usr/bin/env python3
import pandas as pd
import json
from .json_convert_utils import remove_optional_null_values, check_null_values


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
    specimen_name_col: str = "specimen_name",
    panel_name_col: str = "panel_name",
    sequencing_info_name_col: str = None,
    alternate_identifiers_col: str = None,
    experiment_accession_col: str = None,
    fastqs_loc_col: str = None,
    library_prep_plate_name_col: str = None,
    library_prep_plate_col_col: str = None,
    library_prep_plate_row_col: str = None,
    library_prep_plate_position_col: str = None,
    parasite_density_col: str = None,
    parasite_density_method_col: str = None,
    run_accession_col: str = None,
    additional_library_sample_info_cols: list | None = None,
    list_values_library_values: list | None = ["alternate_identifiers"],
    list_values_library_values_delimiter: str = ",",
):
    """
    Convert a DataFrame containing library information into JSON.

    :param contents: input DataFrame containing library data
    :type contents: pd.DataFrame
    :param library_sample_name_col: column name for library sample names. Default: library_sample_name
    :type library_sample_name_col: str
    :param specimen_name_col: column name for specimen names. Default: specimen_name
    :type specimen_name_col: str
    :param panel_name_col: column name for panel names. Default: panel_name
    :type panel_name_col: str
    :param sequencing_info_name_col: column name for sequencing information names
    :type sequencing_info_name_col: str, optional
    :param alternate_identifiers_col: column name for alternate identifiers
    :type alternate_identifiers_col: str, optional
    :param experiment_accession_col: column name for experiment accession information
    :type experiment_accession_col: str, optional
    :param fastqs_loc_col: column name for location of fastqs
    :type fastqs_loc_col: str, optional
    :param library_prep_plate_name_col: column name containing plate name for sequencing
    :type library_prep_plate_name_col: str, optional
    :param library_prep_plate_col_col: column name for the column of the sample on the sequencing plate
    :type library_prep_plate_col_col: str, optional
    :param library_prep_plate_row_col: column name for the row of the sample on the sequencing plate
    :type library_prep_plate_row_col: str, optional
    :param library_prep_plate_position_col: column name for position on the sequencing plate (e.g. A01). Can't be set if library_prep_plate_col_col and library_prep_plate_row_col are specified.
    :type library_prep_plate_position_col: str, optional
    :param parasite_density_col: the parasite density in parasites per microliter
    :type parasite_density_col: str or list of str, optional
    :param parasite_density_method_col: the method of how the density was obtained. If set, parasite_density_col must also be specified.
    :type parasite_density_method_col: str or list of str, optional
    :param run_accession_col: column name for run accession information
    :type run_accession_col: str, optional
    :param additional_library_sample_info_cols: additional column names to include
    :type additional_library_sample_info_cols: list of str, optional
    :param list_values_library_values: columns that contain values that could be a list, delimited by list_values_library_values_delimiter
    :type list_values_library_values: list of str, optional
    :param list_values_library_values_delimiter: delimiter between list_values_library_values. Default: ','
    :type list_values_library_values_delimiter: str
    :return: JSON format where keys are ``library_sample_id`` and values are corresponding row data
    :rtype: dict
    """
    # Check contents is a dataframe
    if not isinstance(contents, pd.DataFrame):
        raise ValueError("contents must be a pandas DataFrame.")

    copy_contents = contents.copy()
    column_mapping = {library_sample_name_col: "library_sample_name"}
    required_columns = list(column_mapping.keys())
    recommended_columns = [specimen_name_col, panel_name_col, sequencing_info_name_col]

    # Add optional columns
    optional_column_mapping = {
        specimen_name_col: "specimen_name",
        panel_name_col: "panel_name",
        sequencing_info_name_col: "sequencing_info_name",
        alternate_identifiers_col: "alternate_identifiers",
        experiment_accession_col: "experiment_accession",
        fastqs_loc_col: "fastqs_loc",
        run_accession_col: "run_accession",
    }
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
            alternate_identifiers_col,
            experiment_accession_col,
            fastqs_loc_col,
        ]
    )
    check_columns_exist(copy_contents, list(column_mapping.keys()))

    # Check for null values in required columns and recommended columns
    columns_to_check = required_columns
    recommended_columns_present = list(
        set(recommended_columns) & set(column_mapping.keys())
    )
    if len(recommended_columns_present) > 0:
        columns_to_check.extend(recommended_columns_present)
    check_null_values(copy_contents, columns_to_check)

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
    meta_json = add_parasite_density_info(
        parasite_density_col,
        parasite_density_method_col,
        meta_json,
        copy_contents,
        "library_sample_name",
        entry_name="parasite_density_info",
    )
    # listify columns that contain values that could be list, are delimited by the argument list_values_library_values_delimiter
    primitives = (int, float, str, bool, complex)
    for col in list_values_library_values:
        if col in copy_contents.columns:
            for lib in meta_json:
                if isinstance(lib[col], str):
                    lib[col] = lib[col].split(list_values_library_values_delimiter)
                elif isinstance(lib[col], list):
                    pass
                elif isinstance(lib[col], primitives):
                    lib[col] = [lib[col]]
                else:
                    raise ValueError(
                        f"Column '{col}' must contain either strings or lists of strings."
                    )

    meta_json = remove_optional_null_values(
        meta_json, list(optional_column_mapping.values())
    )
    return meta_json


def specimen_info_table_to_pmo(
    contents: pd.DataFrame,
    specimen_name_col: str = "specimen_name",
    specimen_taxon_id_col: str = None,
    host_taxon_id_col: str = None,
    collection_date_col: str = None,
    collection_country_col: str = None,
    project_name_col: str = None,
    alternate_identifiers_col: str = None,
    blood_meal_col: str = None,
    drug_usage_col: str = None,
    env_broad_scale_col: str = None,
    env_local_scale_col: str = None,
    env_medium_col: str = None,
    geo_admin1_col: str = None,
    geo_admin2_col: str = None,
    geo_admin3_col: str = None,
    gravid_col: str = None,
    gravidity_col: str = None,
    has_travel_out_six_month_col: str = None,
    host_age_col: str = None,
    host_sex_col: str = None,
    host_subject_id: str = None,
    lat_lon_col: str = None,
    parasite_density_col: str = None,
    parasite_density_method_col: str = None,
    specimen_accession_col: str = None,
    storage_plate_col_col: str = None,
    storage_plate_name_col: str = None,
    storage_plate_row_col: str = None,
    storage_plate_position_col: str = None,
    specimen_collect_device_col: str = None,
    specimen_comments_col: str = None,
    specimen_store_loc_col: str = None,
    specimen_type_col: str = None,
    treatment_status_col: str = None,
    additional_specimen_cols: list | None = None,
    list_values_specimen_values: list | None = [
        "alternate_identifiers",
        "drug_usage",
        "specimen_comments",
        "treatment_status",
        "specimen_taxon_id",
    ],
    list_values_specimen_values_delimiter: str = ",",
):
    """
    Convert a DataFrame containing specimen information into JSON.

    :param contents: the input DataFrame containing specimen data
    :type contents: pd.DataFrame
    :param specimen_name_col: the column name for specimen names. Default: specimen_name
    :type specimen_name_col: str
    :param specimen_taxon_id_col: NCBI taxonomy number of the organism
    :type specimen_taxon_id_col: str, optional
    :param host_taxon_id_col: NCBI taxonomy number of the host
    :type host_taxon_id_col: str, optional
    :param collection_date_col: date of the sample collection
    :type collection_date_col: str, optional
    :param collection_country_col: name of country collected in (admin level 0)
    :type collection_country_col: str, optional
    :param project_name_col: name of the project
    :type project_name_col: str, optional
    :param alternate_identifiers_col: list of optional alternative names for the samples
    :type alternate_identifiers_col: str, optional
    :param blood_meal_col: whether the host specimen has had a recent blood meal
    :type blood_meal_col: str, optional
    :param drug_usage_col: any drug used by the subject and the frequency of usage; can include multiple drugs used
    :type drug_usage_col: str, optional
    :param env_broad_scale_col: the broad environment from which the specimen was collected
    :type env_broad_scale_col: str, optional
    :param env_local_scale_col: the local environment from which the specimen was collected
    :type env_local_scale_col: str, optional
    :param env_medium_col: the environment medium from which the specimen was collected
    :type env_medium_col: str, optional
    :param geo_admin1_col: geographical admin level 1
    :type geo_admin1_col: str, optional
    :param geo_admin2_col: geographical admin level 2
    :type geo_admin2_col: str, optional
    :param geo_admin3_col: geographical admin level 3
    :type geo_admin3_col: str, optional
    :param gravid_col: whether the host specimen is pregnant
    :type gravid_col: str, optional
    :param gravidity_col: the number of previous pregnancies
    :type gravidity_col: str, optional
    :param has_travel_out_six_month_col: whether the host specimen has travelled out from the local region in the last six months
    :type has_travel_out_six_month_col: str, optional
    :param host_age_col: the age in years of the person
    :type host_age_col: str, optional
    :param host_sex_col: if the specimen is from a person, the sex of that person
    :type host_sex_col: str, optional
    :param host_subject_id: ID for the individual a specimen was collected from
    :type host_subject_id: str, optional
    :param lat_lon_col: latitude and longitude of the collection site
    :type lat_lon_col: str, optional
    :param parasite_density_col: the parasite density in parasites per microliter
    :type parasite_density_col: str or list of str, optional
    :param parasite_density_method_col: the method of how the density was obtained. If set, parasite_density_col must also be specified.
    :type parasite_density_method_col: str or list of str, optional
    :param specimen_accession_col: the accession number of the specimen
    :type specimen_accession_col: str, optional
    :param storage_plate_col_col: column the specimen was in on the plate. If set, storage_plate_row_col must also be specified.
    :type storage_plate_col_col: str, optional
    :param storage_plate_name_col: name of the plate the specimen was in
    :type storage_plate_name_col: str, optional
    :param storage_plate_row_col: row the specimen was in on the plate. If set, storage_plate_col_col must also be specified.
    :type storage_plate_row_col: str, optional
    :param storage_plate_position_col: position of the specimen on the plate (e.g. A01). Can't be set if storage_plate_col_col and storage_plate_row_col are specified.
    :type storage_plate_position_col: str, optional
    :param specimen_collect_device_col: the way the specimen was collected
    :type specimen_collect_device_col: str, optional
    :param specimen_comments_col: additional comments about the specimen
    :type specimen_comments_col: str, optional
    :param specimen_store_loc_col: specimen storage site
    :type specimen_store_loc_col: str, optional
    :param specimen_type_col: type of specimen, e.g. negative_control, positive_control, field_sample
    :type specimen_type_col: str, optional
    :param treatment_status_col: if the person has been treated with drugs, what the treatment outcome was
    :type treatment_status_col: str, optional
    :param additional_specimen_cols: additional column names to include
    :type additional_specimen_cols: list of str, optional
    :param list_values_specimen_values: columns that contain values that could be a list, delimited by list_values_specimen_values_delimiter
    :type list_values_specimen_values: list of str, optional
    :param list_values_specimen_values_delimiter: delimiter between list_values_specimen_values. Default: ','
    :type list_values_specimen_values_delimiter: str
    :return: JSON format where keys are ``specimen_name`` and values are corresponding row data
    :rtype: dict
    """
    # Check contents is a dataframe
    if not isinstance(contents, pd.DataFrame):
        raise ValueError("contents must be a pandas DataFrame.")

    copy_contents = contents.copy()

    column_mapping = {specimen_name_col: "specimen_name"}
    required_columns = list(column_mapping.keys())
    recommended_columns = [
        specimen_taxon_id_col,
        host_taxon_id_col,
        collection_date_col,
        collection_country_col,
        project_name_col,
    ]
    optional_column_mapping = {
        specimen_taxon_id_col: "specimen_taxon_id",
        host_taxon_id_col: "host_taxon_id",
        collection_date_col: "collection_date",
        collection_country_col: "collection_country",
        project_name_col: "project_name",
        alternate_identifiers_col: "alternate_identifiers",
        drug_usage_col: "drug_usage",
        blood_meal_col: "blood_meal",
        gravid_col: "gravid",
        gravidity_col: "gravidity",
        has_travel_out_six_month_col: "has_travel_out_six_month",
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
        specimen_accession_col: "specimen_accession",
        specimen_type_col: "specimen_type",
        treatment_status_col: "treatment_status",
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
            specimen_accession_col,
            specimen_type_col,
            treatment_status_col,
            storage_plate_col_col,
            storage_plate_name_col,
            storage_plate_row_col,
            storage_plate_position_col,
            specimen_collect_device_col,
            specimen_comments_col,
            specimen_store_loc_col,
            blood_meal_col,
            gravid_col,
            gravidity_col,
            has_travel_out_six_month_col,
        ]
    )
    check_columns_exist(copy_contents, list(column_mapping.keys()))

    # Check for null values in required columns and recommended columns
    columns_to_check = required_columns
    recommended_columns_present = list(
        set(recommended_columns) & set(column_mapping.keys())
    )
    if len(recommended_columns_present) > 0:
        columns_to_check.extend(recommended_columns_present)
    check_null_values(copy_contents, columns_to_check)

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

    # listify columns that contain values that could be list, are delimited by the argument list_values_specimen_values_delimiter
    primitives = (int, float, str, bool, complex)
    for col in list_values_specimen_values:
        if col in copy_contents.columns:
            for spec in meta_json:
                if isinstance(spec[col], str):
                    spec[col] = spec[col].split(list_values_specimen_values_delimiter)
                elif isinstance(spec[col], list):
                    pass
                elif isinstance(spec[col], primitives):
                    spec[col] = [spec[col]]
                else:
                    raise ValueError(
                        f"Column '{col}' must contain either strings or lists of strings."
                    )

    meta_json = remove_optional_null_values(
        meta_json, list(optional_column_mapping.values())
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
        if plate_col_val is not None and not pd.isna(plate_col_val):
            try:
                plate_col_val = int(plate_col_val)
            except (TypeError, ValueError):
                plate_col_val = plate_col_val
        plate_info = {}
        if plate_name_val:
            plate_info["plate_name"] = plate_name_val
        if plate_row_val:
            plate_info["plate_row"] = plate_row_val
        if plate_col_val is not None and not pd.isna(plate_col_val):
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
