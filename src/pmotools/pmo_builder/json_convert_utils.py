import math


def check_additional_columns_exist(df, additional_column_list):
    if additional_column_list:
        missing_cols = set(additional_column_list) - set(df.columns)
        if missing_cols:
            raise ValueError(f"Missing additional columns: {missing_cols}")


def remove_optional_null_values(json_data, optional_columns):
    """
    Remove empty values from optional fields in a list of dictionaries.

    Empty values include: None, empty strings (''), empty dicts ({}), empty lists ([]), and NaN.
    """
    optional_fields_set = set(optional_columns) if optional_columns else set()

    for item in json_data:
        keys_to_remove = []
        for key, value in item.items():
            if key in optional_fields_set:
                is_empty = (
                    value is None
                    or value == ""
                    or value == {}
                    or value == []
                    or (isinstance(value, float) and math.isnan(value))
                )
                if is_empty:
                    keys_to_remove.append(key)

        for key in keys_to_remove:
            del item[key]

    return json_data


def check_null_values(df, columns):
    """
    Check for null values in a list of columns

    :param df: DataFrame to check
    :param columns: List of column names to check
    :return: None
    """
    null_columns = []
    for col in columns:
        if df[col].isna().any():
            null_columns.append(col)
    if null_columns:
        raise ValueError(f"The following columns contain null values: {null_columns}")
