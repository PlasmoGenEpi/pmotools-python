def check_additional_columns_exist(df, additional_column_list):
    if additional_column_list:
        missing_cols = set(additional_column_list) - set(df.columns)
        if missing_cols:
            raise ValueError(f"Missing additional columns: {missing_cols}")


def remove_optional_null_values(json_data, optional_columns):
    """
    Remove empty values from optional fields in a list of dictionaries.

    :param json_data: List of dictionaries to process
    :param optional_columns: List of optional field names to check for empty values
    :return: List of dictionaries with empty optional fields removed

    Empty values include: None, empty strings (''), empty dicts ({}), and empty lists ([])
    """
    # Convert optional_columns to a set for faster lookup
    optional_fields_set = set(optional_columns) if optional_columns else set()

    for item in json_data:
        # Collect keys to remove to avoid modifying dict while iterating
        keys_to_remove = []
        for key, value in item.items():
            if key in optional_fields_set:
                # Check if value is empty: None, empty string, empty dict, or empty list
                if value is None or value == "" or value == {} or value == []:
                    keys_to_remove.append(key)

        # Remove the empty optional fields
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
