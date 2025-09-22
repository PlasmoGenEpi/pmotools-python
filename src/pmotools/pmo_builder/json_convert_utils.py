def check_additional_columns_exist(df, additional_column_list):
    if additional_column_list:
        missing_cols = set(additional_column_list) - set(df.columns)
        if missing_cols:
            raise ValueError(f"Missing additional columns: {missing_cols}")
