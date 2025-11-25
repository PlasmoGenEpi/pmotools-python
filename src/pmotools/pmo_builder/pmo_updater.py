#!/usr/bin/env python3

import pandas as pd
from pmotools.pmo_engine.pmo_processor import PMOProcessor
from datetime import datetime


class PMOUpdater(object):
    @staticmethod
    def check_if_date_yyyy_mm_or_yyyy_mm_dd(date_string: str) -> bool:
        """
        Checks if a string is in YYYY-MM or YYYY-MM-DD format.
        :param date_string: the string to be checked
        """
        try:
            datetime.strptime(date_string, "%Y-%m-%d")
            return True  # Matches YYYY-MM-DD
        except ValueError:
            try:
                datetime.strptime(date_string, "%Y-%m")
                return True  # Matches YYYY-MM
            except ValueError:
                return False  # Does not match either format

    @staticmethod
    def update_specimen_meta_with_traveler_info(
        pmo,
        traveler_info: pd.DataFrame,
        specimen_name_col: str = "specimen_name",
        travel_country_col: str = "travel_country",
        travel_start_col: str = "travel_start_date",
        travel_end_col: str = "travel_end_date",
        bed_net_usage_col: str = None,
        geo_admin1_col: str = None,
        geo_admin2_col: str = None,
        geo_admin3_col: str = None,
        lat_lon_col: str = None,
        replace_current_traveler_info: bool = False,
    ):
        """
        Update a PMO's specimen's metadata with travel info
        :param pmo: the PMO to update, will directly modify this PMO
        :param traveler_info: the traveler info
        :param specimen_name_col: the specimen name column within the traveler input table
        :param travel_country_col: the column name containing the traveled to country
        :param travel_start_col: the column name containing the traveled start date, format YYYY-MM-DD or YYYY-MM
        :param travel_end_col: the column name containing the traveled end date, format YYYY-MM-DD or YYYY-MM
        :param bed_net_usage_col: (Optional) a number between 0 - 1 for rough frequency of bednet usage while traveling
        :param geo_admin1_col: (Optional) the column name containing the traveled to country admin level 1 info
        :param geo_admin2_col: (Optional) the column name containing the traveled to country admin level 2 info
        :param geo_admin3_col: (Optional) the column name containing the traveled to country admin level 3 info
        :param lat_lon_col: (Optional) the latitude and longitude column name containing the region traveled to latitude and longitude
        :param replace_current_traveler_info: whether to replace current travel info
        :return: a reference to the updated PMO
        """
        required_cols = [
            specimen_name_col,
            travel_country_col,
            travel_start_col,
            travel_end_col,
        ]
        if bed_net_usage_col is not None:
            required_cols.append(bed_net_usage_col)
        if geo_admin1_col is not None:
            required_cols.append(geo_admin1_col)
        if geo_admin2_col is not None:
            required_cols.append(geo_admin2_col)
        if geo_admin3_col is not None:
            required_cols.append(geo_admin3_col)
        if lat_lon_col is not None:
            required_cols.append(lat_lon_col)

        if not set(required_cols).issubset(traveler_info.columns):
            raise Exception(
                "missing traveler_info columns: " + ",".join(required_cols),
                " columns in table: " + ",".join(traveler_info.columns),
            )

        specimen_names_in_pmo = set(PMOProcessor.get_specimen_names(pmo))
        specimen_names_in_traveler_info = set(
            traveler_info[specimen_name_col].astype(str).tolist()
        )

        # check to see if provided traveler info for a specimen that cannot be found in PMO
        missing_traveler_specs = specimen_names_in_traveler_info - specimen_names_in_pmo

        if missing_traveler_specs:
            raise ValueError(
                f"Provided traveler info for the following specimens but they are missing from the PMO: {sorted(missing_traveler_specs)}"
            )
        traveler_info_records = traveler_info[required_cols].to_dict(orient="records")
        spec_indexs = PMOProcessor.get_index_key_of_specimen_names(pmo)

        # prep traveler info lists, clear the list if we are replacing or start an empty list to append to if none exist already
        for specimen_name in specimen_names_in_traveler_info:
            if (
                replace_current_traveler_info
                or "travel_out_six_month"
                not in pmo["specimen_info"][spec_indexs[specimen_name]]
            ):
                pmo["specimen_info"][spec_indexs[specimen_name]][
                    "travel_out_six_month"
                ] = []

        for travel_rec in traveler_info_records:
            specimen_name = str(travel_rec[specimen_name_col])
            # Validate date formats
            for date_col in (travel_start_col, travel_end_col):
                val = travel_rec[date_col]
                if pd.isna(val):
                    raise ValueError(
                        f"Missing required date value in column '{date_col}' for specimen '{specimen_name}'"
                    )
                val_str = str(val)
                if not PMOUpdater.check_if_date_yyyy_mm_or_yyyy_mm_dd(val_str):
                    raise ValueError(
                        f"Invalid date format in '{date_col}' for specimen '{specimen_name}': '{val_str}'. "
                        f"Expected YYYY-MM or YYYY-MM-DD"
                    )
            # add in travel_rec
            travel_rec.pop(specimen_name_col, None)
            pmo["specimen_info"][spec_indexs[specimen_name]][
                "travel_out_six_month"
            ].append(travel_rec)
        return pmo
