#!/usr/bin/env python3

import pandas as pd
from pmotools.pmo_engine.pmo_processor import PMOProcessor
from datetime import datetime
import copy
import logging
from typing import Any

logger = logging.getLogger(__name__)


class PMOUpdater(object):
    @staticmethod
    def _check_if_date_yyyy_mm_or_yyyy_mm_dd(date_string: str) -> bool:
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
        optional_cols = {
            "geo_admin1": geo_admin1_col,
            "geo_admin2": geo_admin2_col,
            "geo_admin3": geo_admin3_col,
            "lat_lon": lat_lon_col,
        }
        PMOUpdater._check_columns(
            traveler_info,
            [
                specimen_name_col,
                travel_country_col,
                travel_start_col,
                travel_end_col,
                bed_net_usage_col,
                *optional_cols.values(),
            ],
            "traveler_info",
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

        for _, row in traveler_info.iterrows():
            specimen_name = str(row[specimen_name_col])
            if PMOUpdater._is_blank(row[travel_country_col]):
                raise ValueError(
                    f"Missing required value in column '{travel_country_col}' for specimen '{specimen_name}'"
                )
            travel_rec = {"travel_country": str(row[travel_country_col])}
            # Validate date formats
            for date_col, key in (
                (travel_start_col, "travel_start_date"),
                (travel_end_col, "travel_end_date"),
            ):
                val = row[date_col]
                if PMOUpdater._is_blank(val):
                    raise ValueError(
                        f"Missing required date value in column '{date_col}' for specimen '{specimen_name}'"
                    )
                val_str = str(val)
                if not PMOUpdater._check_if_date_yyyy_mm_or_yyyy_mm_dd(val_str):
                    raise ValueError(
                        f"Invalid date format in '{date_col}' for specimen '{specimen_name}': '{val_str}'. "
                        f"Expected YYYY-MM or YYYY-MM-DD"
                    )
                travel_rec[key] = val_str
            if bed_net_usage_col and not PMOUpdater._is_blank(row[bed_net_usage_col]):
                travel_rec["bed_net_usage"] = float(row[bed_net_usage_col])
            for key, col in optional_cols.items():
                if col and not PMOUpdater._is_blank(row[col]):
                    travel_rec[key] = str(row[col])
            pmo["specimen_info"][spec_indexs[specimen_name]][
                "travel_out_six_month"
            ].append(travel_rec)
        return pmo

    @staticmethod
    def _is_blank(value) -> bool:
        if value is None:
            return True
        if isinstance(value, str):
            return not value.strip()
        return bool(pd.isna(value))

    @staticmethod
    def _check_columns(table: pd.DataFrame, cols: list, table_label: str) -> None:
        missing = [col for col in cols if col is not None and col not in table.columns]
        if missing:
            raise ValueError(
                f"Missing {table_label} columns: {missing}. Columns in table: {list(table.columns)}"
            )

    @staticmethod
    def _genomic_location_from_row(
        row,
        context: str,
        genome_id_col: str | None,
        chrom_col: str,
        start_col: str,
        end_col: str,
        strand_col: str | None = None,
        ref_seq_col: str | None = None,
        alt_seq_col: str | None = None,
    ) -> dict:
        """
        Build a GenomicLocation dict from a table row. If genome_id_col is None the genome_id defaults to 0.
        """
        required = [chrom_col, start_col, end_col]
        if genome_id_col is not None:
            required.append(genome_id_col)
        for col in required:
            if PMOUpdater._is_blank(row[col]):
                raise ValueError(
                    f"Missing required value in column '{col}' for {context}"
                )
        try:
            location = {
                "genome_id": int(row[genome_id_col]) if genome_id_col else 0,
                "chrom": str(row[chrom_col]),
                "start": int(row[start_col]),
                "end": int(row[end_col]),
            }
        except (TypeError, ValueError) as e:
            raise ValueError(
                f"genome id, start and end must be integers for {context}"
            ) from e
        for key, col in (
            ("strand", strand_col),
            ("ref_seq", ref_seq_col),
            ("alt_seq", alt_seq_col),
        ):
            if col and not PMOUpdater._is_blank(row[col]):
                location[key] = str(row[col])
        return location

    @staticmethod
    def _get_representative_microhaplotype(
        pmo, target_name: str, seq: str, target_indexes: dict
    ) -> dict:
        if target_name not in target_indexes:
            raise ValueError(
                f"Target '{target_name}' has no representative microhaplotypes in the PMO"
            )
        microhaplotypes = pmo["representative_microhaplotypes"]["targets"][
            target_indexes[target_name]
        ]["microhaplotypes"]
        for microhaplotype in microhaplotypes:
            if microhaplotype["seq"] == seq:
                return microhaplotype
        raise ValueError(
            f"No representative microhaplotype with seq '{seq}' for target '{target_name}'"
        )

    @staticmethod
    def update_target_info_with_markers_of_interest(
        pmo,
        markers_info: pd.DataFrame,
        target_name_col: str = "target_name",
        chrom_col: str = "chrom",
        start_col: str = "start",
        end_col: str = "end",
        genome_id_col: str | None = None,
        strand_col: str | None = None,
        ref_seq_col: str | None = None,
        alt_seq_col: str | None = None,
        associations_col: str | None = None,
        associations_delim: str = ",",
        replace_current_markers: bool = False,
    ):
        """
        Update a PMO's target_info with markers of interest, from a table with one row per marker

        :param pmo: the PMO to update, will directly modify this PMO
        :param markers_info: the markers of interest table
        :param target_name_col: the column containing the name of the target the marker is covered by
        :param chrom_col: the column containing the chromosome of the marker
        :param start_col: the column containing the start of the marker, 0-based
        :param end_col: the column containing the end of the marker, 0-based
        :param genome_id_col: (Optional) the column containing the index into targeted_genomes. Default: 0 for every marker
        :param strand_col: (Optional) the column containing the strand of the marker
        :param ref_seq_col: (Optional) the column containing the reference sequence of the marker
        :param alt_seq_col: (Optional) the column containing an alternative sequence of the marker
        :param associations_col: (Optional) the column containing associations with the marker, e.g. SP resistance
        :param associations_delim: the delimiter between associations. Default: ','
        :param replace_current_markers: whether to replace current markers of interest for the targets in the table
        :return: a reference to the updated PMO
        """
        PMOUpdater._check_columns(
            markers_info,
            [
                target_name_col,
                chrom_col,
                start_col,
                end_col,
                genome_id_col,
                strand_col,
                ref_seq_col,
                alt_seq_col,
                associations_col,
            ],
            "markers_info",
        )
        target_indexes = PMOProcessor.get_index_key_of_target_names(pmo)
        target_names_in_table = set(markers_info[target_name_col].astype(str))
        missing_targets = target_names_in_table - set(target_indexes)
        if missing_targets:
            raise ValueError(
                f"Provided markers of interest for the following targets but they are missing from the PMO: {sorted(missing_targets)}"
            )

        for target_name in target_names_in_table:
            target = pmo["target_info"][target_indexes[target_name]]
            if replace_current_markers or not target.get("markers_of_interest"):
                target["markers_of_interest"] = []

        for _, row in markers_info.iterrows():
            target_name = str(row[target_name_col])
            marker = {
                "marker_location": PMOUpdater._genomic_location_from_row(
                    row,
                    f"marker of interest for target '{target_name}'",
                    genome_id_col,
                    chrom_col,
                    start_col,
                    end_col,
                    strand_col,
                    ref_seq_col,
                    alt_seq_col,
                )
            }
            if associations_col and not PMOUpdater._is_blank(row[associations_col]):
                marker["associations"] = [
                    association.strip()
                    for association in str(row[associations_col]).split(
                        associations_delim
                    )
                    if association.strip()
                ]
            pmo["target_info"][target_indexes[target_name]][
                "markers_of_interest"
            ].append(marker)
        return pmo

    @staticmethod
    def update_representative_microhaplotypes_with_seq_variants(
        pmo,
        seq_variants_info: pd.DataFrame,
        target_name_col: str = "target_name",
        seq_col: str = "seq",
        chrom_col: str = "chrom",
        start_col: str = "start",
        end_col: str = "end",
        genome_id_col: str | None = None,
        strand_col: str | None = None,
        ref_seq_col: str | None = None,
        alt_seq_col: str | None = None,
        replace_current_seq_variants: bool = False,
    ):
        """
        Update a PMO's representative microhaplotypes with associated sequence variants, from a table with one row per variant

        :param pmo: the PMO to update, will directly modify this PMO
        :param seq_variants_info: the sequence variants table
        :param target_name_col: the column containing the target name of the microhaplotype
        :param seq_col: the column containing the sequence of the microhaplotype
        :param chrom_col: the column containing the chromosome of the variant
        :param start_col: the column containing the start of the variant, 0-based
        :param end_col: the column containing the end of the variant, 0-based
        :param genome_id_col: (Optional) the column containing the index into targeted_genomes. Default: 0 for every variant
        :param strand_col: (Optional) the column containing the strand of the variant
        :param ref_seq_col: (Optional) the column containing the reference sequence of the variant
        :param alt_seq_col: (Optional) the column containing the alternative sequence of the variant
        :param replace_current_seq_variants: whether to replace current sequence variants for the microhaplotypes in the table
        :return: a reference to the updated PMO
        """
        PMOUpdater._check_columns(
            seq_variants_info,
            [
                target_name_col,
                seq_col,
                chrom_col,
                start_col,
                end_col,
                genome_id_col,
                strand_col,
                ref_seq_col,
                alt_seq_col,
            ],
            "seq_variants_info",
        )
        target_indexes = (
            PMOProcessor.get_index_key_of_target_in_representative_microhaplotypes(pmo)
        )
        microhaplotypes = [
            PMOUpdater._get_representative_microhaplotype(
                pmo, str(row[target_name_col]), str(row[seq_col]), target_indexes
            )
            for _, row in seq_variants_info.iterrows()
        ]

        for microhaplotype in microhaplotypes:
            if replace_current_seq_variants or not microhaplotype.get(
                "associated_seq_variants"
            ):
                microhaplotype["associated_seq_variants"] = []

        for (_, row), microhaplotype in zip(
            seq_variants_info.iterrows(), microhaplotypes
        ):
            microhaplotype["associated_seq_variants"].append(
                PMOUpdater._genomic_location_from_row(
                    row,
                    f"sequence variant for target '{row[target_name_col]}'",
                    genome_id_col,
                    chrom_col,
                    start_col,
                    end_col,
                    strand_col,
                    ref_seq_col,
                    alt_seq_col,
                )
            )
        return pmo

    @staticmethod
    def update_representative_microhaplotypes_with_protein_variants(
        pmo,
        protein_variants_info: pd.DataFrame,
        target_name_col: str = "target_name",
        seq_col: str = "seq",
        transcript_col: str = "transcript",
        protein_start_col: str = "protein_start",
        protein_end_col: str = "protein_end",
        protein_genome_id_col: str | None = None,
        protein_ref_seq_col: str | None = None,
        protein_alt_seq_col: str | None = None,
        gene_name_col: str | None = None,
        alternative_gene_name_col: str | None = None,
        codon_chrom_col: str | None = None,
        codon_start_col: str | None = None,
        codon_end_col: str | None = None,
        codon_genome_id_col: str | None = None,
        codon_strand_col: str | None = None,
        codon_ref_seq_col: str | None = None,
        codon_alt_seq_col: str | None = None,
        replace_current_protein_variants: bool = False,
    ):
        """
        Update a PMO's representative microhaplotypes with associated protein variants, from a table with one row per variant

        :param pmo: the PMO to update, will directly modify this PMO
        :param protein_variants_info: the protein variants table
        :param target_name_col: the column containing the target name of the microhaplotype
        :param seq_col: the column containing the sequence of the microhaplotype
        :param transcript_col: the column containing the transcript name, used as the chrom of the protein location
        :param protein_start_col: the column containing the start of the variant within the protein, 0-based
        :param protein_end_col: the column containing the end of the variant within the protein, 0-based
        :param protein_genome_id_col: (Optional) the column containing the index into targeted_genomes for the protein location. Default: 0
        :param protein_ref_seq_col: (Optional) the column containing the reference amino acid(s)
        :param protein_alt_seq_col: (Optional) the column containing the alternative amino acid(s)
        :param gene_name_col: (Optional) the column containing the gene name
        :param alternative_gene_name_col: (Optional) the column containing an alternative gene name
        :param codon_chrom_col: (Optional) the column containing the chromosome of the codon. Must be set with codon_start_col and codon_end_col
        :param codon_start_col: (Optional) the column containing the genomic start of the codon, 0-based
        :param codon_end_col: (Optional) the column containing the genomic end of the codon, 0-based
        :param codon_genome_id_col: (Optional) the column containing the index into targeted_genomes for the codon location. Default: 0
        :param codon_strand_col: (Optional) the column containing the strand of the codon
        :param codon_ref_seq_col: (Optional) the column containing the reference sequence of the codon
        :param codon_alt_seq_col: (Optional) the column containing the alternative sequence of the codon
        :param replace_current_protein_variants: whether to replace current protein variants for the microhaplotypes in the table
        :return: a reference to the updated PMO
        """
        codon_location_cols = [codon_chrom_col, codon_start_col, codon_end_col]
        if any(codon_location_cols) and not all(codon_location_cols):
            raise ValueError(
                "If any of codon_chrom_col, codon_start_col or codon_end_col are set, all three must be."
            )
        PMOUpdater._check_columns(
            protein_variants_info,
            [
                target_name_col,
                seq_col,
                transcript_col,
                protein_start_col,
                protein_end_col,
                protein_genome_id_col,
                protein_ref_seq_col,
                protein_alt_seq_col,
                gene_name_col,
                alternative_gene_name_col,
                *codon_location_cols,
                codon_genome_id_col,
                codon_strand_col,
                codon_ref_seq_col,
                codon_alt_seq_col,
            ],
            "protein_variants_info",
        )
        target_indexes = (
            PMOProcessor.get_index_key_of_target_in_representative_microhaplotypes(pmo)
        )
        microhaplotypes = [
            PMOUpdater._get_representative_microhaplotype(
                pmo, str(row[target_name_col]), str(row[seq_col]), target_indexes
            )
            for _, row in protein_variants_info.iterrows()
        ]

        for microhaplotype in microhaplotypes:
            if replace_current_protein_variants or not microhaplotype.get(
                "associated_protein_variants"
            ):
                microhaplotype["associated_protein_variants"] = []

        for (_, row), microhaplotype in zip(
            protein_variants_info.iterrows(), microhaplotypes
        ):
            context = f"protein variant for target '{row[target_name_col]}'"
            variant = {
                "protein_location": PMOUpdater._genomic_location_from_row(
                    row,
                    context,
                    protein_genome_id_col,
                    transcript_col,
                    protein_start_col,
                    protein_end_col,
                    ref_seq_col=protein_ref_seq_col,
                    alt_seq_col=protein_alt_seq_col,
                )
            }
            for key, col in (
                ("gene_name", gene_name_col),
                ("alternative_gene_name", alternative_gene_name_col),
            ):
                if col and not PMOUpdater._is_blank(row[col]):
                    variant[key] = str(row[col])
            if codon_chrom_col and not all(
                PMOUpdater._is_blank(row[col]) for col in codon_location_cols
            ):
                variant[
                    "codon_genomic_location"
                ] = PMOUpdater._genomic_location_from_row(
                    row,
                    f"codon of {context}",
                    codon_genome_id_col,
                    codon_chrom_col,
                    codon_start_col,
                    codon_end_col,
                    codon_strand_col,
                    codon_ref_seq_col,
                    codon_alt_seq_col,
                )
            microhaplotype["associated_protein_variants"].append(variant)
        return pmo

    @staticmethod
    def merge_dicts_by_key(
        main_list: list[dict],
        update_list: list[dict],
        key_field: str,
        replace: bool = False,
        ignore_fields: list[str] | None = None,
    ) -> list[dict]:
        """
        Merge two lists of dicts by a shared key field.

        The first list is treated as the main/base data source. The second list
        provides updates that are applied on top. Both input lists are left
        untouched (deep copies are used internally).

        Args:
            main_list:     The primary list of dicts (source of truth).
            update_list:   The list of dicts whose values will be merged in.
            key_field:     The dict key used to match records across lists.
            replace:       If True, existing values in main are overwritten by
                           update values. If False, a conflict raises a ValueError.
            ignore_fields: Optional list of field names to skip entirely during
                           the merge (they are never read from update_list).

        Returns:
            A new list of dicts with updates applied.

        Raises:
            ValueError: If either list contains duplicate values for key_field.
            KeyError:   If any dict in either list is missing key_field.
            KeyError:   If update_list contains a key_field value that does not
                        exist in main_list.
            ValueError: If replace=False and an update would overwrite an
                        existing field.
        """
        ignore_fields = set(ignore_fields or [])

        # check to see if any of the input (the main or the update lists) have missing key_field
        def _check_missing_key(lst: list[dict], label: str) -> None:
            bad = [i for i, d in enumerate(lst) if key_field not in d]
            if bad:
                raise KeyError(f"{label} is missing '{key_field}' at index(es): {bad}")

        _check_missing_key(main_list, "main_list")
        _check_missing_key(update_list, "update_list")

        # check if there are duplicate key_field values
        def _check_duplicates(lst: list[dict], label: str) -> None:
            seen: set = set()
            dupes: set = set()
            for d in lst:
                val = d[key_field]
                (dupes if val in seen else seen).add(val)
            if dupes:
                raise ValueError(
                    f"{label} contains duplicate '{key_field}' values: {sorted(dupes)}"
                )

        _check_duplicates(main_list, "main_list")
        _check_duplicates(update_list, "update_list")

        # Build lookup from deep copies
        main_map: dict[Any, dict] = {d[key_field]: copy.deepcopy(d) for d in main_list}
        update_map: dict[Any, dict] = {
            d[key_field]: copy.deepcopy(d) for d in update_list
        }

        # update keys must exist in main
        extra_keys = set(update_map) - set(main_map)
        if extra_keys:
            raise KeyError(
                f"update_list contains '{key_field}' values not found in "
                f"main_list: {sorted(extra_keys)}"
            )

        # Warn if any of the main keys absent from update, this way can update some of the values but
        # not necessary to update all of them
        missing_from_update = set(main_map) - set(update_map)
        if missing_from_update:
            logger.warning(
                "The following '%s' values are in main_list but not in "
                "update_list (skipping): %s",
                key_field,
                sorted(missing_from_update),
            )

        # now merge
        for key, update_dict in update_map.items():
            main_dict = main_map[key]
            for field, value in update_dict.items():
                if field == key_field or field in ignore_fields:
                    continue
                if field in main_dict:
                    if not replace:
                        raise ValueError(
                            f"Field '{field}' already exists in record "
                            f"'{key_field}={key}' and replace=False."
                        )
                    main_dict[field] = value
                else:
                    main_dict[field] = value
        return list(main_map.values())
