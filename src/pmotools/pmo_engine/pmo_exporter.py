#!/usr/bin/env python3
import copy
import json
import os
from collections import defaultdict
from typing import NamedTuple
import pandas as pd
from openpyxl.utils import get_column_letter
from dataclasses import dataclass

from pmotools.pmo_engine.pmo_checker import PMOChecker
from pmotools.pmo_engine.pmo_processor import PMOProcessor

from pmotools import __version__ as __pmotools_version__


class BedLoc(NamedTuple):
    """
    A single BED-format genomic location.

    Used when extracting target / panel insert locations out of a PMO so they
    can be written to a BED file.

    :ivar chrom: chromosome / contig name
    :ivar start: 0-based start position
    :ivar end: end position (exclusive)
    :ivar name: target name
    :ivar score: BED score column; here the insert length (``end - start``)
    :ivar strand: ``+`` or ``-``
    :ivar ref_seq: reference sequence for the insert, empty string if not loaded
    :ivar extra_info: free-text key/value annotation, e.g. genome name/version
    """

    chrom: str
    start: int
    end: int
    name: str
    score: float
    strand: str
    ref_seq: str
    extra_info: str


class PMOExporter(object):
    """
    A collection of functions to export information out of a PMO
    """

    @staticmethod
    def _is_primitive(x) -> bool:
        """
        Check whether a value is a primitive that can be written directly to a table cell.

        :param x: the value to check
        :return: True if ``x`` is a str, int, float, bool, or None
        """
        return isinstance(x, (str, int, float, bool)) or x is None

    @staticmethod
    def _is_primitive_list(x) -> bool:
        """
        Check whether a value is a list or tuple containing only primitives.

        :param x: the value to check
        :return: True if ``x`` is a list/tuple and every element is a primitive
            (see :meth:`is_primitive`)
        """
        return isinstance(x, (list, tuple)) and all(
            PMOExporter._is_primitive(i) for i in x
        )

    @staticmethod
    def _is_exportable(x) -> bool:
        """
        Check whether a value can be exported to a flat table.

        A value is exportable if it is a primitive or a list/tuple of primitives;
        complex nested objects (e.g. TravelInfo, parasite densities) are not.

        :param x: the value to check
        :return: True if ``x`` is a primitive or a primitive list
        """
        return PMOExporter._is_primitive(x) or PMOExporter._is_primitive_list(x)

    @staticmethod
    def export_specimen_travel_meta_table(
        pmodata, separator: str = ","
    ) -> pd.DataFrame:
        """
        Export the specimen meta information of a PMO to a dataframe
        Currently avoiding exporting values of complex object types like TravelInfo or Parasite densities, best to export such values in their own tables

        :param pmodata: the pmo export the information from
        :param separator: the separator to use for list values
        :return: a pandas dataframe of the specimen metadata
        """
        rows = []
        for specimen in pmodata["specimen_info"]:
            if "travel_out_six_month" in specimen:
                for travel_meta in specimen["travel_out_six_month"]:
                    export_row = {"specimen_name": specimen["specimen_name"]}
                    for key, value in travel_meta.items():
                        if PMOExporter._is_primitive(value):
                            export_row[key] = value
                        elif PMOExporter._is_primitive_list(value):
                            export_row[key] = separator.join(str(v) for v in value)
                    rows.append(export_row)
        return pd.DataFrame(rows)

    @staticmethod
    def export_specimen_meta_table(pmodata, separator: str = ",") -> pd.DataFrame:
        """
        Export the specimen meta information of a PMO to a dataframe
        Currently avoiding exporting values of complex object types like TravelInfo or Parasite densities, best to export such values in their own tables

        :param pmodata: the pmo export the information from
        :param separator: the separator to use for list values
        :return: a pandas dataframe of the specimen metadata
        """
        rows = []
        for specimen in pmodata["specimen_info"]:
            export_row = {}
            for key, value in specimen.items():
                if "project_id" == key and "project_info" in pmodata:
                    export_row["project_name"] = pmodata["project_info"][value][
                        "project_name"
                    ]
                elif "project_id" == key:
                    # project_info is optional as of schema v1.1.0; keep the raw id
                    # rather than failing to resolve a name that isn't available
                    export_row[key] = value
                elif PMOExporter._is_primitive(value):
                    export_row[key] = value
                elif PMOExporter._is_primitive_list(value):
                    export_row[key] = separator.join(str(v) for v in value)
            rows.append(export_row)
        return pd.DataFrame(rows)

    @staticmethod
    def export_library_sample_meta_table(pmodata, separator: str = ",") -> pd.DataFrame:
        """
        Export the library_sample meta information of a PMO to a dataframe

        :param pmodata: the pmo export the information from
        :param separator: the separator to use for list values
        :return: a pandas dataframe of the library_sample metadata
        """
        rows = []
        for library_sample in pmodata["library_sample_info"]:
            export_row = {}
            for key, value in library_sample.items():
                if "sequencing_info_id" == key and "sequencing_info" in pmodata:
                    export_row["sequencing_info_name"] = pmodata["sequencing_info"][
                        value
                    ]["sequencing_info_name"]
                elif "sequencing_info_id" == key:
                    # sequencing_info is optional as of schema v1.1.0; keep the raw id
                    # rather than failing to resolve a name that isn't available
                    export_row[key] = value
                elif "specimen_id" == key:
                    export_row["specimen_name"] = pmodata["specimen_info"][value][
                        "specimen_name"
                    ]
                elif "panel_id" == key:
                    export_row["panel_name"] = pmodata["panel_info"][value][
                        "panel_name"
                    ]
                elif PMOExporter._is_primitive(value):
                    export_row[key] = value
                elif PMOExporter._is_primitive_list(value):
                    export_row[key] = separator.join(str(v) for v in value)
            rows.append(export_row)
        return pd.DataFrame(rows)

    @staticmethod
    def export_bioinformatics_run_info_meta_table(
        pmodata, separator: str = ","
    ) -> pd.DataFrame:
        """
        Export the bioinformatics_run_info meta information of a PMO to a dataframe

        :param pmodata: the pmo export the information from
        :param separator: the separator to use for list values
        :return: a pandas dataframe of the library_sample metadata
        """
        rows = []
        if "bioinformatics_run_info" not in pmodata.keys():
            raise ValueError("no bioinformatics_run_info found in input PMO")
        run_id = 0
        for bioinformatics_run_info in pmodata["bioinformatics_run_info"]:
            export_row = {}
            export_row["run_id"] = run_id
            run_id += 1
            for key, value in bioinformatics_run_info.items():
                if PMOExporter._is_primitive(value):
                    export_row[key] = value
                elif PMOExporter._is_primitive_list(value):
                    export_row[key] = separator.join(str(v) for v in value)
            rows.append(export_row)
        return pd.DataFrame(rows)

    @staticmethod
    def export_bioinformatics_methods_info_meta_table(
        pmodata, separator: str = ","
    ) -> pd.DataFrame:
        """
        Export the bioinformatics_methods_info meta information of a PMO to a dataframe

        :param pmodata: the pmo export the information from
        :param separator: the separator to use for list values
        :return: a pandas dataframe of the library_sample metadata
        """
        # @todo write pytest export_bioinformatics_methods_info_meta_table
        rows = []
        if "bioinformatics_methods_info" not in pmodata.keys():
            raise ValueError("no bioinformatics_methods_info found in input PMO")
        bioinformatics_methods_id = 0
        for bioinformatics_methods_info in pmodata["bioinformatics_methods_info"]:
            bioinformatics_methods_id += 1
            export_row = {}
            for key, value in bioinformatics_methods_info.items():
                export_row["bioinformatics_methods_id"] = bioinformatics_methods_id
                if PMOExporter._is_primitive(value):
                    export_row[key] = value
                elif PMOExporter._is_primitive_list(value):
                    export_row[key] = separator.join(str(v) for v in value)

            method_count = 0
            for method in bioinformatics_methods_info["methods"]:
                method_count += 1
                method_export_row = copy.deepcopy(export_row)
                method_export_row["method_id"] = method_count
                for method_key in method.keys():
                    method_export_row[method_key] = method[method_key]
                rows.append(method_export_row)
        return pd.DataFrame(rows)

    @staticmethod
    def export_sequencing_info_meta_table(
        pmodata, separator: str = ","
    ) -> pd.DataFrame:
        """
        Export the sequencing_info meta information of a PMO to a dataframe

        :param pmodata: the pmo export the information from
        :param separator: the separator to use for list values
        :return: a pandas dataframe of the sequencing_info metadata
        """
        # check to see sequencing_info is loaded
        if "sequencing_info" not in pmodata.keys():
            raise ValueError("no sequencing_info found in input PMO")
        rows = []
        for sequencing_info in pmodata["sequencing_info"]:
            export_row = {}
            for key, value in sequencing_info.items():
                if PMOExporter._is_primitive(value):
                    export_row[key] = value
                elif PMOExporter._is_primitive_list(value):
                    export_row[key] = separator.join(str(v) for v in value)
            rows.append(export_row)
        return pd.DataFrame(rows)

    @staticmethod
    def export_project_info_meta_table(pmodata, separator: str = ",") -> pd.DataFrame:
        """
        Export the project_info meta information of a PMO to a dataframe

        :param pmodata: the pmo export the information from
        :param separator: the separator to use for list values
        :return: a pandas dataframe of the project_info metadata
        """
        # check to see sequencing_info is loaded
        if "project_info" not in pmodata.keys():
            raise ValueError("no project_info found in input PMO")
        rows = []
        for project_info in pmodata["project_info"]:
            export_row = {}
            for key, value in project_info.items():
                if PMOExporter._is_primitive(value):
                    export_row[key] = value
                elif PMOExporter._is_primitive_list(value):
                    export_row[key] = separator.join(str(v) for v in value)
            rows.append(export_row)
        return pd.DataFrame(rows)

    @staticmethod
    def export_panel_info_meta_table(pmodata, separator: str = ",") -> pd.DataFrame:
        """
        Export the panel meta information of a PMO to a dataframe

        :param pmodata: the pmo export the information from
        :param separator: the separator to use for list values
        :return: a pandas dataframe of the panel metadata
        """
        rows = []
        for panel_info in pmodata["panel_info"]:
            export_row = {}
            for key, value in panel_info.items():
                if PMOExporter._is_primitive(value):
                    export_row[key] = value
                elif PMOExporter._is_primitive_list(value):
                    export_row[key] = separator.join(str(v) for v in value)
            reactions_for_target = defaultdict(list)
            for reaction in panel_info["reactions"]:
                for target_id in reaction["panel_targets"]:
                    reactions_for_target[
                        pmodata["target_info"][target_id]["target_name"]
                    ].append(reaction["reaction_name"])
            for target, reactions in reactions_for_target.items():
                export_row_per_target = copy.deepcopy(export_row)
                export_row_per_target["target_name"] = target
                export_row_per_target["reaction_name"] = separator.join(reactions)
                rows.append(export_row_per_target)
        return pd.DataFrame(rows)

    @staticmethod
    def export_target_info_meta_table(pmodata, separator: str = ",") -> pd.DataFrame:
        """
        Export the target meta information of a PMO to a dataframe

        :param pmodata: the pmo export the information from
        :param separator: the separator to use for list values
        :return: a pandas dataframe of the panel metadata
        """
        rows = []
        for target_info in pmodata["target_info"]:
            export_row = {}
            for key, value in target_info.items():
                if "forward_primer" == key:
                    export_row["forward_primer_seq"] = target_info["forward_primer"][
                        "seq"
                    ]
                    if "location" in target_info["forward_primer"]:
                        for primer_loc_key in target_info["forward_primer"][
                            "location"
                        ].keys():
                            export_row[
                                "forward_primer_" + primer_loc_key
                            ] = target_info["forward_primer"]["location"][
                                primer_loc_key
                            ]
                elif "reverse_primer" == key:
                    export_row["reverse_primer_seq"] = target_info["reverse_primer"][
                        "seq"
                    ]
                    if "location" in target_info["reverse_primer"]:
                        for primer_loc_key in target_info["reverse_primer"][
                            "location"
                        ].keys():
                            export_row[
                                "reverse_primer_" + primer_loc_key
                            ] = target_info["reverse_primer"]["location"][
                                primer_loc_key
                            ]
                elif "insert_location" == key:
                    for insert_key in value.keys():
                        export_row["insert_" + insert_key] = value[insert_key]
                elif PMOExporter._is_primitive(value):
                    export_row[key] = value
                elif PMOExporter._is_primitive_list(value):
                    export_row[key] = separator.join(str(v) for v in value)
            rows.append(export_row)

        df = pd.DataFrame(rows)

        priority_cols = ["target_name", "forward_primer_seq", "reverse_primer_seq"]
        leading = [c for c in priority_cols if c in df.columns]
        rest = sorted(c for c in df.columns if c not in priority_cols)

        return df[leading + rest]

    @staticmethod
    def export_pmo_header_table(pmodata, separator: str = ",") -> pd.DataFrame:
        """
        Export the pmo header meta information of a PMO to a dataframe

        :param pmodata: the pmo export the information from
        :param separator: the separator to use for list values
        :return: a pandas dataframe of the genomes metadata
        """
        rows = []

        if "pmo_header" not in pmodata.keys():
            raise ValueError("no pmo_header found in input PMO")
        export_row = {}
        for key, value in pmodata["pmo_header"].items():
            if "generation_method" == key:
                export_row["generation_method.program_version"] = value[
                    "program_version"
                ]
                export_row["generation_method.program_name"] = value["program_name"]
            elif PMOExporter._is_primitive(value):
                export_row[key] = value
            elif PMOExporter._is_primitive_list(value):
                export_row[key] = separator.join(str(v) for v in value)
        rows.append(export_row)
        df = pd.DataFrame(rows)
        priority_cols = ["pmo_version"]
        leading = [c for c in priority_cols if c in df.columns]
        rest = sorted(c for c in df.columns if c not in priority_cols)

        return df[leading + rest]

    @staticmethod
    def export_targeted_genomes_meta_table(
        pmodata, separator: str = ","
    ) -> pd.DataFrame:
        """
        Export the targeted genomes meta information of a PMO to a dataframe

        :param pmodata: the pmo export the information from
        :param separator: the separator to use for list values
        :return: a pandas dataframe of the genomes metadata
        """
        rows = []
        genome_id = 0
        if "targeted_genomes" not in pmodata.keys():
            raise ValueError("no targeted_genomes found in input PMO")
        for genome_info in pmodata["targeted_genomes"]:
            export_row = {}
            export_row["genome_id"] = genome_id
            genome_id += 1
            for key, value in genome_info.items():
                if PMOExporter._is_primitive(value):
                    export_row[key] = value
                elif PMOExporter._is_primitive_list(value):
                    export_row[key] = separator.join(str(v) for v in value)
            rows.append(export_row)

        df = pd.DataFrame(rows)
        priority_cols = ["name", "genome_version", "taxon_id", "genome_id", "url"]
        leading = [c for c in priority_cols if c in df.columns]
        rest = sorted(c for c in df.columns if c not in priority_cols)

        return df[leading + rest]

    @staticmethod
    def write_bed_locs(bed_locs: list[BedLoc], fnp, add_header: bool = False):
        """
        Write out a list of BedLoc to a file, will auto overwrite it

        :param bed_locs: a list of BedLoc
        :param fnp: output file path, will be overwritten if it exists
        :param add_header: add header of #chrom,start end,name,score,strand,ref_seq,extra_info, starts with comment so tools will treat it as a comment line
        """
        with open(fnp, "w") as f:
            if add_header:
                f.write(
                    "\t".join(
                        [
                            "#chrom",
                            "start",
                            "end",
                            "name",
                            "score",
                            "strand",
                            "ref_seq",
                            "extra_info",
                        ]
                    )
                )
                f.write("\n")
            for bed_loc in bed_locs:
                f.write(
                    "\t".join(
                        [
                            bed_loc.chrom,
                            str(bed_loc.start),
                            str(bed_loc.end),
                            bed_loc.name,
                            str(bed_loc.score),
                            bed_loc.strand,
                            str(bed_loc.ref_seq),
                            bed_loc.extra_info,
                        ]
                    )
                )
                f.write("\n")

    @staticmethod
    def extract_targets_insert_bed_loc(
        pmodata, select_target_ids: list[int] = None, sort_output: bool = True
    ):
        """
        Extract out of a PMO the insert location for targets, will add ref seq if loaded into PMO

        :param pmodata: the PMO to extract from
        :param select_target_ids: a list of target ids to select, if None will select all targets
        :param sort_output: whether to sort output by genomic location
        :return: a list of target inserts, with named tuples with fields: chrom, start, end, name, score, strand, ref_seq, extra_info
        """
        bed_loc_out = []
        if select_target_ids is None:
            select_target_ids = list(range(len(pmodata["target_info"])))
        for target_id in select_target_ids:
            tar = pmodata["target_info"][target_id]
            if "insert_location" not in tar:
                raise Exception(
                    "no insert_location in pmodata for target id "
                    + str(target_id)
                    + " target_name "
                    + str(tar["target_name"])
                    + ", cannot extract insert_location"
                )
            genome_info = pmodata["targeted_genomes"][
                tar["insert_location"]["genome_id"]
            ]
            genome_name_version = (
                genome_info["name"] + "_" + genome_info["genome_version"]
            )
            extra_info = (
                str("[") + str("genome_name_version=") + genome_name_version + ";]"
            )
            strand = (
                "+"
                if "strand" not in tar["insert_location"]
                else tar["insert_location"]["strand"]
            )
            ref_seq = (
                ""
                if "ref_seq" not in tar["insert_location"]
                else tar["insert_location"]["ref_seq"]
            )
            bed_loc_out.append(
                BedLoc(
                    tar["insert_location"]["chrom"],
                    tar["insert_location"]["start"],
                    tar["insert_location"]["end"],
                    tar["target_name"],
                    tar["insert_location"]["end"] - tar["insert_location"]["start"],
                    strand,
                    ref_seq,
                    extra_info,
                )
            )
        if sort_output:
            return sorted(bed_loc_out, key=lambda bed: (bed.chrom, bed.start, bed.end))
        return bed_loc_out

    @staticmethod
    def extract_panels_insert_bed_loc(
        pmodata, select_panel_ids: list[int] = None, sort_output: bool = True
    ):
        """
        Extract out of a PMO the insert location for panels, will add ref seq if loaded into PMO

        :param pmodata: the PMO to extract from
        :param select_panel_ids: a list of panels ids to select, if None will select all panels
        :param sort_output: whether to sort output by genomic location
        :return: a list of target inserts, with named tuples with fields: chrom, start, end, name, score, strand, ref_seq, extra_info
        """
        bed_loc_out = []
        if select_panel_ids is None:
            select_panel_ids = list(range(len(pmodata["panel_info"])))
        for panel_id in select_panel_ids:
            for reaction_id in range(len(pmodata["panel_info"][panel_id]["reactions"])):
                for target_id in pmodata["panel_info"][panel_id]["reactions"][
                    reaction_id
                ]["panel_targets"]:
                    tar = pmodata["target_info"][target_id]
                    if "insert_location" not in tar:
                        raise Exception(
                            "no insert_location in pmodata for target id "
                            + str(target_id)
                            + " target_name "
                            + str(tar["target_name"])
                            + ", cannot extract insert_location"
                        )
                    genome_info = pmodata["targeted_genomes"][
                        tar["insert_location"]["genome_id"]
                    ]
                    genome_name_version = (
                        genome_info["name"] + "_" + genome_info["genome_version"]
                    )
                    extra_info = (
                        str("[")
                        + "genome_name_version="
                        + genome_name_version
                        + ";"
                        + "panel="
                        + pmodata["panel_info"][panel_id]["panel_name"]
                        + ";"
                        + "reaction="
                        + pmodata["panel_info"][panel_id]["reactions"][reaction_id][
                            "reaction_name"
                        ]
                        + ";"
                        + "]"
                    )
                    strand = (
                        "+"
                        if "strand" not in tar["insert_location"]
                        else tar["insert_location"]["strand"]
                    )
                    ref_seq = (
                        ""
                        if "ref_seq" not in tar["insert_location"]
                        else tar["insert_location"]["ref_seq"]
                    )
                    bed_loc_out.append(
                        BedLoc(
                            tar["insert_location"]["chrom"],
                            tar["insert_location"]["start"],
                            tar["insert_location"]["end"],
                            tar["target_name"],
                            tar["insert_location"]["end"]
                            - tar["insert_location"]["start"],
                            strand,
                            ref_seq,
                            extra_info,
                        )
                    )
        if sort_output:
            return sorted(bed_loc_out, key=lambda bed: (bed.chrom, bed.start, bed.end))
        return bed_loc_out

    @staticmethod
    def extract_alleles_per_sample_table(
        pmodata,
        additional_specimen_info_fields: list[str] = None,
        additional_library_sample_info_fields: list[str] = None,
        additional_microhap_fields: list[str] = None,
        additional_representative_info_fields: list[str] = None,
        default_base_col_names: list[str] = [
            "library_sample_name",
            "target_name",
            "seq",
        ],
        jsonschema_fnp=os.path.join(
            os.path.dirname(
                os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            ),
            "schemas/",
            f"portable_microhaplotype_object_v{__pmotools_version__}.schema.json",
        ),
        validate_pmo: bool = False,
    ) -> pd.DataFrame:
        """
        Create a pd.Dataframe of sample, target and allele. Can optionally add on any other additional fields

        :param pmodata: the data to write from
        :param additional_specimen_info_fields: any additional fields to write from the specimen_info object
        :param additional_library_sample_info_fields: any additional fields to write from the library_samples object
        :param additional_microhap_fields: any additional fields to write from the microhap object
        :param additional_representative_info_fields: any additional fields to write from the representative_microhaplotype_sequences object
        :param default_base_col_names: The default column name for the library_sample_name, target_name and seq
        :param jsonschema_fnp: path to the jsonschema schema file to validate the PMO against
        :param validate_pmo: whether to validate the PMO with a jsonschema
        :return: pandas dataframe
        """

        # check input
        if validate_pmo:
            with open(jsonschema_fnp) as f:
                checker = PMOChecker(json.load(f))
                checker.validate_pmo_json(pmodata)

        # Check to see if at least 1 sample has supplied meta field
        # samples without this meta field will have NA
        if additional_specimen_info_fields is not None:
            # Find meta fields that have at least some data
            meta_fields_with_data = {
                metafield
                for metafield in additional_specimen_info_fields
                for specimen_data in pmodata["specimen_info"]
                if metafield in specimen_data
            }

            # Determine meta fields with no samples having data
            meta_fields_with_no_samples = (
                set(additional_specimen_info_fields) - meta_fields_with_data
            )

            if meta_fields_with_no_samples:
                raise Exception(
                    f"No specimen_info have data for fields: {', '.join(meta_fields_with_no_samples)}"
                )
        # Check to see if at least 1 sample has supplied meta field
        # samples without this meta field will have NA
        if additional_library_sample_info_fields is not None:
            # Find meta fields that have at least some data
            meta_fields_with_data = {
                metafield
                for metafield in additional_library_sample_info_fields
                for library_data in pmodata["library_sample_info"]
                if metafield in library_data
            }
            # Determine meta fields with no samples having data
            meta_fields_with_no_samples = (
                set(additional_library_sample_info_fields) - meta_fields_with_data
            )

            if meta_fields_with_no_samples:
                raise Exception(
                    f"No library_sample_info have data for fields: {', '.join(meta_fields_with_no_samples)}"
                )

        # Check to see if at least 1 haplotype has this field
        # samples without this meta field will have NA
        if additional_microhap_fields is not None:
            # Find meta fields that have at least some data
            additional_microhap_fields_with_data = {
                additional_microhap_field
                for additional_microhap_field in additional_microhap_fields
                for detected_microhaplotypes in pmodata["detected_microhaplotypes"]
                for library_samples_data in detected_microhaplotypes["library_samples"]
                for target_data in library_samples_data["target_results"]
                for microhap_data in target_data["mhaps"]
                if additional_microhap_field in microhap_data
            }
            # Determine meta fields with no samples having data
            additional_microhap_fields_with_no_samples = (
                set(additional_microhap_fields) - additional_microhap_fields_with_data
            )

            if additional_microhap_fields_with_no_samples:
                raise Exception(
                    f"No detected_microhaplotypes have data for fields: {', '.join(additional_microhap_fields_with_no_samples)}"
                )
        # Check to see if at least 1 haplotype has this field
        # samples without this meta field will have NA
        if additional_representative_info_fields is not None:
            # Find meta fields that have at least some data
            # not add seq as it's being added by default, so don't output twice
            additional_microhap_fields_with_data = {
                additional_microhap_field
                for additional_microhap_field in additional_representative_info_fields
                for target_data in pmodata["representative_microhaplotypes"]["targets"]
                for microhap_data in target_data["microhaplotypes"]
                if additional_microhap_field in microhap_data
                and additional_microhap_field != "seq"
            }
            # Determine meta fields with no samples having data
            additional_microhap_fields_with_no_samples = (
                set(additional_representative_info_fields)
                - additional_microhap_fields_with_data
            )

            if additional_microhap_fields_with_no_samples:
                raise Exception(
                    f"No representative_microhaplotype_sequences have data for fields: {', '.join(additional_microhap_fields_with_no_samples)}"
                )

        if len(default_base_col_names) != 3:
            raise Exception(
                "Must have 3 default columns for allele counts, not {}".format(
                    len(default_base_col_names)
                )
            )

        rows = []
        specimen_info = pmodata["specimen_info"]
        target_info = pmodata["target_info"]
        library_sample_info = pmodata["library_sample_info"]
        detected_microhaps = pmodata["detected_microhaplotypes"]
        rep_haps = pmodata["representative_microhaplotypes"]["targets"]
        bioinformatics_run_names = None
        if "bioinformatics_run_info" in pmodata:
            bioinformatics_run_names = PMOProcessor.get_bioinformatics_run_names(
                pmodata
            )
        detected_microhaplotypes_count = 0
        for bio_run_for_detected_microhaps in detected_microhaps:
            bioinformatics_run_id = None
            if "bioinformatics_run_id" in bio_run_for_detected_microhaps:
                bioinformatics_run_id = bio_run_for_detected_microhaps[
                    "bioinformatics_run_id"
                ]
            for sample_data in bio_run_for_detected_microhaps["library_samples"]:
                library_sample_id = sample_data["library_sample_id"]
                specimen_id = library_sample_info[library_sample_id]["specimen_id"]
                library_meta = library_sample_info[library_sample_id]
                specimen_meta = specimen_info[specimen_id]
                for target_data in sample_data["target_results"]:
                    target_name = target_info[
                        rep_haps[target_data["mhaps_target_id"]]["target_id"]
                    ]["target_name"]
                    for microhap_data in target_data["mhaps"]:
                        allele_id = microhap_data["mhap_id"]
                        # print(rep_haps[target_data["mhaps_target_id"]])
                        rep_hap_meta = rep_haps[target_data["mhaps_target_id"]][
                            "microhaplotypes"
                        ][allele_id]

                        row = {
                            default_base_col_names[0]: library_meta[
                                "library_sample_name"
                            ],
                            default_base_col_names[1]: target_name,
                            default_base_col_names[2]: rep_hap_meta["seq"],
                        }
                        if (
                            bioinformatics_run_names is not None
                            and bioinformatics_run_id is not None
                        ):
                            row["bioinformatics_run_name"] = bioinformatics_run_names[
                                bioinformatics_run_id
                            ]
                        else:
                            row[
                                "bioinformatics_run_name"
                            ] = f"detected_microhaplotypes_count_idx_{detected_microhaplotypes_count}"
                        if additional_library_sample_info_fields is not None:
                            for field in additional_library_sample_info_fields:
                                row[field] = library_meta.get(field, "NA")
                        if additional_specimen_info_fields is not None:
                            for field in additional_specimen_info_fields:
                                row[field] = specimen_meta.get(field, "NA")
                        if additional_microhap_fields is not None:
                            for field in additional_microhap_fields:
                                row[field] = microhap_data.get(field, "NA")
                        if additional_representative_info_fields is not None:
                            for field in additional_representative_info_fields:
                                row[field] = rep_hap_meta.get(field, "NA")
                        rows.append(row)
            detected_microhaplotypes_count += 1
        # Build and return DataFrame
        return pd.DataFrame(rows)

    @staticmethod
    def list_library_sample_names_per_specimen_name(
        pmodata,
        select_specimen_ids: list[int] = None,
        select_specimen_names: list[str] = None,
    ) -> pd.DataFrame:
        """
        List all the library_sample_names per specimen_name

        :param pmodata: the PMO
        :param select_specimen_ids: a list of specimen_ids to select, if None, all specimen_ids are used
        :param select_specimen_names: a list of specimen_names to select, if None, all specimen_names are used
        :return: a pandas dataframe with 3 columns, specimen_id, library_sample_id, and library_sample_id_count(the number of library_sample_ids per specimen_id)
        """
        if select_specimen_ids is not None and select_specimen_names is not None:
            raise ValueError(
                "Cannot specify both select_specimen_ids and select_specimen_names"
            )
        lib_samples_per_spec = defaultdict(list[str])
        if select_specimen_names is not None:
            select_specimen_ids = PMOProcessor.get_index_of_specimen_names(
                pmodata, select_specimen_names
            )
        for lib_sample in pmodata["library_sample_info"]:
            if (
                select_specimen_ids is None
                or lib_sample["specimen_id"] in select_specimen_ids
            ):
                lib_samples_per_spec[
                    pmodata["specimen_info"][lib_sample["specimen_id"]]["specimen_name"]
                ].append(lib_sample["library_sample_name"])

        specimens_not_list = []
        for specimen in pmodata["specimen_info"]:
            if specimen["specimen_name"] not in lib_samples_per_spec:
                specimens_not_list.append(specimen["specimen_name"])

        # Prepare the data for DataFrame creation
        data = []
        for specimen_name, library_sample_names in lib_samples_per_spec.items():
            for library_sample_name in library_sample_names:
                data.append(
                    {
                        "specimen_name": specimen_name,
                        "library_sample_name": library_sample_name,
                        "library_sample_count": len(library_sample_names),
                    }
                )

        # Create the DataFrame
        df = pd.DataFrame(
            data,
            columns=["specimen_name", "library_sample_name", "library_sample_count"],
        )
        return df

    @staticmethod
    def _write_sheet(writer: pd.ExcelWriter, config: "PMOExporter.SheetConfig") -> None:
        """Write a single DataFrame to an Excel sheet and autofit its columns."""
        config.df.to_excel(writer, sheet_name=config.sheet_name, index=False)
        PMOExporter._autofit_columns(
            writer,
            config.sheet_name,
            config.df,
            specific_cols=config.specific_cols,
            max_row_check=config.max_row_check,
        )

    @staticmethod
    def _autofit_columns(
        writer: pd.ExcelWriter,
        sheet_name: str,
        df: pd.DataFrame,
        specific_cols: list[str] | None = None,
        max_row_check: int | None = None,
    ) -> None:
        """
        Auto-adjusts column widths in an Excel worksheet based on content length.

        Args:
            writer:        The active ExcelWriter instance (post df.to_excel call).
            sheet_name:    The name of the worksheet to adjust.
            df:            The DataFrame that was written to the sheet.
            specific_cols: Optional list of column names to adjust. If None, all columns are adjusted.
            max_row_check: Optional max number of rows to sample when calculating width.
                           Always includes the header regardless of this value.
                           If None, all rows are checked.
        """
        worksheet = writer.sheets[sheet_name]
        columns = specific_cols if specific_cols is not None else list(df.columns)

        for column in columns:
            col_idx = df.columns.get_loc(column) + 1
            sample = (
                df[column] if max_row_check is None else df[column].iloc[:max_row_check]
            )
            max_length = max(
                sample.astype(str).map(len).max() if len(sample) > 0 else 0,
                len(str(column)),
            )
            col_letter = get_column_letter(col_idx)
            worksheet.column_dimensions[col_letter].width = max_length + 1

    @dataclass
    class SheetConfig:
        """Configuration for writing a DataFrame to an Excel sheet."""

        sheet_name: str
        df: pd.DataFrame
        max_row_check: int | None = None
        specific_cols: list[str] | None = None

    @staticmethod
    def _build_pmo_sheet_configs(pmo) -> "list[PMOExporter.SheetConfig]":
        """
        Build the ordered list of SheetConfigs to export from a PMO object.
        Optional sheets are included only if their key is present in pmo.
        """
        sheet_conf = PMOExporter.SheetConfig

        sheets = [
            sheet_conf("PMO Header", PMOExporter.export_pmo_header_table(pmo)),
            sheet_conf(
                "Required Panel Targets", PMOExporter.export_target_info_meta_table(pmo)
            ),
            sheet_conf(
                "Required Panel Info", PMOExporter.export_panel_info_meta_table(pmo)
            ),
        ]

        if "targeted_genomes" in pmo:
            sheets.append(
                sheet_conf(
                    "Optional GenomeInfo",
                    PMOExporter.export_targeted_genomes_meta_table(pmo),
                )
            )

        sheets.append(
            sheet_conf(
                "Required Microhaplotype",
                # @todo add in the optional fields of the detected_microhaplotypes and representative_microhaplotypes
                PMOExporter.extract_alleles_per_sample_table(
                    pmo, additional_microhap_fields=["reads"]
                ),
                max_row_check=10,
            )
        )
        sheets.append(
            sheet_conf(
                "Optional Specimen Level",
                PMOExporter.export_specimen_meta_table(pmo),
                max_row_check=10,
            )
        )

        sheets.append(
            sheet_conf(
                "Optional LibrarySampleInfo",
                PMOExporter.export_library_sample_meta_table(pmo),
                max_row_check=10,
            )
        )

        if "project_info" in pmo:
            sheets.append(
                sheet_conf(
                    "Optional ProjectInfo",
                    PMOExporter.export_project_info_meta_table(pmo),
                )
            )
        if "sequencing_info" in pmo:
            sheets.append(
                sheet_conf(
                    "Optional SequencingInfo",
                    PMOExporter.export_sequencing_info_meta_table(pmo),
                )
            )
        if "bioinformatics_methods_info" in pmo:
            sheets.append(
                sheet_conf(
                    "Optional Bioinformatics Methods",
                    PMOExporter.export_bioinformatics_methods_info_meta_table(pmo),
                )
            )
        if "bioinformatics_run_info" in pmo:
            sheets.append(
                sheet_conf(
                    "Optional Bioinformatics Run",
                    PMOExporter.export_bioinformatics_run_info_meta_table(pmo),
                )
            )

        return sheets

    @staticmethod
    def export_to_excel(pmo, output_path: str) -> None:
        """
        Export a PMO object to a multi-sheet Excel file.

        Args:
            pmo:         The PMO object to export.
            output_path: The path to write the Excel file to.
        """
        with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
            for config in PMOExporter._build_pmo_sheet_configs(pmo):
                PMOExporter._write_sheet(writer, config)
