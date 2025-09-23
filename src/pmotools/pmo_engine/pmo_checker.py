#!/usr/bin/env python3


from jsonschema import Draft7Validator, validate


class PMOChecker:
    """
    A class to house utilities to help check the formatting of read in PMO files.
    """

    def __init__(self, pmo_jsonschema: dict):
        """
        Constructor for PMOChecker with the json read from the json schema file

        for example:
        with open("portable_microhaplotype_object_v0.1.0.schema.json") as f: pmo_jsonschema_data = json.load(f)
        PMOChecker checker(pmo_jsonschema_data)
        or use loader
        PMOChecker checker(load_schema("portable_microhaplotype_object_v0.1.0.schema.json")
        """
        self.pmo_jsonschema = pmo_jsonschema
        self.pmo_validator = Draft7Validator(pmo_jsonschema)
        # below assumes the jsonschema loaded is a specific pmo jsonschema and assumes these fields exist
        # might be a challenge to validate the validating schema
        self.all_required_base_fields = self.pmo_jsonschema["required"]
        # to find the required fields of a specific class use the following:
        # self.pmo_jsonschema["$defs"]["CLASS_NAME"]["required"],
        # e.g. for SpecimenInfo class self.pmo_jsonschema["$defs"]["SpecimenInfo"]["required"]

    def get_required_fields_for_pmo_class(self, pmo_class):
        """
        Get the required fields for the pmo_class from the pmo_jsonschema
        :param pmo_class: the class to get a required fields for, will throw an exception if class is not found within the schema
        :return: the required fields for the pmo_class
        """

        if pmo_class not in self.pmo_jsonschema["$defs"]:
            raise Exception(
                f"PMO class {pmo_class} is not found in pmo_jsonschema, available fields are {', '.join(self.pmo_jsonschema['$defs'].keys())}"
            )
        return self.pmo_jsonschema["$defs"][pmo_class]["required"]

    def validate_pmo_json(self, pmo_json):
        """
        Validate the PMO json file with loaded schema
        """
        validate(pmo_json, self.pmo_jsonschema)

    def check_for_required_base_fields(self, pmo_object):
        """
        Check that all required base fields are present in a pmo object

        :param pmo_object: the pmo object to check
        :return: return void if passes, otherwise raises an exception
        """
        missing_base_fields = []
        for base_field in self.all_required_base_fields:
            if base_field not in pmo_object:
                missing_base_fields.append(base_field)
        if len(missing_base_fields) > 0:
            raise Exception(
                "Missing required base fields: {}".format(missing_base_fields)
            )
