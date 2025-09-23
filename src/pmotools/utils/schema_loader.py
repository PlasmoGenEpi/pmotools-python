import json
import importlib.resources as resources


def load_schema(name: str) -> dict:
    """
    Load a JSON schema from the pmotools.schemas package.

    Parameters
    ----------
    name : str
        The filename of the schema (e.g. "pmo_schema.json").

    Returns
    -------
    dict
        Parsed JSON schema as a Python dictionary.

    Raises
    ------
    FileNotFoundError
        If the schema file does not exist.
    json.JSONDecodeError
        If the schema file is not valid JSON.
    """
    with resources.files("pmotools.schemas").joinpath(name).open(
        "r", encoding="utf-8"
    ) as f:
        return json.load(f)
