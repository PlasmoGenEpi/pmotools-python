#!/usr/bin/env python3

import argparse
from dataclasses import dataclass
from typing import Callable, Dict, Tuple

from pmotools import __version__
from pmotools.utils.color_text import ColorText as CT

# convertors_to_pmo
from pmotools.scripts.convertors_to_pmo.text_meta_to_json_meta import (
    text_meta_to_json_meta,
)
from pmotools.scripts.convertors_to_pmo.excel_meta_to_json_meta import (
    excel_meta_to_json_meta,
)
from pmotools.scripts.convertors_to_pmo.microhaplotype_table_to_json_file import (
    microhaplotype_table_to_json_file,
)
from pmotools.scripts.convertors_to_pmo.terra_amp_output_to_json import (
    terra_amp_output_to_json,
)

# extractors_from_pmo
from pmotools.scripts.extractors_from_pmo.extract_pmo_with_selected_meta import (
    extract_pmo_with_selected_meta,
)
from pmotools.scripts.extractors_from_pmo.extract_pmo_with_select_specimen_names import (
    extract_pmo_with_select_specimen_names,
)
from pmotools.scripts.extractors_from_pmo.extract_pmo_with_select_library_sample_names import (
    extract_pmo_with_select_library_sample_names,
)
from pmotools.scripts.extractors_from_pmo.extract_pmo_with_select_targets import (
    extract_pmo_with_select_targets,
)
from pmotools.scripts.extractors_from_pmo.extract_pmo_with_read_filter import (
    extract_pmo_with_read_filter,
)
from pmotools.scripts.extractors_from_pmo.extract_allele_table import (
    extract_for_allele_table,
)

# pmo_utils
from pmotools.scripts.pmo_utils.combine_pmos import combine_pmos
from pmotools.scripts.pmo_utils.validate_pmo import validate_pmo

# extract_info_from_pmo
from pmotools.scripts.extract_info_from_pmo.list_library_sample_names_per_specimen_name import (
    list_library_sample_names_per_specimen_name,
)
from pmotools.scripts.extract_info_from_pmo.list_specimen_meta_fields import (
    list_specimen_meta_fields,
)
from pmotools.scripts.extract_info_from_pmo.list_bioinformatics_run_names import (
    list_bioinformatics_run_names,
)
from pmotools.scripts.extract_info_from_pmo.count_specimen_meta import (
    count_specimen_meta,
)
from pmotools.scripts.extract_info_from_pmo.count_targets_per_library_sample import (
    count_targets_per_library_sample,
)
from pmotools.scripts.extract_info_from_pmo.count_library_samples_per_target import (
    count_library_samples_per_target,
)

# panel info subset
from pmotools.scripts.extract_info_from_pmo.extract_insert_of_panels import (
    extract_insert_of_panels,
)
from pmotools.scripts.extract_info_from_pmo.extract_refseq_of_inserts_of_panels import (
    extract_refseq_of_inserts_of_panels,
)


@dataclass(frozen=True)
class PmoCommand:
    func: Callable[[], None]
    help: str


REGISTRY: Dict[str, Dict[str, PmoCommand]] = {
    "convertors_to_json": {
        "text_meta_to_json_meta": PmoCommand(
            text_meta_to_json_meta, "Convert text file meta to JSON Meta"
        ),
        "excel_meta_to_json_meta": PmoCommand(
            excel_meta_to_json_meta, "Convert Excel file meta to JSON Meta"
        ),
        "microhaplotype_table_to_json_file": PmoCommand(
            microhaplotype_table_to_json_file,
            "Convert microhaplotype table to a JSON file",
        ),
        "terra_amp_output_to_json": PmoCommand(
            terra_amp_output_to_json, "Convert Terra output to JSON sequence table"
        ),
    },
    "extractors_from_pmo": {
        "extract_pmo_with_selected_meta": PmoCommand(
            extract_pmo_with_selected_meta,
            "Extract samples + haplotypes using selected meta",
        ),
        "extract_pmo_with_select_specimen_names": PmoCommand(
            extract_pmo_with_select_specimen_names,
            "Extract specific samples from the specimens table",
        ),
        "extract_pmo_with_select_library_sample_names": PmoCommand(
            extract_pmo_with_select_library_sample_names,
            "Extract experiment sample names from experiment_info table",
        ),
        "extract_pmo_with_select_targets": PmoCommand(
            extract_pmo_with_select_targets, "Extract specific targets"
        ),
        "extract_pmo_with_read_filter": PmoCommand(
            extract_pmo_with_read_filter, "Extract with a read filter"
        ),
        "extract_allele_table": PmoCommand(
            extract_for_allele_table,
            "Extract allele tables for tools like dcifer or moire",
        ),
        "extract_insert_of_panels": PmoCommand(
            extract_insert_of_panels, "Extract inserts of panels from a PMO"
        ),
        "extract_refseq_of_inserts_of_panels": PmoCommand(
            extract_refseq_of_inserts_of_panels,
            "Extract ref_seq of panel inserts from a PMO",
        ),
    },
    "working_with_multiple_pmos": {
        "combine_pmos": PmoCommand(
            combine_pmos, "Combine multiple PMOs of the same panel"
        ),
    },
    "extract_basic_info_from_pmo": {
        "list_library_sample_names_per_specimen_name": PmoCommand(
            list_library_sample_names_per_specimen_name,
            "List experiment_sample_ids per specimen_id",
        ),
        "list_specimen_meta_fields": PmoCommand(
            list_specimen_meta_fields,
            "List specimen meta fields in the specimen_info section",
        ),
        "list_bioinformatics_run_names": PmoCommand(
            list_bioinformatics_run_names,
            "List all tar_amp_bioinformatics_info_ids in a PMO",
        ),
        "count_specimen_meta": PmoCommand(
            count_specimen_meta, "Count values of selected specimen meta fields"
        ),
        "count_targets_per_library_sample": PmoCommand(
            count_targets_per_library_sample, "Count number of targets per sample"
        ),
        "count_library_samples_per_target": PmoCommand(
            count_library_samples_per_target, "Count number of samples per target"
        ),
    },
    "validation": {
        "validate_pmo": PmoCommand(
            validate_pmo, "Validate a PMO file against a JSON Schema"
        )
    },
}


def _iter_all_commands():
    for group, commands in REGISTRY.items():
        for name, cmd in commands.items():
            yield group, name, cmd.help


def _print_catalog_plain():
    """
    Print commands in a machine-friendly, no-color format:
    '<command>\t<group>\t<help>'
    One per line; used by bash completion.
    """
    import sys

    for group, name, cmdhelp in _iter_all_commands():
        sys.stdout.write(f"{name}\t{group}\t{cmdhelp}\n")


def _print_catalog() -> None:
    """Print all groups and their commands like your previous version."""
    import sys

    sys.stdout.write(
        f"pmotools-python v{__version__} - A suite of tools for interacting with "
        + CT.boldGreen("Portable Microhaplotype Object (PMO)")
        + " file format\n\n"
    )
    sys.stdout.write("Available functions organized by groups are\n")
    for group, commands in REGISTRY.items():
        sys.stdout.write(CT.boldBlue(group) + "\n")
        for name, cmd in commands.items():
            sys.stdout.write(f"\t{name} - {cmd.help}\n")
        sys.stdout.write("\n")


def _print_group(group: str) -> int:
    """Print a single group's commands (blue header) if it exists."""
    import sys

    if group not in REGISTRY:
        sys.stdout.write(
            CT.boldRed("Did not find group ") + CT.boldWhite(group) + "\n\n"
        )
        _print_catalog()
        return 2

    sys.stdout.write(CT.boldBlue(group) + "\n")
    for name, cmd in REGISTRY[group].items():
        sys.stdout.write(f"\t{name} - {cmd.help}\n")
    sys.stdout.write("\n")
    return 0


def _print_bash_completion():
    # NOTE: this uses --list-plain to avoid ANSI color parsing and be stable.
    script = r"""# bash completion for pmotools-python
# add the below to your ~/.bash_completion

_pmotools_python_complete()
{
    local cur prev
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # 1) Completing the command name (1st arg): list all commands
    if [[ ${COMP_CWORD} -eq 1 ]]; then
        # Our CLI prints machine-friendly list via --list-plain:
        # "<command>\t<group>\t<help>"
        local lines cmds
        lines="$(${COMP_WORDS[0]} --list-plain 2>/dev/null)"
        cmds="$(printf '%s\n' "${lines}" | awk -F'\t' '{print $1}')"
        COMPREPLY=( $(compgen -W "${cmds}" -- "${cur}") )
        return 0
    fi

    # 2) Completing flags for a leaf command: scrape leaf -h
    if [[ "${cur}" == -* ]]; then
        local helps opts
        helps="$(${COMP_WORDS[0]} ${COMP_WORDS[1]} -h 2>/dev/null)"
        # Pull out flag tokens and split comma-separated forms
        opts="$(printf '%s\n' "${helps}" \
            | sed -n 's/^[[:space:]]\{0,\}\(-[-[:alnum:]][-[:alnum:]]*\)\(, *-[[:alnum:]][-[:alnum:]]*\)\{0,\}.*/\1/p' \
            | sed 's/, / /g')"
        COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
        return 0
    fi

    # 3) Otherwise, fall back to filename completion for positional args
    COMPREPLY=( $(compgen -f -- "${cur}") )
    return 0
}

complete -F _pmotools_python_complete pmotools-python
"""
    import sys

    sys.stdout.write(script)


def _build_parser() -> (
    Tuple[argparse.ArgumentParser, Dict[str, Tuple[str, PmoCommand]]]
):
    """
    Build a flat CLI:
      pmotools-python <command> [args...]
    Returns the parser and an index mapping command_name -> (group, PmoCommand)
    """
    description = (
        f"pmotools-python v{__version__} – A suite of tools for interacting with "
        f"{CT.boldGreen('Portable Microhaplotype Object (PMO)')} files"
    )
    parser = argparse.ArgumentParser(
        prog="pmotools-python",
        description=description,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--list-plain",
        action="store_true",
        help=argparse.SUPPRESS,  # keep it hidden; for completion script
    )
    parser.add_argument(
        "--bash-completion",
        action="store_true",
        help="Print bash completion script for pmotools-python",
    )
    parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {__version__}"
    )
    parser.add_argument(
        "--list",
        nargs="?",
        const="__ALL__",
        metavar="[group]",
        help="List all commands, or only those within a specific group",
    )

    subparsers = parser.add_subparsers(
        title="Commands", dest="command", metavar="<command>"
    )

    command_index: Dict[str, Tuple[str, PmoCommand]] = {}

    for group, commands in REGISTRY.items():
        for cmd_name, cmd in commands.items():
            if cmd_name in command_index:
                # Hard fail early if duplicate command names exist across groups
                raise RuntimeError(
                    f"Duplicate command name detected: '{cmd_name}'. "
                    f"Please rename one of the commands or add an alias."
                )
            sp = subparsers.add_parser(
                cmd_name,
                help=f"{cmd.help}  [{group}]",
                description=f"{cmd.help} (group: {group})",
                add_help=False,
            )
            sp.set_defaults(_handler=cmd.func, _group=group, _cmd_name=cmd_name)
            command_index[cmd_name] = (group, cmd)

    return parser, command_index


def main(argv: list[str] | None = None) -> int:
    parser, command_index = _build_parser()
    args, unknown = parser.parse_known_args(argv)

    if getattr(args, "bash_completion", False):
        _print_bash_completion()
        return 0

    if getattr(args, "list_plain", False):
        _print_catalog_plain()
        return 0

    if getattr(args, "list", None):
        group = args.list
        if group == "__ALL__":
            _print_catalog()
            return 0
        else:
            return _print_group(group)

    # No command provided: show the catalog
    if not getattr(args, "command", None):
        _print_catalog()
        return 0

    # Dispatch to the leaf and forward remaining args to its own argparse
    handler = getattr(args, "_handler", None)
    if handler is None:
        parser.error("No handler bound for this command (internal error).")

    import sys

    leaf_prog = f"pmotools-python {getattr(args, '_cmd_name', 'unknown')}"
    old_argv = sys.argv[:]
    try:
        sys.argv = [leaf_prog, *unknown]
        handler()
    finally:
        sys.argv = old_argv

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
