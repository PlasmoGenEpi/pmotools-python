#!/usr/bin/env python3
"""
Generate Sphinx .rst pages documenting every pmotools-python CLI command,
using the sphinx-argparse `.. argparse::` directive.

Walks REGISTRY from the main CLI module, finds the parser-builder function
for each command, and writes one .rst page per command group.

Run from the project root:
    python man/gen_cli_docs.py

Then make sure man/source/index.rst includes `commands/index` in a toctree.
"""

from __future__ import annotations

import importlib
from pathlib import Path

# --- adjust this import to wherever REGISTRY actually lives ------------------
from pmotools.cli import REGISTRY
# ----------------------------------------------------------------------------

# Where to write the generated pages (mirrors your man/source layout)
OUTPUT_DIR = Path(__file__).parent / "source" / "commands"

# How the tool is invoked on the command line
PROG = "pmotools-python"


def parser_func_candidates(command_name: str) -> list[str]:
    """
    Candidate names for the parser-building function inside a leaf module,
    tried in order; the first that exists wins. Lets you migrate commands
    to a uniform ``get_parser`` convention one at a time.
    """
    return ["get_parser", f"get_parser_{command_name}", "build_parser"]


def find_parser_func(module_name: str, command_name: str) -> str | None:
    """Return the name of the parser-builder in module_name, or None."""
    module = importlib.import_module(module_name)
    for candidate in parser_func_candidates(command_name):
        if hasattr(module, candidate):
            return candidate
    return None


def rst_header(text: str, char: str) -> str:
    return f"{text}\n{char * len(text)}"


def group_title(group_key: str) -> str:
    # "convertors_to_json" -> "Convertors To Json"
    return group_key.replace("_", " ").title()


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    group_pages: list[str] = []
    missing: list[str] = []

    for group, commands in REGISTRY.items():
        lines = [rst_header(group_title(group), "="), ""]

        for name, cmd in commands.items():
            module_name = cmd.func.__module__
            func_name = find_parser_func(module_name, name)

            if func_name is None:
                missing.append(f"{name} ({module_name})")
                continue

            lines.append(rst_header(name, "-"))
            lines.append("")
            lines.append(".. argparse::")
            lines.append(f"   :module: {module_name}")
            lines.append(f"   :func: {func_name}")
            lines.append(f"   :prog: {PROG} {name}")
            lines.append("")

        (OUTPUT_DIR / f"{group}.rst").write_text("\n".join(lines))
        group_pages.append(group)

    # write the index that ties the group pages together
    index_lines = [
        rst_header("Command-line reference", "="),
        "",
        ".. toctree::",
        "   :maxdepth: 2",
        "",
    ]
    index_lines += [f"   {g}" for g in group_pages]
    index_lines.append("")
    (OUTPUT_DIR / "index.rst").write_text("\n".join(index_lines))

    print(f"Wrote {len(group_pages)} group pages to {OUTPUT_DIR}")
    if missing:
        print("\nNo parser-builder found for these commands (skipped):")
        for m in missing:
            print(f"  - {m}")
        print("\nAdd a get_parser() to each, then re-run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
