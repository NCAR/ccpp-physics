#!/usr/bin/env python3
"""Generate HTML argument tables from CCPP metadata files.

This is a self-contained replacement for the metadata2html.py script of the
retired CCPP prebuild framework. It reads every ``.meta`` file tracked in a
ccpp-physics checkout and writes one ``<section name>.html`` table for each
metadata section that declares variables. The rendered values follow the
conventions of the former script so the published tables stay unchanged.
"""

import argparse
from html import escape
from pathlib import Path
import re
import subprocess


ATTRIBUTES = (
    "local_name",
    "standard_name",
    "long_name",
    "units",
    "type",
    "dimensions",
    "kind",
    "intent",
)

# Host-specific variants share section names with the SCM variants. As in the
# former SCM prebuild configuration, only the SCM variants are rendered.
HOST_VARIANT_SUFFIXES = (".fv3.meta", ".mpas.meta", ".neptune.meta")

FORTRAN_ID = r"[A-Za-z][A-Za-z0-9_]*"
TABLE_START = re.compile(r"(?i)\[\s*ccpp-table-properties\s*\]$")
SECTION_START = re.compile(r"(?i)\[\s*ccpp-arg-table\s*\]$")
VARIABLE_START = re.compile(r"\[\s*(" + FORTRAN_ID + r"(?:\s*\([^]]*\))?)\s*\]$")
ARRAY_REFERENCE = re.compile(r"(" + FORTRAN_ID + r")\s*\((.*)\)$")
REAL_NUMBER = re.compile(r"(.*\d)p(\d.*)")
DOUBLE_PRECISION = re.compile(r"(?i)double\s*precision$")
INTRINSIC_TYPES = ("integer", "real", "logical", "complex", "character")


class MetadataError(Exception):
    """A metadata file cannot be converted to HTML."""


def read_lines(path):
    """Yield (line number, line) for non-blank lines, joining continuations."""
    pending = ""
    start = None
    with path.open(encoding="utf-8") as source:
        for number, line in enumerate(source, start=1):
            line = line.strip()
            if start is None:
                start = number
            if line.endswith("\\"):
                pending += line[:-1]
                continue
            line = pending + line
            pending = ""
            if line and not line.startswith(("#", ";")):
                yield start, line
            start = None
    if pending:
        yield start, pending


def parse_properties(line, where):
    """Parse ``key = value`` pairs, several of which may share a line."""
    properties = []
    for item in line.split("|"):
        key, separator, value = item.partition("=")
        if not separator or not key.strip():
            raise MetadataError(f"{where}: expected 'key = value', found {item.strip()!r}")
        properties.append((key.strip().lower(), value.strip()))
    return properties


def parse_dimensions(value, where):
    """Return dimension ranges as the former framework reported them."""
    value = value.strip()
    if not (value.startswith("(") and value.endswith(")")):
        raise MetadataError(f"{where}: dimensions must be a parenthesized list")
    dimensions = []
    for dimension in value[1:-1].split(","):
        dimension = dimension.strip()
        if not dimension:
            continue
        if ":" not in dimension:
            dimension = "ccpp_constant_one:" + dimension
        dimensions.append(dimension)
    return dimensions


def is_intrinsic(type_name):
    type_name = type_name.lower()
    return type_name in INTRINSIC_TYPES or DOUBLE_PRECISION.match(type_name)


def default_long_name(standard_name):
    """Derive a long name from a standard name, e.g. 0p55mu becomes 0.55mu."""
    if not standard_name:
        return ""
    long_name = standard_name[0].upper() + standard_name[1:].replace("_", " ")
    match = REAL_NUMBER.match(long_name)
    while match:
        long_name = match.group(1) + "." + match.group(2)
        match = REAL_NUMBER.match(long_name)
    return long_name


def default_kind(type_name):
    type_name = type_name.lower()
    if type_name in ("real", "complex") or DOUBLE_PRECISION.match(type_name):
        return "kind_phys"
    return ""


def finish_variable(variable, where):
    """Apply the defaults and normalizations of the former framework."""
    type_name = variable.get("type")
    if type_name is not None and not is_intrinsic(type_name):
        # Derived types are reported with the type name as their kind.
        variable["kind"] = type_name
        variable.setdefault("units", "")
    elif type_name is not None and "kind" not in variable:
        variable["kind"] = default_kind(type_name)
    if "long_name" not in variable and "standard_name" in variable:
        variable["long_name"] = default_long_name(variable["standard_name"])

    reference = ARRAY_REFERENCE.match(variable["local_name"])
    if reference:
        indices = [index.strip() for index in reference.group(2).split(",")]
        rank = len(variable.get("dimensions", ()))
        if indices.count(":") != rank:
            raise MetadataError(
                f"{where}: {variable['local_name']} does not match its rank {rank}"
            )
        variable["local_name"] = f"{reference.group(1)}({', '.join(indices)})"
    return variable


def parse_metadata_file(path):
    """Return a list of (section name, variables) for a metadata file."""
    sections = []
    section = None
    variable = None
    state = None
    for number, line in read_lines(path):
        where = f"{path}:{number}"
        if TABLE_START.match(line):
            state, section, variable = "table", None, None
            continue
        if SECTION_START.match(line):
            if state is None:
                raise MetadataError(f"{where}: [ccpp-arg-table] before [ccpp-table-properties]")
            state, variable = "section", None
            section = {"name": None, "variables": []}
            sections.append(section)
            continue
        start = VARIABLE_START.match(line)
        if start:
            if state not in ("section", "variable"):
                raise MetadataError(f"{where}: variable outside of [ccpp-arg-table]")
            state = "variable"
            variable = {"local_name": start.group(1)}
            section["variables"].append((variable, where))
            continue
        if state is None:
            raise MetadataError(f"{where}: property outside of a metadata table")
        for key, value in parse_properties(line, where):
            if state == "section" and key == "name":
                section["name"] = value
            elif state == "variable":
                if key == "dimensions":
                    value = parse_dimensions(value, where)
                elif key == "standard_name":
                    value = value.lower()
                variable[key] = value

    tables = []
    for section in sections:
        if not section["name"]:
            raise MetadataError(f"{path}: [ccpp-arg-table] without a name")
        variables = [finish_variable(*entry) for entry in section["variables"]]
        tables.append((section["name"], variables))
    return tables


def render_table(name, variables):
    """Return the HTML document for one metadata section."""
    header = "<tr>" + "".join(f"<th>{attribute}</th>" for attribute in ATTRIBUTES) + "</tr>\n"
    rows = []
    for variable in variables:
        cells = []
        for attribute in ATTRIBUTES:
            value = variable.get(attribute)
            if attribute == "dimensions" and value is not None:
                value = "(" + ", ".join(value) + ")"
            elif value is None:
                value = "n/a"
            cells.append(f"<td>{escape(value, quote=False)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>\n")
    return (
        "\n<html>\n<head>\n"
        f"<title>{name} argument table</title>\n"
        '<meta charset="UTF-8">\n</head>\n<body>\n<table border="1">\n'
        + header
        + "".join(rows)
        + "</table>\n</body>\n</html>\n"
    )


def generate(physics, output):
    """Generate one HTML table for every metadata section with variables."""
    output.mkdir(parents=True, exist_ok=True)
    for old_table in output.glob("*.html"):
        old_table.unlink()

    files = subprocess.check_output(
        ["git", "ls-files", "--recurse-submodules", "-z", "--", "*.meta"],
        cwd=physics,
    ).decode().split("\0")
    generated = {}
    for filename in sorted(filter(None, files)):
        if filename.endswith(HOST_VARIANT_SUFFIXES):
            continue
        source = physics / filename
        for name, variables in parse_metadata_file(source):
            if not variables:
                continue
            if not re.fullmatch(FORTRAN_ID, name):
                raise MetadataError(f"{source}: invalid section name {name!r}")
            if name in generated:
                raise MetadataError(
                    f"Duplicate section {name}: {generated[name]} and {source}"
                )
            (output / f"{name}.html").write_text(
                render_table(name, variables), encoding="utf-8"
            )
            generated[name] = source

    if not generated:
        raise MetadataError(f"No metadata tables found in {physics}")
    print(f"Generated {len(generated)} metadata HTML tables in {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--physics", required=True, type=Path,
                        help="ccpp-physics checkout to read .meta files from")
    parser.add_argument("--output", required=True, type=Path,
                        help="directory to write the HTML tables to")
    arguments = parser.parse_args()
    try:
        generate(arguments.physics.resolve(), arguments.output.resolve())
    except MetadataError as error:
        raise SystemExit(f"error: {error}")


if __name__ == "__main__":
    main()
