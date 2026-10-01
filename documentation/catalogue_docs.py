"""Build the output variable catalogue page of the user guide.

The page is made from three source files, so it cannot drift from the code:

* ``output_catalogue.yaml``: what each output variable is (name, units, long name,
  dimensions, group, module, restart name, ...).
* ``cable_output_bindings*.F90``: which model variable each name is bound to, how it
  is converted for output, and when it is available.
* ``cable_checks.F90``: the valid range of each variable.

`render_catalogue` returns Markdown. It is called from ``main.py`` (the hook the
mkdocs ``macros`` plugin loads) and can also be run on its own to preview the page:

    python catalogue_docs.py > catalogue.md
"""
import re
import sys
from pathlib import Path

import yaml

# Where the sources live relative to the directory containing mkdocs.yml.
CATALOGUE_PATH = "../src/util/output/output_catalogue.yaml"
BINDINGS_PATHS = ["../src/offline/cable_output_bindings.F90", "../src/offline/cable_output_bindings_casa.F90"]
RANGES_PATH = "../src/offline/cable_checks.F90"

# Headings for the modules named in the catalogue, in the order they are shown.
MODULE_TITLES = {
    "forcing": "Meteorological forcing",
    "biogeophysics": "Biogeophysics",
    "biogeochemistry": "Biogeochemistry",
    "parameters": "Parameters",
}

# What the dimension names in the catalogue mean.
DIMENSION_MEANING = {
    "patch": "one value per tile (patch) in each grid cell",
    "soil": "soil layers",
    "snow": "snow layers",
    "rad": "radiation bands (visible and near-infrared)",
    "plant_carbon_pools": "plant carbon pools (leaf, wood, fine root)",
    "soil_carbon_pools": "soil carbon pools (microbial, slow, passive)",
    "land_global": "land points of the whole grid",
}

# What the conditions under which a variable exists mean. The names are the
# arguments of the binding functions.
AVAILABILITY_MEANING = {
    "use_groundwater_model": "the groundwater model is on (`cable_user%gw_model`)",
    "use_popluc": "land-use change with POPLUC is on (`cable_user%POPLUC`)",
    "calculate_soil_albedo": "soil albedo is calculated (`calcsoilalbedo`)",
    "casa_enabled": "CASA-CNP is on (`l_casacnp`)",
}


def load_catalogue(path):
    """The list of catalogue entries, in file order."""
    with open(path) as handle:
        return yaml.safe_load(handle)["variables"]


def _split_top_level(text):
    """Split on commas that are not inside parentheses, brackets or quotes."""
    parts, depth, current, quote = [], 0, [], None
    for char in text:
        if quote:
            current.append(char)
            if char == quote:
                quote = None
            continue
        if char in "\"'":
            quote = char
        elif char in "([":
            depth += 1
        elif char in ")]":
            depth -= 1
        if char == "," and depth == 0:
            parts.append("".join(current).strip())
            current = []
        else:
            current.append(char)
    if "".join(current).strip():
        parts.append("".join(current).strip())
    return parts


def _calls(text, marker):
    """The argument text of every ``marker( ... )`` call, with continuation marks removed."""
    start = 0
    while True:
        start = text.find(marker, start)
        if start < 0:
            return
        index = start + len(marker)
        depth, quote = 1, None
        while depth:
            char = text[index]
            if quote:
                if char == quote:
                    quote = None
            elif char in "\"'":
                quote = char
            elif char == "(":
                depth += 1
            elif char == ")":
                depth -= 1
            index += 1
        yield re.sub(r"\s*&\s*\n\s*", " ", text[start + len(marker):index - 1])
        start = index


def _strip_comment(line):
    """The line without a trailing Fortran comment (a ``!`` that is not inside a string)."""
    quote = None
    for index, char in enumerate(line):
        if quote:
            if char == quote:
                quote = None
        elif char in "\"'":
            quote = char
        elif char == "!":
            return line[:index]
    return line


def parse_bindings(paths):
    """What each binding says about its variable, and the declared type of every argument.

    Returns ``(bindings, argument_types)``. ``bindings`` maps a catalogue name to a dict
    with ``source`` (the model variable), ``scale_by``, ``divide_by``, ``offset_by``,
    ``range`` and ``available`` (each a string or None). ``argument_types`` maps the
    name of a binding function argument (such as ``canopy``) to its Fortran type.
    """
    bindings, argument_types = {}, {}
    for path in paths:
        text = Path(path).read_text()
        # Comments would confuse the search for calls, so drop them first.
        text = "\n".join(_strip_comment(line) for line in text.splitlines())
        for match in re.finditer(r"\b(type\((\w+)\)|integer|real|logical)[^:\n]*::\s*(\w+)", text, re.I):
            declared = match.group(2) or match.group(1)
            argument_types[match.group(3)] = declared
        for call in _calls(text, "cable_output_binding_t("):
            arguments = {}
            for part in _split_top_level(call):
                key, _, value = part.partition("=")
                arguments[key.strip()] = value.strip()
            # Entries without an aggregator are the placeholders for variables that
            # do not exist (CASA switched off); the real entry is elsewhere in the file.
            if "aggregator" not in arguments or "name" not in arguments:
                continue
            source = re.fullmatch(r"new_aggregator\((.*)\)", arguments["aggregator"], re.S)
            bindings[arguments["name"].strip('"')] = {
                "source": source.group(1).strip() if source else arguments["aggregator"],
                "scale_by": arguments.get("scale_by"),
                "divide_by": arguments.get("divide_by"),
                "offset_by": arguments.get("offset_by"),
                "range": arguments.get("range"),
                "available": arguments.get("available"),
            }
    return bindings, argument_types


def parse_ranges(path):
    """The valid range of each variable, from the defaults of ``ranges_type``.

    Returns a dict from the (lower case) variable name to the text ``low to high``.
    """
    text = Path(path).read_text()
    start = text.upper().index("TYPE RANGES_TYPE")
    block = text[start:text.upper().index("END TYPE RANGES_TYPE")]
    ranges = {}
    for match in re.finditer(r"^\s*(\w+)\s*=\s*\[([^\]]*)\]", block, re.M):
        low, _, high = (value.strip() for value in match.group(2).partition(","))
        ranges[match.group(1).lower()] = f"{low} to {high}"
    return ranges


def _code(text):
    return f"`{text}`"


def _escape(text):
    """Make text safe inside a Markdown table cell."""
    return str(text).replace("|", "\\|").replace("\n", " ")


def _conversion(binding):
    """A short description of the unit conversion, or an empty string if there is none."""
    steps = []
    if binding["scale_by"]:
        steps.append(f"× {_code(binding['scale_by'])}")
    if binding["divide_by"]:
        steps.append(f"÷ {_code(binding['divide_by'])}")
    if binding["offset_by"]:
        steps.append(f"+ {_code(binding['offset_by'])}")
    return " ".join(steps)


def _valid_range(binding, ranges):
    """The numeric valid range if it can be found, otherwise the Fortran expression."""
    if not binding["range"]:
        return ""
    name = binding["range"].split("%")[-1].lower()
    return ranges.get(name, _code(binding["range"]))


def _bound_to(binding, argument_types):
    """The model variable, with the type of the structure it belongs to."""
    root = re.match(r"[A-Za-z_]\w*", binding["source"]).group(0)
    kind = argument_types.get(root)
    text = _code(binding["source"])
    return f"{text} ({_code(kind)})" if kind else text


def _notes(entry, binding, ranges, default_condition):
    """Everything else worth knowing about a variable, as short phrases."""
    notes = []
    if entry.get("parameter"):
        notes.append("parameter: written once, no time axis")
    if entry.get("native_frequency"):
        notes.append(f"updated by the model {entry['native_frequency']}: cannot be written more often")
    if entry.get("distributed") is False:
        notes.append("not distributed across processes")
    if binding:
        conversion = _conversion(binding)
        if conversion:
            notes.append(f"converted for output: {conversion}")
        valid = _valid_range(binding, ranges)
        if valid:
            notes.append(f"valid range: {valid}")
        condition = binding["available"] or default_condition
        if condition:
            notes.append("only available when " + AVAILABILITY_MEANING.get(condition, _code(condition)))
    if entry.get("restart_name"):
        notes.append(f"restart file name: {_code(entry['restart_name'])}")
    return "; ".join(notes)


def _table(entries, bindings, argument_types, ranges, conditions, columns):
    """A Markdown table with one row per entry."""
    lines = ["| " + " | ".join(title for title, _ in columns) + " |", "|" + "---|" * len(columns)]
    for entry in entries:
        binding = bindings.get(entry["name"])
        context = {
            "entry": entry, "binding": binding, "argument_types": argument_types, "ranges": ranges,
            "condition": conditions.get(entry["name"]),
        }
        cells = [_escape(function(context)) for _, function in columns]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def _name_cell(context):
    name = context["entry"]["name"]
    return f'<span id="{name.lower()}"></span>**{name}**'


def _long_name(context):
    return context["entry"].get("metadata", {}).get("long_name", "")


def _units(context):
    units = context["entry"].get("metadata", {}).get("units")
    return _code(units) if units else ""


def _dimensions(context):
    dimensions = context["entry"].get("dimensions", [])
    return ", ".join(_code(d) for d in dimensions) if dimensions else "scalar"


def _bound(context):
    binding = context["binding"]
    return _bound_to(binding, context["argument_types"]) if binding else "(no binding)"


def _other_notes(context):
    return _notes(context["entry"], context["binding"], context["ranges"], context["condition"])


def _group(context):
    group = context["entry"].get("group")
    return _code(group) if group else ""


def _restart(context):
    name = context["entry"].get("restart_name")
    return _code(name) if name else ""


def render_catalogue(project_dir):
    """The Markdown for the whole catalogue page, built from the source files.

    `project_dir` is the directory that contains ``mkdocs.yml``.
    """
    base = Path(project_dir)
    entries = load_catalogue(base / CATALOGUE_PATH)
    bindings, argument_types = parse_bindings([base / path for path in BINDINGS_PATHS])
    ranges = parse_ranges(base / RANGES_PATH)

    # Every variable in the CASA bindings file is only there when CASA is on, so give
    # those the CASA condition unless the binding has a more specific one.
    casa_names = set(parse_bindings([base / BINDINGS_PATHS[1]])[0])
    conditions = {name: "casa_enabled" for name in casa_names}

    output_entries = [e for e in entries if e.get("output", True)]
    restart_only = [e for e in entries if not e.get("output", True)]
    missing = [e["name"] for e in entries if e["name"] not in bindings]

    out = [
        f"This page lists the {len(output_entries)} variables that can be written to output files, and "
        f"the {len(restart_only)} that are only written to restart files. It is generated from the "
        "variable catalogue (`output_catalogue.yaml`) and the bindings to model variables "
        "(`cable_output_bindings*.F90`) every time the documentation is built, so it always matches the code.",
        "",
        "Use a variable's **Name** in the output configuration file. A group or module name selects every "
        "variable listed under it that is available in the model configuration being run.",
        "",
        "## How to read the tables",
        "",
        "| Column | Meaning |",
        "|---|---|",
        "| Name | The name to use in the configuration file. It is also the default name of the variable in the NetCDF file. |",
        "| Long name, Units | Written to the NetCDF file as the `long_name` and `units` attributes. |",
        "| Dimensions | The shape of the model data. Time is added for variables that change in time. |",
        "| Group | The group that selects this variable (see the configuration reference). |",
        "| Bound to | The model variable the output is taken from, with the type of the structure it lives in. |",
        "| Notes | Unit conversion applied for output, valid range used by range checking, when the variable exists, and other properties. |",
        "",
        "Dimensions:",
        "",
    ]
    out += [f"- {_code(name)}: {meaning}" for name, meaning in DIMENSION_MEANING.items()]
    out.append("")

    columns = [
        ("Name", _name_cell), ("Long name", _long_name), ("Units", _units), ("Dimensions", _dimensions),
        ("Group", _group), ("Bound to", _bound), ("Notes", _other_notes),
    ]
    modules = [name for name in MODULE_TITLES if any(e.get("module") == name for e in output_entries)]
    modules += sorted({e.get("module") for e in output_entries if e.get("module")} - set(MODULE_TITLES))
    for module in modules:
        members = [e for e in output_entries if e.get("module") == module]
        title = MODULE_TITLES.get(module, module.replace("_", " ").title())
        out += [f"## {title}", "", f"Module `{module}`: {len(members)} variables.", "",
                _table(members, bindings, argument_types, ranges, conditions, columns), ""]

    ungrouped = [e for e in output_entries if not e.get("module")]
    if ungrouped:
        out += ["## Not in any module", "",
                "These can only be requested by name.", "",
                _table(ungrouped, bindings, argument_types, ranges, conditions, columns), ""]

    if restart_only:
        restart_columns = [
            ("Name", _name_cell), ("Long name", _long_name), ("Units", _units), ("Dimensions", _dimensions),
            ("Bound to", _bound), ("Restart file name", _restart),
        ]
        out += ["## Restart-only variables", "",
                "These are model state saved at the end of a run so that it can be restarted. They cannot "
                "be selected in the output configuration. They are written with the model's own numeric type.", "",
                _table(restart_only, bindings, argument_types, ranges, conditions, restart_columns), ""]

    if missing:
        # Shown rather than hidden: it means the catalogue and the bindings disagree.
        out += ["!!! warning", "    These catalogue entries have no binding: " + ", ".join(_code(n) for n in missing), ""]
    return "\n".join(out)


if __name__ == "__main__":
    print(render_catalogue(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parent))
