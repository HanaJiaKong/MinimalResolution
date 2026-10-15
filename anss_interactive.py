#!/usr/bin/env python3
"""Export MinimalResolution tables to an offline, interactive HTML chart.

    python3 anss_interactive.py 30 -d output/alpha1-t30-l8 --all

Uses the same parser and cyclic-tower grouping as anss_chart.py. No packages,
server, or internet connection are required to generate or open the result.
"""

import argparse
from collections import Counter
from html import escape
import json
import re
from pathlib import Path
import sys

from anss_chart import (
    CLASS_RE, MULT_RE, NAME, PRIME, build_marks, looks_characteristic_p, parse_mult_detail,
    parse_run_info, parse_table, pretty,
)

ROOT = Path(__file__).resolve().parent
SOURCE_COMMIT = "https://github.com/yqkewang/ANSS_prime3/commit/f995e30f327ce83304eda57768d2c8b4260b1fe3"
SUBJECTS = {
    "sphere": "Sphere", "alpha_1": "Cα₁", "mod_p": "S/3",
    "mod_p_v1": "V(1)", "alpha_1_mod_p": "Cα₁ ∧ S/3", "triv_01": "S ∨ S¹",
}


def run_prefix(directory, degree, subject):
    """Match the static chart's legacy sphere fallback, without changing files."""
    prefix = directory / f"{degree}_{subject}BP"
    legacy = directory / f"{degree}_BP"
    if subject == "sphere" and not Path(str(prefix) + "AANSS_table.txt").exists():
        prefix = legacy
    return str(prefix)


def load_run(directory, degree, subject, height_override=None):
    base = run_prefix(directory, degree, subject)
    table = Path(base + "AANSS_table.txt")
    if not table.exists():
        raise ValueError(f"Missing {table}. Compute this comodule first, or choose --dir.")
    # parse_table historically skips malformed lines. An interactive export should
    # never quietly display an incomplete page as though all data was consumed.
    for number, line in enumerate(table.read_text().splitlines(), 1):
        if line and not CLASS_RE.fullmatch(line):
            raise ValueError(f"Unrecognized class at {table}:{number}: {line!r}")
    classes, differentials = parse_table(table)
    if not classes:
        raise ValueError(f"No classes in {table}.")
    info = parse_run_info(base + "run_info.txt") or {}
    warnings = []
    if height_override is not None:
        height = height_override
        characteristic = PRIME if height else 0
    elif "height" in info:
        height = int(info["height"])
        characteristic = int(info.get("characteristic", PRIME if height else 0))
    elif looks_characteristic_p(classes):
        raise ValueError(f"{base}run_info.txt is missing; use --height to identify the base ring.")
    else:
        height, characteristic = 0, 0
        warnings.append("Run metadata is missing. Assumed height 0 from the v₀ classes.")
    if height < 0 or characteristic != (PRIME if height else 0):
        raise ValueError(f"Inconsistent height/characteristic in {base}run_info.txt.")
    if "max_degree" in info and int(info["max_degree"]) != degree:
        raise ValueError(f"Run metadata does not match requested internal degree {degree}.")
    length = int(info["resolution_length"]) if "resolution_length" in info else None
    products = {}
    for kind, suffix in [("alpha1", "h0"), ("scalar", f"a{height}")]:
        path = base + f"AANSS_{suffix}.txt"
        if Path(path).exists():
            for number, line in enumerate(Path(path).read_text().splitlines(), 1):
                if not line:
                    continue
                match = MULT_RE.fullmatch(line)
                if not match:
                    raise ValueError(f"Unrecognized product at {path}:{number}: {line!r}")
                targets = [t for t in match["targets"].split("+") if t]
                if any(not re.fullmatch(NAME, t) for t in targets):
                    raise ValueError(f"Invalid product target at {path}:{number}.")
                if any(name not in classes for name in [match["src"], *targets]):
                    raise ValueError(f"Product references an absent class at {path}:{number}.")
        products[kind] = parse_mult_detail(path)
        if products[kind] is None:
            warnings.append(f"{Path(path).name} is missing; these products are unknown.")
    return classes, differentials, products, height, characteristic, length, warnings


def make_charts(directory, degree, subject, height_override=None):
    classes, diffs, products, height, char, length, warnings = load_run(
        directory, degree, subject, height_override)
    scalar = products["scalar"]
    known_towers = char != 0 or scalar is not None
    marks, owner = build_marks(classes, None if char else scalar,
                               ext0_torsion_free=not char)
    per_cell = Counter((m["stem"], m["s"]) for m in marks)
    nodes = []
    for mark in marks:
        node = dict(mark)
        node.update(id=mark["name"], label=pretty(mark["name"]), x=mark["stem"],
                    y=mark["s"], t=mark["stem"] + mark["s"], killed=False,
                    cohabitants=1 if char else per_cell[mark["stem"], mark["s"]])
        if not known_towers and mark["s"] > 0:
            # No scalar table means no evidence for any finite group order.
            node.update(group="Unknown (multiplication by 3 table missing)",
                        kind="circles", truncated=True)
        nodes.append(node)
    title = SUBJECTS.get(subject, subject)
    metadata = dict(subject=subject, subject_label=title, max_degree=degree,
                    resolution_length=length, height=height, characteristic=char,
                    warnings=warnings, scalar_label="3" if height == 0 else f"v{height}",
                    source_prefix=Path(run_prefix(directory, degree, subject)).name)

    def edges_for(node_ids, owners=None):
        edges, seen = [], set()
        for kind, mapping in products.items():
            if owners is not None and kind == "scalar":
                continue  # 3-towers are encoded in the cyclic-summand glyph.
            for source, targets in (mapping or {}).items():
                for target in targets:
                    if owners is not None:
                        if source not in owners or target not in owners:
                            continue
                        src, dst = marks[owners[source]]["name"], marks[owners[target]]["name"]
                    else:
                        src, dst = source, target
                    key = src, dst, kind
                    if src in node_ids and dst in node_ids and src != dst and key not in seen:
                        seen.add(key)
                        edges.append(dict(source=src, target=dst, type=kind,
                                          label="α₁" if kind == "alpha1" else metadata["scalar_label"]))
        return edges

    anss = dict(metadata, id=f"{subject}:anss", title=f"{title} · Adams–Novikov E₂",
                grading="anss", nodes=nodes,
                edges=edges_for({n["id"] for n in nodes}, owner))
    if char:
        # v_n is a meaningful product in characteristic 3, even on the ANSS
        # view; unlike 3 at height 0 it is not absorbed into the glyphs.
        for edge in edges_for({n["id"] for n in nodes}):
            if edge["type"] == "scalar":
                anss["edges"].append(edge)
    alg_nodes = [dict(c, id=c["name"], label=pretty(c["name"]), x=c["stem"],
                      y=c["algnov"], t=c["stem"] + c["s"], kind="dot", n=1,
                      group="Associated-graded basis class", truncated=False,
                      members=[c["name"]], cohabitants=1) for c in classes.values()]
    ids = {n["id"] for n in alg_nodes}
    alg_edges = edges_for(ids)
    alg_edges.extend(dict(source=source, target=target, type="differential", label=f"d{r}", r=r)
                     for source, target, r in diffs if source in ids and target in ids)
    algnov = dict(metadata, id=f"{subject}:algnov", title=f"{title} · Algebraic Novikov",
                  grading="algnov", nodes=alg_nodes, edges=alg_edges)
    return [anss, algnov]


def discover_subjects(directory, degree):
    suffix, subjects = "BPAANSS_table.txt", []
    for table in sorted(directory.glob(f"{degree}_*{suffix}")):
        name = table.name[len(str(degree)) + 1:-len(suffix)]
        subjects.append(name or "sphere")
    return list(dict.fromkeys(subjects))


def render_html(data):
    # JSON lives in an inert script element. Escape '<' so user-supplied names
    # cannot terminate the element, even in a locally opened HTML file.
    payload = json.dumps(data, ensure_ascii=False).replace("<", "\\u003c").replace(
        "\u2028", "\\u2028").replace("\u2029", "\\u2029")
    template = (ROOT / "visualizer" / "template.html").read_text()
    license_text = (ROOT / "visualizer" / "SEQSEE_LICENSE.txt").read_text()
    page = template.replace("__PAGE_TITLE__", escape("MinimalResolution · Interactive charts"))
    page = page.replace("__CHART_DATA__", payload)
    return "<!-- SeqSee attribution and license\n" + license_text + "\n-->\n" + page


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("degree", type=int, help="maximum internal degree t, matching the run filenames")
    parser.add_argument("-d", "--dir", type=Path, default=Path("."), help="directory of completed resolution tables")
    parser.add_argument("-c", "--comodule", action="append", help="comodule name; repeat to include several (default: sphere)")
    parser.add_argument("--all", action="store_true", help="include every matching completed comodule in --dir")
    parser.add_argument("--height", type=int, help="explicit base-ring height, useful for older runs without metadata")
    parser.add_argument("-o", "--output", type=Path, help="HTML output (default: <dir>/<degree>_interactive.html)")
    parser.add_argument("--json", type=Path, help="also write the chart dataset as JSON")
    args = parser.parse_args(argv)
    if args.degree <= 0 or (args.height is not None and args.height < 0):
        parser.error("degree must be positive and height must be nonnegative")
    if args.all and args.comodule:
        parser.error("choose --all or --comodule")
    subjects = discover_subjects(args.dir, args.degree) if args.all else (args.comodule or ["sphere"])
    if not subjects:
        parser.error(f"no completed {args.degree}_*BPAANSS_table.txt runs in {args.dir}")
    data = dict(format="minimalresolution-chart-v1", prime=PRIME,
                attribution=dict(name="SeqSee / ANSS_prime3 adaptation", source=SOURCE_COMMIT), charts=[])
    try:
        for subject in dict.fromkeys(subjects):
            data["charts"].extend(make_charts(args.dir, args.degree, subject, args.height))
        html = render_html(data)
        output = args.output or args.dir / f"{args.degree}_interactive.html"
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(html, encoding="utf-8")
        if args.json:
            args.json.parent.mkdir(parents=True, exist_ok=True)
            args.json.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    except (ValueError, OSError) as error:
        parser.exit(1, f"error: {error}\n")
    print(f"Wrote {output} ({len(data['charts'])} views; open directly in a browser)")
    for chart in data["charts"]:
        print(f"  {chart['title']}: {len(chart['nodes'])} marks, {len(chart['edges'])} edges")
    return 0


if __name__ == "__main__":
    sys.exit(main())
