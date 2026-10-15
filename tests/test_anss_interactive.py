"""Data-contract checks for the offline viewer (stdlib unittest only)."""

import json
from pathlib import Path
import re
import tempfile
import unittest

import anss_interactive as viewer


class InteractiveChartTests(unittest.TestCase):
    def test_reference_alpha1_groups_and_differentials(self):
        charts = viewer.make_charts(viewer.ROOT / "example_data", 30, "alpha_1")
        anss, algnov = charts
        self.assertEqual(len(anss["nodes"]), 12)
        self.assertEqual(sum(n["kind"] == "box" for n in anss["nodes"]), 2)
        self.assertEqual(sum(n["truncated"] for n in anss["nodes"]), 1)
        self.assertEqual(len(anss["edges"]), 3)
        self.assertTrue(all(e["type"] == "alpha1" for e in anss["edges"]))
        self.assertFalse(any(n["killed"] for n in anss["nodes"]))
        self.assertEqual(len(algnov["nodes"]), 45)
        differential_edges = [e for e in algnov["edges"] if e["type"] == "differential"]
        self.assertEqual(len(differential_edges), 5)
        lookup = {n["id"]: n for n in algnov["nodes"]}
        for edge in differential_edges:
            source, target = lookup[edge["source"]], lookup[edge["target"]]
            self.assertEqual(target["x"], source["x"] - 1)
            self.assertEqual(target["y"], source["y"] + edge["r"])
            self.assertTrue(source["killed"] and target["killed"])
        for chart in charts:
            ids = {n["id"] for n in chart["nodes"]}
            self.assertTrue(all(e["source"] in ids and e["target"] in ids for e in chart["edges"]))
            for node in chart["nodes"]:
                self.assertEqual(node["t"], node["x"] + node["s"])

    def write_run(self, path, height=1):
        base = path / "16_testBP"
        Path(str(base) + "AANSS_table.txt").write_text("[0-0]\t|deg=(0,0)\nv1^1[0-0]\t|deg=(4,1)\n")
        Path(str(base) + "run_info.txt").write_text(
            f"height {height}\ncharacteristic {3 if height else 0}\nmax_degree 16\nresolution_length 4\n")
        Path(str(base) + f"AANSS_a{height}.txt").write_text(
            "[0-0]\t->\tv1^1[0-0]+o\nv1^1[0-0]\t->\to\n")
        Path(str(base) + "AANSS_h0.txt").write_text("[0-0]\t->\to\n")
        return base

    def test_characteristic_three_preserves_vn_products_without_boxes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            self.write_run(path)
            anss, _ = viewer.make_charts(path, 16, "test")
            self.assertEqual(len(anss["nodes"]), 2)
            self.assertTrue(all(n["kind"] == "dot" and n["group"] == "Z/3" for n in anss["nodes"]))
            self.assertEqual([e["type"] for e in anss["edges"]], ["scalar"])
            self.assertEqual(anss["edges"][0]["label"], "v1")

    def test_missing_ring_metadata_and_bad_rows_fail_clearly(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            base = self.write_run(path)
            Path(str(base) + "run_info.txt").unlink()
            with self.assertRaisesRegex(ValueError, "use --height"):
                viewer.make_charts(path, 16, "test")
            self.assertEqual(viewer.make_charts(path, 16, "test", 1)[0]["height"], 1)
            Path(str(base) + "AANSS_table.txt").write_text("broken row\n")
            with self.assertRaisesRegex(ValueError, "Unrecognized class"):
                viewer.make_charts(path, 16, "test", 1)

    def test_absent_scalar_table_does_not_invent_finite_group_orders(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            base = self.write_run(path, height=0)
            Path(str(base) + "AANSS_table.txt").write_text("[1-0]\t|deg=(3,1)\n")
            Path(str(base) + "AANSS_h0.txt").write_text("[1-0]\t->\to\n")
            Path(str(base) + "AANSS_a0.txt").unlink()
            chart = viewer.make_charts(path, 16, "test")[0]
            self.assertIn("Unknown", chart["nodes"][0]["group"])
            self.assertTrue(chart["warnings"])

    def test_malformed_and_dangling_products_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            base = self.write_run(path)
            product = Path(str(base) + "AANSS_h0.txt")
            product.write_text("bad product\n")
            with self.assertRaisesRegex(ValueError, "Unrecognized product"):
                viewer.make_charts(path, 16, "test")
            product.write_text("[0-0]\t->\t[9-9]+o\n")
            with self.assertRaisesRegex(ValueError, "absent class"):
                viewer.make_charts(path, 16, "test")

    def test_discovery_legacy_sphere_and_safe_html_roundtrip(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            base = self.write_run(path)
            legacy = path / "16_BP"
            for file in path.iterdir():
                file.rename(Path(str(file).replace(str(base), str(legacy))))
            self.assertEqual(viewer.discover_subjects(path, 16), ["sphere"])
            self.assertEqual(viewer.make_charts(path, 16, "sphere")[0]["subject"], "sphere")
        payload = {"format": "minimalresolution-chart-v1", "charts": [],
                   "unsafe_label": '</script><img src=x onerror="boom">'}
        html = viewer.render_html(payload)
        self.assertNotIn(payload["unsafe_label"], html)
        embedded = re.search(r'<script[^>]+id="chart-data"[^>]*>(.*?)</script>', html, re.S)
        self.assertIsNotNone(embedded)
        self.assertEqual(json.loads(embedded.group(1)), payload)
        self.assertNotRegex(html, r'<(?:script|link)\b[^>]*(?:src|href)="https?://')


if __name__ == "__main__":
    unittest.main()
