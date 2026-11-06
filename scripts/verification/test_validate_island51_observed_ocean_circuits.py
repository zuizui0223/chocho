"""Run the portable Julia boundary against hand-solved distributed circuits.

Set JULIA_BINARY and CIRCUITSCAPE_PROJECT to an installed pinned runtime.
The tests exercise real Circuitscape; the ordinary Python suite need not install it.
"""
import csv
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("validate_island51_observed_ocean_circuits.jl")
JULIA = os.environ.get("JULIA_BINARY")
PROJECT = os.environ.get("CIRCUITSCAPE_PROJECT")


@unittest.skipUnless(JULIA and PROJECT, "Pinned Julia/Circuitscape environment is optional")
class OfficialDistributedCircuits(unittest.TestCase):
    def fixture(self, root, *, wrong_reference=False):
        networks, results, out = (root / key for key in ("networks", "results", "out"))
        networks.mkdir()
        results.mkdir()
        for field, conductance in [("uniform", 1), ("ais", 2), ("wind", 4)]:
            (networks / f"{field}_conductances.txt").write_text(
                f"1 2 {conductance}\n2 3 {conductance}\n3 4 {conductance}\n")
        # Duplicate rows for node1 must add, not overwrite each other. The
        # distributed chain energy is .5^2+1^2+.5^2=1.5; endpoint energy is3.
        (networks / "island_node_weights.csv").write_text(
            "landmass_id,network_node_id,source_graph_node_id,weight\n"
            "9,1,0,0.25\n9,1,0,0.25\n9,2,1,0.5\n"
            "15,3,2,0.5\n15,4,3,0.5\n715,1,0,1\n1577,4,3,1\n")
        first = 4. if wrong_reference else 1.5
        (results / "island51_observed_circuit_pairs.csv").write_text(
            "island_a,island_b,uniform_distributed_port_resistance,"
            "ais_distributed_port_resistance,wind_distributed_port_resistance\n"
            f"9,15,{first},0.75,0.375\n715,1577,3,1.5,0.75\n")
        return networks, results, out

    def run_julia(self, directories):
        return subprocess.run([JULIA, f"--project={PROJECT}", str(SCRIPT),
                               *(str(path) for path in directories)],
                              text=True, capture_output=True, timeout=120)

    def test_cli_aggregates_port_currents_and_matches_six_analytic_values(self):
        with tempfile.TemporaryDirectory() as temp:
            directories = self.fixture(Path(temp))
            completed = self.run_julia(directories)
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            with (directories[2] / "official_checks.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 6)
            expected = {("uniform", "9"): 1.5, ("ais", "9"): .75, ("wind", "9"): .375,
                        ("uniform", "715"): 3., ("ais", "715"): 1.5, ("wind", "715"): .75}
            for row in rows:
                self.assertAlmostEqual(float(row["official_resistance"]),
                                       expected[row["field"], row["island_a"]], places=10)
                self.assertLess(float(row["relative_Kirchhoff_residual"]), 1e-10)
            self.assertTrue((directories[2] / "verification_summary.toml").is_file())

    def test_cli_rejects_a_matrix_reference_inconsistent_with_actual_network(self):
        with tempfile.TemporaryDirectory() as temp:
            directories = self.fixture(Path(temp), wrong_reference=True)
            completed = self.run_julia(directories)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("resistance disagrees with the reference pair table", completed.stderr)
            self.assertFalse((directories[2] / "verification_summary.toml").exists())


if __name__ == "__main__":
    unittest.main()
