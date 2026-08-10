"""Regression tests for the Python 3 LDSC integration."""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from postgwas.modules.ldsc.ldsc_runner import run_ldsc


class LDSCRuntimeTests(unittest.TestCase):
    def test_runner_uses_active_python_and_distinct_weight_scores(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sumstats = root / "input.tsv"
            hm3 = root / "w_hm3.snplist"
            reference = root / "reference"
            weights = root / "weights"
            output = root / "results" / "study"
            sumstats.write_text("SNP\tN\tZ\nrs1\t1000\t1\n", encoding="utf-8")
            hm3.write_text("SNP\tA1\tA2\nrs1\tA\tG\n", encoding="utf-8")
            reference.mkdir()
            weights.mkdir()
            output.parent.mkdir()
            Path(f"{output}.sumstats.gz").touch()

            with patch(
                "postgwas.modules.ldsc.ldsc_runner._installed_script",
                side_effect=["/env/bin/munge_sumstats.py", "/env/bin/ldsc.py"],
            ), patch(
                "postgwas.modules.ldsc.ldsc_runner.run_subprocess"
            ) as run:
                run_ldsc(
                    sumstats_tsv=str(sumstats),
                    out_prefix=str(output),
                    hm3_snplist=str(hm3),
                    ldscore_dir=str(reference),
                    weight_ldscore_dir=str(weights),
                )

            munge_command = run.call_args_list[0].args[0]
            h2_command = run.call_args_list[1].args[0]
            self.assertEqual(munge_command[:2], [sys.executable, "/env/bin/munge_sumstats.py"])
            self.assertEqual(h2_command[:2], [sys.executable, "/env/bin/ldsc.py"])
            self.assertEqual(
                h2_command[h2_command.index("--ref-ld-chr") + 1],
                str(reference.resolve()) + "/",
            )
            self.assertEqual(
                h2_command[h2_command.index("--w-ld-chr") + 1],
                str(weights.resolve()) + "/",
            )

    def test_dockerfile_uses_cbiit_ldsc_in_postgwas_environment(self):
        dockerfile = (Path(__file__).parents[1] / "Dockerfile").read_text(encoding="utf-8")
        self.assertIn("https://github.com/CBIIT/ldsc.git", dockerfile)
        self.assertIn("LDSC_COMMIT=6c673952cee74bd5c57aef1555a03b1c015399a0", dockerfile)
        self.assertIn("pip install --no-deps --no-cache-dir /opt/ldsc", dockerfile)
        self.assertNotIn("python=2.7", dockerfile)
        self.assertNotIn("micromamba create -y -n ldsc", dockerfile)
        self.assertNotIn("github.com/bulik/ldsc", dockerfile)


if __name__ == "__main__":
    unittest.main()
