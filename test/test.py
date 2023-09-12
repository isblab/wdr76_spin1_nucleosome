#!/usr/bin/env python
import unittest
import os
import shutil
import sys
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), "../"))


class Tests(unittest.TestCase):
    def test_sims(self):
        """Test WDR76-SPIN1-nucleosome simulations"""
        # os.chdir(os.path.join(TOPDIR, "scripts"))
        if os.path.exists("output"):
            shutil.rmtree("output")

        p = subprocess.check_call(
            ["python", "../../scripts/modeling.py", "test", "test"]
        )

        # require that the number of frames is present
        total_num_lines_stat_files = 0
        for i in range(1):
            with open("run_test/stat." + str(i) + ".out", "r") as statf:
                total_num_lines_stat_files += len(statf.readlines())
        self.assertEqual(total_num_lines_stat_files, 51)

        # require that output files were produced
        for i in range(1):
            os.unlink("run_test/rmfs/" + str(i) + ".rmf3")
            os.unlink("run_test/stat." + str(i) + ".out")
            os.unlink("run_test/stat_replica." + str(i) + ".out")


if __name__ == "__main__":
    unittest.main()
