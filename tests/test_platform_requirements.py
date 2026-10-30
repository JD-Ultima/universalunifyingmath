"""Unsupported platforms must fail before the heavyweight audits start work."""
import pathlib
import runpy
import sys
import unittest
from unittest import mock


ROOT = pathlib.Path(__file__).resolve().parents[1]
TRANSPORT = ROOT / "computation/chain-theory/alpha/transport.py"
LEAN_BUILD = ROOT / "computation/chain-theory/compendium/lean_build.py"


class PlatformRequirements(unittest.TestCase):
    def test_transport_without_fork_exits_before_numerical_setup(self):
        # This also works in a clean Python installation without NumPy/SymPy.
        with mock.patch("multiprocessing.get_all_start_methods", return_value=["spawn"]):
            with self.assertRaisesRegex(SystemExit, "Linux/WSL"):
                runpy.run_path(str(TRANSPORT), run_name="platform_check")

    def test_lean_build_on_windows_exits_before_launching_or_copying(self):
        build = runpy.run_path(str(LEAN_BUILD), run_name="platform_check")
        args = [str(LEAN_BUILD), "--lean", "missing", "--src", "missing", "--work", "missing"]
        with mock.patch.object(sys, "argv", args), mock.patch.object(sys, "platform", "win32"):
            with mock.patch("subprocess.run") as run, mock.patch.object(pathlib.Path, "mkdir") as mkdir:
                with self.assertRaisesRegex(SystemExit, "Linux/WSL"):
                    build["main"]()
                run.assert_not_called()
                mkdir.assert_not_called()

    def test_lean_build_help_remains_available_on_windows(self):
        build = runpy.run_path(str(LEAN_BUILD), run_name="platform_check")
        with mock.patch.object(sys, "argv", [str(LEAN_BUILD), "--help"]):
            with mock.patch.object(sys, "platform", "win32"), mock.patch("sys.stdout"):
                with self.assertRaises(SystemExit) as caught:
                    build["main"]()
                self.assertEqual(caught.exception.code, 0)


if __name__ == "__main__":
    unittest.main()
