"""Raw source pins distinguish canonical LF bytes from a CRLF checkout."""
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("research_verify_pin", ROOT / "scripts/verify.py")
VERIFY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERIFY)


class AlphaSourcePin(unittest.TestCase):
    def setUp(self):
        parent = (ROOT / ".verification-output").resolve()
        if ROOT.resolve() not in parent.parents:
            raise RuntimeError("fixture directory must stay inside the repository")
        parent.mkdir(exist_ok=True)
        temporary = tempfile.TemporaryDirectory(prefix="alpha-pin-", dir=parent)
        self.directory = Path(temporary.name).resolve()
        if parent not in self.directory.parents:
            raise RuntimeError("fixture escaped the verification directory")
        self.addCleanup(temporary.cleanup)
        self.source = self.directory / "papers/alpha/lean/ChainAlpha.lean"
        self.source.parent.mkdir(parents=True)
        self.lf = b"-- canonical source fixture\nexample : True := True.intro\n"
        pin = mock.patch.object(VERIFY, "LEAN_SOURCE_SHA256", hashlib.sha256(self.lf).hexdigest())
        pin.start()
        self.addCleanup(pin.stop)

    def test_exact_lf_bytes_pass(self):
        self.source.write_bytes(self.lf)
        VERIFY.check_alpha_source_bytes(self.directory)

    def test_equivalent_crlf_text_is_rejected_with_checkout_help(self):
        self.source.write_bytes(self.lf.replace(b"\n", b"\r\n"))
        with self.assertRaisesRegex(ValueError, "core.autocrlf=false"):
            VERIFY.check_alpha_source_bytes(self.directory)


if __name__ == "__main__":
    unittest.main()
