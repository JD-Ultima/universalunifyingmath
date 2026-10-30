"""An exact-source audit must reject extra inputs even when Git ignores them."""
import importlib.util
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("research_verify", ROOT / "scripts/verify.py")
VERIFY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERIFY)


@unittest.skipUnless(shutil.which("git"), "Git is required to validate a source snapshot")
class SourceSnapshot(unittest.TestCase):
    def setUp(self):
        fixtures = (ROOT / ".verification-output").resolve()
        if not fixtures.is_relative_to(ROOT.resolve()):
            raise RuntimeError("Source fixtures must stay inside the repository")
        fixtures.mkdir(exist_ok=True)
        temporary = tempfile.TemporaryDirectory(prefix="source-snapshot-", dir=fixtures)
        self.addCleanup(self.cleanup_source, temporary, fixtures)
        self.source = Path(temporary.name)
        self.git("init", "--quiet")
        index = self.source / "compendium/index/statements.tsv"
        index.parent.mkdir(parents=True)
        index.write_text("number\ttype\n", encoding="utf-8")
        (self.source / ".gitignore").write_text("ignored_sources/\n", encoding="utf-8")
        self.git("add", ".gitignore", "compendium/index/statements.tsv")
        self.git("-c", "user.name=Research regression", "-c", "user.email=regression@example.invalid",
                 "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "Source fixture")
        expected_commit = self.git("rev-parse", "HEAD").strip()
        self.expected = mock.patch.object(VERIFY, "SOURCE_COMMIT", expected_commit)
        self.expected.start()
        self.addCleanup(self.expected.stop)

    @staticmethod
    def cleanup_source(temporary, fixtures):
        if not Path(temporary.name).resolve().is_relative_to(fixtures):
            raise RuntimeError("Refusing to clean up a source fixture outside its directory")
        temporary.cleanup()

    def git(self, *arguments):
        return subprocess.run(["git", "-C", str(self.source), *arguments], check=True,
                              capture_output=True, text=True, encoding="utf-8", timeout=30).stdout

    def test_pristine_recorded_commit_is_accepted(self):
        self.assertEqual(len(VERIFY.source_commands(self.source)), 4)

    def test_untracked_tex_source_is_rejected(self):
        extra = self.source / "compendium/tex/parts/extra.tex"
        extra.parent.mkdir(parents=True)
        extra.write_text("Additional text changes the recursively scanned source.\n", encoding="utf-8")
        # HEAD and all tracked files still match the fixture commit.
        self.assertEqual(self.git("status", "--porcelain", "--untracked-files=no"), "")
        with self.assertRaisesRegex(ValueError, "pristine.*untracked files"):
            VERIFY.source_commands(self.source)

    def test_ignored_lean_source_is_rejected(self):
        extra = self.source / "compendium/lean/ignored_sources/Extra.lean"
        extra.parent.mkdir(parents=True)
        extra.write_text("theorem extra : True := True.intro\n", encoding="utf-8")
        # Even the ordinary clean-checkout check hides this recursively scanned input.
        self.assertEqual(self.git("status", "--porcelain", "--untracked-files=all"), "")
        with self.assertRaisesRegex(ValueError, "pristine.*ignored files"):
            VERIFY.source_commands(self.source)


if __name__ == "__main__":
    unittest.main()
