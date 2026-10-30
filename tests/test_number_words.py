"""Regression cases for the printed-number errors disclosed in round 7."""
import pathlib
import sys
import unittest


COMPENDIUM = pathlib.Path(__file__).resolve().parents[1] / "computation/chain-theory/compendium"
sys.path.insert(0, str(COMPENDIUM))
import lean_resolve
import lean_tables
import number_words


class PrintedNumbers(unittest.TestCase):
    def test_cardinal_phrases(self):
        cases = {
            "four one-pair acts": {4, 1},
            "two hundred twenty-eight one-act reads": {228, 1},
            "four, thirty-two": {4, 32},
            "four; thirty-two": {4, 32},
            "two hundred and twenty-eight": {228},
            "two thousand and one": {2001},
            "twenty-four": {24},
            "zero one": {0, 1},
        }
        for text, expected in cases.items():
            with self.subTest(text=text):
                self.assertEqual(number_words.cardinals(text), expected)

    def test_ordinals_and_tex_labels(self):
        self.assertEqual(number_words.printed_numbers("twenty-fourth"), ({24}, {24}))
        self.assertEqual(number_words.printed_numbers("tenth and ninth"), ({10, 9}, {10, 9}))
        self.assertEqual(number_words.printed_numbers("TWENTY FOURTH"), ({24}, {24}))
        self.assertEqual(number_words.printed_numbers("24 twenty-fourth"), ({24}, set()))
        self.assertEqual(number_words.numbers_in(r"\label{box:999} four \textit{one-pair} acts"), {4, 1})

    def test_both_audits_share_the_corrected_parser(self):
        for text in ("four one-pair acts", "four, thirty-two", "twenty-fourth"):
            with self.subTest(text=text):
                self.assertEqual(lean_resolve.numbers_in(text), lean_tables.printed_numbers(text)[0])

    def test_historical_comparison_remains_explicit(self):
        self.assertEqual(number_words.legacy_round6_numbers_in("four one-pair acts"), {5})
        self.assertEqual(number_words.numbers_in("four one-pair acts"), {4, 1})


if __name__ == "__main__":
    unittest.main()
