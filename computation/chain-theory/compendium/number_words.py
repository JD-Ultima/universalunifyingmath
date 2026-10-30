"""Integer words used in the compendium's printed-number audits.

Both audits use the same cardinal/ordinal parser. The legacy parser is retained
only to reproduce round 7's historical comparison with the original round 6.
"""
import re


UNITS = {w: i for i, w in enumerate("zero one two three four five six seven eight nine ten eleven twelve thirteen "
                                    "fourteen fifteen sixteen seventeen eighteen nineteen".split())}
TENS = {w: 10 * (i + 2) for i, w in enumerate("twenty thirty forty fifty sixty seventy eighty ninety".split())}
ORD = {w: i + 1 for i, w in enumerate("first second third fourth fifth sixth seventh eighth ninth tenth eleventh twelfth "
                                      "thirteenth fourteenth fifteenth sixteenth seventeenth eighteenth nineteenth".split())}
ORD.update({"twentieth": 20, "thirtieth": 30, "fortieth": 40, "fiftieth": 50, "sixtieth": 60, "seventieth": 70,
            "eightieth": 80, "ninetieth": 90, "hundredth": 100})


def cardinals(text):
    """Read integer digits and cardinal words, separating adjacent number phrases.

    Punctuation ends a phrase ('four, thirty-two' is 4 and 32). A unit/teen
    after a unit/teen also starts a new one ('four one-pair' is 4 and 1).
    """
    t = re.sub(r"\\label\{[^}]*\}", " ", text)
    t = re.sub(r"\\[a-zA-Z]+", " ", t)
    nums = {int(x) for x in re.findall(r"(?<![\w.])(\d+)(?![\w.]*\d)", t)}
    for seg in re.split(r"[,;:.()$\[\]{}]", t.replace("-", " ")):
        words = re.findall(r"[A-Za-z]+", seg)
        total = cur = None
        last = None

        def flush():
            nonlocal total, cur
            if cur is not None or total is not None:
                nums.add((total or 0) + (cur or 0))
            total = cur = None

        for w in words:
            w = w.lower()
            if w in UNITS:
                if last in ("unit", "teen"):
                    flush()
                cur = (cur or 0) + UNITS[w]
                last = "teen" if UNITS[w] >= 10 else "unit"
            elif w in TENS:
                if last in ("unit", "teen", "tens"):
                    flush()
                cur = (cur or 0) + TENS[w]
                last = "tens"
            elif w == "hundred" and cur is not None and last in ("unit", "teen"):
                cur *= 100
                last = "hundred"
            elif w == "thousand" and cur is not None:
                total = (total or 0) + cur * 1000
                cur = None
                last = "thousand"
            elif w == "and" and last in ("hundred", "thousand"):
                continue
            else:
                flush()
                last = None
        flush()
    return nums


def printed_numbers(text):
    """Return (all integers, integers added only by ordinal words).

    Includes simple ordinals and compound tens ('twenty-fourth' -> 24).
    """
    t = re.sub(r"\\label\{[^}]*\}", " ", text)
    ords = set()

    def comp(m):
        v = TENS[m.group(1).lower()] + ORD[m.group(2).lower()]
        ords.add(v)
        return f" {v} "

    t = re.sub(r"\b(twenty|thirty|forty|fifty|sixty|seventy|eighty|ninety)[- ](first|second|third|fourth|fifth|sixth|"
               r"seventh|eighth|ninth)\b", comp, t, flags=re.I)

    def simple(m):
        v = ORD[m.group(1).lower()]
        ords.add(v)
        return f" {v} "

    t = re.sub(r"\b(" + "|".join(sorted(ORD, key=len, reverse=True)) + r")\b", simple, t, flags=re.I)
    base = cardinals(text)
    return cardinals(t), ords - base


def numbers_in(text):
    """Shared public parser for both Lean coverage audits."""
    return printed_numbers(text)[0]


def legacy_round6_numbers_in(text):
    """Original round-6 parser, for the historical round-7 comparison only."""
    t = re.sub(r"\\label\{[^}]*\}", " ", text)
    t = re.sub(r"\\[a-zA-Z]+", " ", t)
    nums = {int(x) for x in re.findall(r"(?<![\w.])(\d+)(?![\w.]*\d)", t)}
    words = re.findall(r"[A-Za-z]+", t.replace("-", " "))
    i = 0
    while i < len(words):
        w = words[i].lower()
        if w in UNITS or w in TENS:
            total, cur, j = 0, 0, i
            while j < len(words):
                v = words[j].lower()
                if v in UNITS:
                    cur += UNITS[v]
                elif v in TENS:
                    cur += TENS[v]
                elif v == "hundred" and cur:
                    cur *= 100
                elif v == "thousand" and cur:
                    total += cur * 1000
                    cur = 0
                elif v == "and" and j + 1 < len(words) and (words[j + 1].lower() in UNITS or words[j + 1].lower() in TENS) and cur >= 100:
                    pass
                else:
                    break
                j += 1
            nums.add(total + cur)
            i = j
        else:
            i += 1
    return nums
