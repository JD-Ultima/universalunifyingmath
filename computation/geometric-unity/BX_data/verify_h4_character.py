"""Decompose the exact rational H4 dual-curvature weight character."""

from functools import lru_cache
from itertools import product
import json

from verify_h4_an import HERE


@lru_cache(None)
def kostant(a, b):
    if a < 0:
        return 0
    # n1(e1-e2)+n2 e2+n3 e1+n4(e1+e2)=(a,b)
    return sum(b + n1 - n4 >= 0 for n1 in range(a + 1) for n4 in range(a - n1 + 1))


def mult(highest, weight):
    x = (2 * highest[0] + 3, 2 * highest[1] + 1)
    y = (2 * weight[0] + 3, 2 * weight[1] + 1)
    value = 0
    for order in ((0, 1), (1, 0)):
        for s, t in product((-1, 1), repeat=2):
            a, b = s * x[order[0]] - y[0], t * x[order[1]] - y[1]
            if a % 2 == b % 2 == 0:
                value += ((1 if order == (0, 1) else -1) * s * t) * kostant(a // 2, b // 2)
    return value


def dimension(a, b):
    return (2 * a + 3) * (2 * b + 1) * (a + b + 2) * (a - b + 1) // 6


def main():
    receipt = json.loads((HERE / "h4_dual_curvature_manifest.json").read_text())
    remaining = {}
    for row in receipt["blocks"]:
        a, b = row["weight"]
        remaining[a, b] = remaining[-a, -b] = row["columns"] - row["exact_rank"]
    assert sum(remaining.values()) == 10388
    components = []
    while any(remaining.values()):
        a, b = max(w for w, c in remaining.items() if c)
        assert a >= b >= 0
        copies = remaining[a, b]
        dist = {w: mult((a, b), w) for w in remaining}
        assert all(x >= 0 for x in dist.values()) and sum(dist.values()) == dimension(a, b)
        for w, v in dist.items():
            remaining[w] -= copies * v
            assert remaining[w] >= 0
        components.append(dict(weight=[a, b], copies=copies, dimension=dimension(a, b), combined=copies * dimension(a, b)))
    assert sum(x["combined"] for x in components) == 10388
    summary = dict(dimension=10388, irreducible_types=len(components),
                   number_of_copies=sum(b["copies"] for b in components), components=components)
    recorded = json.loads((HERE / "h4_character.json").read_text(encoding="utf-8"))
    assert summary == recorded
    print(f"Verified {summary['number_of_copies']} copies of "
          f"{summary['irreducible_types']} SO(5) types, dimension {summary['dimension']}.")


if __name__ == "__main__":
    main()
