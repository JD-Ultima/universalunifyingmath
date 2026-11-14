"""Verify rank-four H4 scalarity and top-grade feasibility over Q(sqrt(22))."""

from fractions import Fraction
from itertools import combinations

import CD_clifford_hook_identities as CD


def main():
    external = range(4)
    light = range(4, 10)
    heavy = range(10, 14)
    volume = (1 << 4) - 1
    first = [{} for _ in range(14)]
    third = [{} for _ in range(14)]
    for t, (a, b) in zip(light, combinations(external, 2)):
        CD.BC.plus(third[a], (volume ^ (1 << b)) | (1 << t), (-1)**b)
        CD.BC.plus(third[b], (volume ^ (1 << a)) | (1 << t), (-1)**a)

    signs = [(1, 1, -1, -1), (1, -1, 1, -1), (1, -1, -1, 1)]
    light_axes = iter(light)
    for i, diagonal in enumerate(signs):
        form = {sum(1 << axis for j, axis in enumerate(heavy) if j != i): 2}
        for j in range(3):
            if j != i:
                mask = sum(1 << axis for k, axis in enumerate(heavy[:3]) if k != j)
                form[mask | (1 << next(light_axes))] = 1
        for a, sign in enumerate(diagonal):
            first[a] = CD.add(first[a], CD.scale(CD.wedge(a, form), sign))

    for tensor in (first, third):
        assert not CD.add(*(CD.wedge(a, blade) for a, blade in enumerate(tensor)))
        assert not CD.add(*(CD.interior(a, blade) for a, blade in enumerate(tensor)))
        assert all(all(mask >> a & 1 for mask in blade) for a, blade in enumerate(tensor))

    clifford = CD.BC.Clifford([1] * 14)
    q1, q3 = clifford.quadratic(first), clifford.quadratic(third)
    total = clifford.quadratic([CD.add(a, b) for a, b in zip(first, third)])
    rational = [CD.add(a, CD.scale(b, 22)) for a, b in zip(q1, q3)]
    mixed = [CD.add(c, CD.scale(a, -1), CD.scale(b, -1)) for c, a, b in zip(total, q1, q3)]
    norm = sum(value**2 for blade in first for value in blade.values())
    norm += 22 * sum(value**2 for blade in third for value in blade.values())
    assert norm == 336
    for a in range(14):
        assert CD.grade(rational[a], 1) == {1 << a: -Fraction(6, 7) * norm}
        assert not CD.grade(mixed[a], 1)
        assert not CD.grade(rational[a], 9) and not CD.grade(mixed[a], 9)
    assert sum(len(CD.grade(value, 5)) for value in rational) == 270
    assert sum(len(CD.grade(value, 5)) for value in mixed) == 228

    def inner(left, right):
        return sum(value * right.get(mask, 0) for mask, value in left.items())

    for a in range(14):
        for b in range(14):
            assert not inner(first[a], third[b])
            actual = inner(first[a], first[b]) + 22 * inner(third[a], third[b])
            assert actual == (84 if a == b and a in external else 0)
    print("Exact H4 rank-four hook verified. Complete scalarity and grade nine pass. Grade five is nonzero.")


if __name__ == "__main__":
    main()
