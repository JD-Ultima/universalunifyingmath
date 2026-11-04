"""Complete exact derivative certificate on the positive 14-dimensional hook.

No repository imports. Each input is one hook basis tensor and one derivative
covector; linearity therefore covers all first jets at a normal-frame point.
"""
from itertools import combinations
from time import perf_counter

N = 14


def product(left, right):
    parity = sum((right & ((1 << j) - 1)).bit_count()
                 for j in range(N) if left & (1 << j))
    return left ^ right, (-1 if parity % 2 else 1)


def insert(result, key, value):
    result[key] = result.get(key, 0) + value
    if not result[key]:
        del result[key]


def tensor(entries):
    out = [{} for _ in range(N)]
    for a, b, c, value in entries:
        assert b != c
        insert(out[a], (1 << b) | (1 << c), value if b < c else -value)
    return out


def basis():
    for a, b, c in combinations(range(N), 3):
        yield tensor([(a, b, c, 1), (b, a, c, 1)])
        yield tensor([(b, a, c, 1), (c, a, b, 1)])
    for b in range(N):
        rows = [a for a in range(N) if a != b]
        reference = rows[0]
        for a in rows[1:]:
            yield tensor([(a, a, b, 1), (reference, reference, b, -1)])


def hook_traces(B):
    alt, ctr = {}, {}
    for a, values in enumerate(B):
        for mask, value in values.items():
            key, sign = product(1 << a, mask)
            insert(ctr if mask & (1 << a) else alt, key, value * sign)
    return alt, ctr


def check(B, xi):
    actual, expected = [{} for _ in range(N)], [{} for _ in range(N)]
    for a in range(N):
        # Direct alternating-three-form commutator on H_bc. All nonzero
        # derivative pairs have one index equal to xi.
        for c in range(N):
            if c == xi or a in (xi, c):
                continue
            b0, c0 = sorted((xi, c))
            triple, sign = product(1 << a, 1 << b0)
            triple, extra = product(triple, 1 << c0)
            sign *= extra * (1 if xi < c else -1)
            for mask, coefficient in B[c].items():
                key, lr = product(triple, mask)
                other, rl = product(mask, triple)
                assert other == key
                insert(actual[a], key, coefficient * sign * (lr - rl))
        # Reduced formula 2(xi wedge B_a - e_a wedge B_xi).
        for mask, coefficient in B[a].items():
            if not mask & (1 << xi):
                key, sign = product(1 << xi, mask)
                insert(expected[a], key, 2 * coefficient * sign)
        for mask, coefficient in B[xi].items():
            if not mask & (1 << a):
                key, sign = product(1 << a, mask)
                insert(expected[a], key, -2 * coefficient * sign)
    assert actual == expected
    contracted = {}
    for a, values in enumerate(actual):
        for mask, coefficient in values.items():
            if mask & (1 << a):
                key, sign = product(1 << a, mask)
                insert(contracted, key, coefficient * sign)
    assert contracted == {key: -22 * value for key, value in B[xi].items()}


def main():
    start = perf_counter()
    tensors = list(basis())
    assert len(tensors) == 896
    for B in tensors:
        assert hook_traces(B) == ({}, {})
        for xi in range(N):
            check(B, xi)
    # The disjoint triple and repeated-index supports establish linear
    # independence and exhaust all 1274 - 364 - 14 hook dimensions.
    assert 2 * 364 + 14 * 12 == 14 * 91 - 364 - 14 == 896
    print('896 explicit integer basis tensors satisfy both hook identities.')
    print('12544/12544 derivative columns match the full commutator formula.')
    print('All contracted columns give -22 times the external divergence.')
    print('This certifies a linear identity on all first hook jets, not fields.')
    print(f'Runtime {perf_counter() - start:.3f} seconds.')


if __name__ == '__main__':
    main()
