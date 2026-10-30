"""TH-001 section 5: H1/H2 realised in semirings whose addition is NOT ordinary +.

log2-semiring: a (+) b = log2(2^a + 2^b), a (x) b = a + b   (iso to (R_{>0}, +, *) via x -> 2^x)
tropical:      a (+) b = max(a, b),        a (x) b = a + b
Exact checks with sympy; distributivity checked symbolically.
"""
import sympy as sp

a, b, c = sp.symbols('a b c', real=True)
oplus = lambda p, q: sp.log(2**p + 2**q, 2)
otimes = lambda p, q: p + q

print("sympy", sp.__version__)
print("log2-semiring:")
print("  H1  1(x)1 =", otimes(1, 1), "; H2  2(x)2 =", otimes(2, 2))
print("  1(+)1 =", sp.simplify(oplus(1, 1)), "; 2(+)2 =", sp.simplify(oplus(2, 2)))
nbar = lambda n: sp.simplify(sp.log(n * 2**1, 2))        # n-fold (+)-sum of 1 = log2(n * 2)
print("  numerals nbar = n-fold (+)-sum of 1 = 1 + log2 n:", {n: nbar(n) for n in (1, 2, 4, 8)})
print("  2bar (x) 2bar =", sp.simplify(otimes(nbar(2), nbar(2))), "= 8bar:", sp.simplify(otimes(nbar(2), nbar(2)) - nbar(8)) == 0,
      "; 4bar =", nbar(4))
print("  (x) distributes over (+):", sp.simplify(2**(otimes(a, oplus(b, c))) - 2**(oplus(otimes(a, b), otimes(a, c)))) == 0)
print("  iso x -> 2^x: 2^(a(x)b) = 2^a*2^b:", sp.simplify(2**otimes(a, b) - 2**a * 2**b) == 0,
      "; 2^(a(+)b) = 2^a + 2^b:", sp.simplify(2**oplus(a, b) - (2**a + 2**b)) == 0)
print("tropical:")
print("  1 max 1 =", sp.Max(1, 1), "; distributivity a+max(b,c) = max(a+b,a+c) on a grid:",
      all(x + max(y, z) == max(x + y, x + z) for x in range(-3, 4) for y in range(-3, 4) for z in range(-3, 4)))
