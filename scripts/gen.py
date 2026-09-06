#!/usr/bin/env python3
"""Generate benchmark unions as .ine files.

  gen.py boxes  D K SEED   -> K random integer boxes in R^D under random unimodular shears
  gen.py grid   D K        -> K^D unit cubes tiling [0,K]^D shifted (many coincident facets)
  gen.py cones  D M K SEED -> K random polytopes each with M facets around a common vertex (degenerate)
  gen.py rand   D M K SEED -> K random polytopes: box ∩ M random halfspaces
"""
import sys, random

def emit(polys, d):
    out = []
    for i, cons in enumerate(polys):
        out.append(f"poly_{i}\nH-representation\nbegin\n{len(cons)} {d+1} integer")
        for a, b in cons:
            out.append(" ".join(str(x) for x in [b] + [-x for x in a]))
        out.append("end")
    print("\n".join(out))

def box(lo, hi):
    d = len(lo)
    cons = []
    for i in range(d):
        a = [0]*d; a[i] = 1; cons.append((a, hi[i]))
        a = [0]*d; a[i] = -1; cons.append((a, -lo[i]))
    return cons

def shear(cons, d, rng, k=2):
    # apply x -> M x with M unimodular: constraint a·x ≤ b becomes (a M^{-1})·y ≤ b
    # build M^{-1} as a product of elementary shears
    minv = [[1 if i == j else 0 for j in range(d)] for i in range(d)]
    for _ in range(k):
        i, j = rng.sample(range(d), 2)
        s = rng.randint(-2, 2)
        # minv <- minv * (I + s E_ij)
        for r in range(d):
            minv[r][j] += s * minv[r][i]
    out = []
    for a, b in cons:
        a2 = [sum(a[i]*minv[i][j] for i in range(d)) for j in range(d)]
        out.append((a2, b))
    return out

def main():
    kind = sys.argv[1]
    if kind == "boxes":
        d, k, seed = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
        rng = random.Random(seed)
        polys = []
        for _ in range(k):
            lo = [rng.randint(-6, 4) for _ in range(d)]
            hi = [l + rng.randint(2, 6) for l in lo]
            polys.append(shear(box(lo, hi), d, rng))
        emit(polys, d)
    elif kind == "grid":
        d, k = int(sys.argv[2]), int(sys.argv[3])
        import itertools
        polys = []
        for idx in itertools.product(range(k), repeat=d):
            lo = list(idx); hi = [x+1 for x in idx]
            polys.append(box(lo, hi))
        emit(polys, d)
    elif kind == "cones":
        d, m, k, seed = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        rng = random.Random(seed)
        polys = []
        for _ in range(k):
            # M facets through the origin-ish apex: a·x ≤ 0 with random a having positive last coord, plus cap x_d ≤ 1 ... apex at origin
            cons = []
            for _ in range(m):
                a = [rng.randint(-9, 9) for _ in range(d-1)] + [-rng.randint(1, 9)]
                cons.append((a, 0))
            cap = [0]*(d-1) + [1]
            cons.append((cap, rng.randint(1, 3)))
            polys.append(shear(cons, d, rng, 1))
        emit(polys, d)
    elif kind == "rand":
        d, m, k, seed = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        rng = random.Random(seed)
        polys = []
        for _ in range(k):
            lo = [rng.randint(-6, 4) for _ in range(d)]
            hi = [l + rng.randint(3, 8) for l in lo]
            cons = box(lo, hi)
            c = [(l+h)/2 for l, h in zip(lo, hi)]
            for _ in range(m):
                a = [rng.randint(-5, 5) for _ in range(d)]
                if all(x == 0 for x in a): continue
                # keep the centre feasible with margin
                b = int(sum(x*y for x, y in zip(a, c))) + rng.randint(1, 8)
                cons.append((a, b))
            polys.append(cons)
        emit(polys, d)
    else:
        print(__doc__); sys.exit(1)

main()
