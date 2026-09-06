#!/usr/bin/env python3
"""Cross-check nefvol against lrs on every convex piece of an .ine union.

Optional tool: the test suite uses the recorded values in tests/lrs_golden.txt
(regenerate with `lrs_check.py --golden inputs/*.ine`); lrs itself is not needed.

For each begin…end block: run lrs to get the vertices (H → V), run
`lrs … volume` on that V-representation, and compare with `nefvol volume`.
"""
import subprocess, sys, re, os, tempfile
from fractions import Fraction as F

LRS = os.environ.get("LRS", "/tmp/lrs/lrslib-073/lrsmp")
NEFVOL = os.environ.get("NEFVOL", "./target/release/nefvol")

def blocks(text):
    lines = [l.strip() for l in text.splitlines() if l.strip() and not l.startswith("*")]
    i = 0
    out = []
    while i < len(lines):
        if lines[i].lower() == "begin":
            m, n = lines[i+1].split()[:2]
            m, n = int(m), int(n)
            nums = []
            j = i + 2
            while len(nums) < m*n:
                nums += [F(x) for x in lines[j].split()]
                j += 1
            assert lines[j].lower() == "end"
            rows = [nums[r*n:(r+1)*n] for r in range(m)]
            out.append(rows)
            i = j + 1
        else:
            i += 1
    return out

def ine(rows, extra=""):
    n = len(rows[0])
    s = f"p\nH-representation\nbegin\n{len(rows)} {n} rational\n"
    for r in rows:
        s += " ".join(str(x) for x in r) + "\n"
    return s + "end\n" + extra

def run_lrs(text):
    with tempfile.NamedTemporaryFile("w", suffix=".ine", delete=False) as f:
        f.write(text); name = f.name
    out = subprocess.run([LRS, name], capture_output=True, text=True).stdout
    os.unlink(name)
    return out

def vertices(rows):
    out = run_lrs(ine(rows))
    vs = []
    for l in out.splitlines():
        t = l.split()
        if t and t[0] == "1" and len(t) == len(rows[0]):
            vs.append([F(x) for x in t[1:]])
    return vs

def lrs_volume(rows):
    # H → V with lrs, then `volume` on the V-representation (for H-representations
    # lrs reports the volume of the polar instead)
    vs = vertices(rows)
    n = len(rows[0])
    s = f"p\nV-representation\nbegin\n{len(vs)} {n} rational\n"
    for v in vs:
        s += "1 " + " ".join(str(x) for x in v) + "\n"
    s += "end\nvolume\n"
    out = run_lrs(s)
    m = re.search(r"olume=\s*([-0-9/]+)", out)
    return F(m.group(1)) if m else None

def main():
    if sys.argv[1] == "--golden":
        # print `instance piece volume` lines for tests/lrs_golden.txt
        for path in sys.argv[2:]:
            name = os.path.basename(path).rsplit(".", 1)[0]
            for k, rows in enumerate(blocks(open(path).read())):
                print(name, k, lrs_volume(rows))
        return
    path = sys.argv[1]
    bl = blocks(open(path).read())
    ok = True
    for k, rows in enumerate(bl):
        lv = lrs_volume(rows)
        text = ine(rows)
        with tempfile.NamedTemporaryFile("w", suffix=".ine", delete=False) as f:
            f.write(text); name = f.name
        nv = subprocess.run([NEFVOL, "volume", "--input", name], capture_output=True, text=True).stdout.strip()
        os.unlink(name)
        nv = F(nv) if nv else None
        status = "OK" if lv == nv else "MISMATCH"
        if lv != nv: ok = False
        print(f"piece {k}: lrs={lv} nefvol={nv} {status}")
    print("ALL OK" if ok else "FAILURES")
    sys.exit(0 if ok else 1)

main()
