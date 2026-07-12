#!/usr/bin/env python3
"""
Export the Erdos #273 encoding (see encode.py docstring) as a DIMACS CNF file,
for a standalone proof-logging solver run (cadical --lrat).

Variable numbering and clause order are IDENTICAL to encode.py
build_and_add_clauses(): vars enumerated modulus-ascending then residue
0..n-1; clauses = pairwise AtMostOne per modulus (ascending), then one
coverage clause per residue r = 0..L-1.

A sidecar .meta.json records the pool, L, and the var map derivation so the
CNF is reproducible and auditable without re-reading this script.
"""
import argparse
import hashlib
import json
from math import lcm
from pathlib import Path

from encode import erdos_pool, reduce_pool_by_lcm

OUT_DIR = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, required=True)
    ap.add_argument("--subset-lcm-cap", type=int, default=5_000_000)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    full_pool = erdos_pool(args.B, "full")
    mods = reduce_pool_by_lcm(full_pool, args.subset_lcm_cap)
    if mods != full_pool:
        raise SystemExit(
            f"pool truncated by lcm cap ({len(mods)}/{len(full_pool)}) -- "
            f"a certificate for a PARTIAL pool proves nothing; aborting"
        )

    L = 1
    for n in mods:
        L = lcm(L, n)

    var = {}
    counter = 1
    for n in mods:
        for a in range(n):
            var[(n, a)] = counter
            counter += 1
    n_vars = counter - 1

    # coverage clause for residue r is exactly one literal per modulus:
    # var(n, r mod n) -- same clause content/order as encode.py's stride fill
    # (per-r lits appended modulus-ascending), but built streaming so L in the
    # tens of millions doesn't hold 14*L ints in RAM.
    n_clauses = sum(n * (n - 1) // 2 for n in mods) + L
    out_path = Path(args.out) if args.out else OUT_DIR / f"erdos273_B{args.B}.cnf"
    with open(out_path, "w") as f:
        f.write(f"c Erdos 273: covering system, distinct moduli n=p-1<=B, p prime>=5\n")
        f.write(f"c B={args.B} moduli={mods} L={L}\n")
        f.write(f"p cnf {n_vars} {n_clauses}\n")
        for n in mods:
            vs = [var[(n, a)] for a in range(n)]
            for i in range(len(vs)):
                for j in range(i + 1, len(vs)):
                    f.write(f"-{vs[i]} -{vs[j]} 0\n")
        base = {n: var[(n, 0)] for n in mods}
        buf = []
        for r in range(L):
            buf.append(" ".join(str(base[n] + r % n) for n in mods) + " 0")
            if len(buf) >= 100_000:
                f.write("\n".join(buf) + "\n")
                buf.clear()
        if buf:
            f.write("\n".join(buf) + "\n")

    sha = hashlib.sha256(out_path.read_bytes()).hexdigest()
    meta = {
        "B": args.B,
        "moduli": mods,
        "L": L,
        "n_vars": n_vars,
        "n_clauses": n_clauses,
        "var_map": "var(n,a) = 1 + sum(n' for n' in moduli if n' < n) + a",
        "clause_order": "pairwise AtMostOne per modulus (ascending), then coverage clause per residue 0..L-1",
        "cnf_sha256": sha,
    }
    meta_path = out_path.with_suffix(".meta.json")
    meta_path.write_text(json.dumps(meta, indent=2))
    print(f"written {out_path} ({out_path.stat().st_size/1e6:.1f} MB) sha256={sha[:16]}…")
    print(f"n_vars={n_vars} n_clauses={n_clauses} L={L}")
    print(f"meta {meta_path}")


if __name__ == "__main__":
    main()
