#!/usr/bin/env python3
"""
SAT attack on Erdos Problem #273.

FORMALIZATION (cross-checked against https://www.erdosproblems.com/273 and
formal-conjectures/FormalConjectures/ErdosProblems/273.lean, StrictCoveringSystem):

  Question: does there exist a finite covering system of Z -- congruences
  a_i (mod n_i), i=1..k -- such that:
    (1) every integer is congruent to some a_i mod n_i ("covering"),
    (2) the moduli n_i are pairwise DISTINCT ("strict" covering system),
    (3) each n_i = p-1 for some prime p >= 5?
  (The p>=3 variant is a separate, already-solved sub-theorem using Selfridge's
  divisors-of-360 example; NOT attacked here since p=3 gives modulus 2, excluded.)

ENCODING. Fix a finite pool M of candidate moduli (all p-1 <= B, p prime, p>=5).
Since a congruence "a mod n" only makes sense relative to a fixed period L that
n divides, restrict to a sub-pool S subset M with lcm(S) = L <= subset-lcm-cap
(greedy: add moduli in increasing order while lcm stays under cap -- this is
exactly the "highly composite L, moduli dividing it" shape of known covering
systems, e.g. Selfridge's divisors of 360).

Boolean variable x_{n,a} for n in S, a in 0..n-1: "congruence a (mod n) is used".
  - AtMostOne(x_{n,*}) for each n in S  <-- encodes "distinct moduli": each n
    is used with at most one residue (n unused <=> all x_{n,*} false).
  - For each residue r in 0..L-1: OR_{n in S} x_{n, r mod n}  <-- "r is covered".
Any satisfying model directly gives a covering system with moduli subset of S
(all of the required form p-1, p>=5, all distinct). Coverage mod L implies
coverage of Z since every n in S divides L.

VERIFICATION. The model is NEVER trusted directly. After extracting the
selected (n,a) pairs from the SAT model, an independent check re-derives
L_verify = lcm(selected n's) and marks residues 0..L_verify-1 covered by
literally walking every residue class of every selected congruence (bytearray
sieve, not reusing solver internals) -- this is the "honest cycle" check.
"""
import argparse
import json
import time
from math import lcm
from pathlib import Path

from pysat.solvers import Cadical195

OUT_DIR = Path(__file__).resolve().parent


# ---------- pool construction ----------

def primes_up_to(n: int) -> list[int]:
    if n < 2:
        return []
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
    return [i for i, v in enumerate(sieve) if v]


def erdos_pool(B: int, variant: str = "full") -> list[int]:
    """variant=full: n = p-1 <= B (the actual problem). variant=half: n =
    (p-1)/2 <= B -- a heuristic structural lever, NOT the problem itself
    (see main() note on why a half-solution doesn't lift automatically)."""
    if variant == "full":
        ps = primes_up_to(B + 2)
        return sorted(p - 1 for p in ps if p >= 5 and p - 1 <= B)
    ps = primes_up_to(2 * B + 2)
    return sorted((p - 1) // 2 for p in ps if p >= 5 and (p - 1) // 2 <= B)


def reduce_pool_by_lcm(mods: list[int], cap: int) -> list[int]:
    """Greedily keep moduli (ascending) while lcm(kept) stays <= cap."""
    L = 1
    kept = []
    for n in mods:
        newL = lcm(L, n)
        if newL <= cap:
            L = newL
            kept.append(n)
    return kept


# ---------- CNF construction ----------

def build_and_add_clauses(solver: Cadical195, mods: list[int], L: int) -> tuple[dict, int, int]:
    """Add AtMostOne (per modulus) + coverage (per residue mod L) clauses
    directly to `solver`. Returns (var_map, n_vars, n_clauses)."""
    var = {}
    counter = 1
    for n in mods:
        for a in range(n):
            var[(n, a)] = counter
            counter += 1
    n_vars = counter - 1

    n_clauses = 0
    # AtMostOne per modulus (pairwise)
    for n in mods:
        vs = [var[(n, a)] for a in range(n)]
        for i in range(len(vs)):
            for j in range(i + 1, len(vs)):
                solver.add_clause([-vs[i], -vs[j]])
                n_clauses += 1

    # Coverage: build clause literal lists residue-by-residue, filled
    # modulus-by-modulus using stride slicing (fast).
    clause_lits: list[list[int]] = [[] for _ in range(L)]
    for n in mods:
        for a in range(n):
            v = var[(n, a)]
            for r in range(a, L, n):
                clause_lits[r].append(v)
    for lits in clause_lits:
        solver.add_clause(lits)
        n_clauses += 1

    return var, n_vars, n_clauses


# ---------- solve ----------
# NB: this pysat build's Cadical195.interrupt()/solve_limited() raises
# NotImplementedError ("Limited solve is currently unsupported by CaDiCaL"),
# so an in-process soft timeout is not available. The 300s solo-call timeout
# is instead enforced externally (OS-level kill of the whole process, e.g.
# the `timeout` param of the invoking shell) -- if killed, no result file is
# written and the caller records TIMEOUT with the observed wall-clock cutoff.

def solve(solver: Cadical195):
    t0 = time.time()
    status = solver.solve()
    elapsed = time.time() - t0
    return ("SAT" if status else "UNSAT"), elapsed


# ---------- independent verification ----------

def verify_covering(system: list[tuple[int, int]]) -> tuple[bool, int]:
    """Honest re-check: does the given set of (n,a) congruences cover every
    residue mod lcm(all n)? Returns (covered_fully, L_verify)."""
    if not system:
        return False, 1
    L_verify = 1
    for n, _ in system:
        L_verify = lcm(L_verify, n)
    covered = bytearray(L_verify)
    for n, a in system:
        covered[a::n] = bytearray(b"\x01") * len(covered[a::n])
    return covered.count(0) == 0, L_verify


# ---------- one experiment ----------

def run_case(name: str, mods: list[int], pool_size: int | None = None, partial: bool = False) -> dict:
    L = 1
    for n in mods:
        L = lcm(L, n)
    print(f"[{name}] pool_used={len(mods)} moduli={mods} L={L}", flush=True)

    solver = Cadical195()
    t_enc0 = time.time()
    var, n_vars, n_clauses = build_and_add_clauses(solver, mods, L)
    t_enc = time.time() - t_enc0
    print(f"[{name}] encoded: n_vars={n_vars} n_clauses={n_clauses} encode_time={t_enc:.2f}s", flush=True)

    status, solve_time = solve(solver)
    print(f"[{name}] solve status={status} time={solve_time:.2f}s", flush=True)

    result = {
        "name": name,
        "pool_size": pool_size if pool_size is not None else len(mods),
        "partial": partial,  # True <=> S is a strict subset of the full pool
                              # (subset-lcm-cap truncated) -- UNSAT on a partial
                              # pool is NOT a full bound for that B.
        "moduli_used": mods,
        "L": L,
        "n_vars": n_vars,
        "n_clauses": n_clauses,
        "encode_time_sec": round(t_enc, 3),
        "status": status,
        "solve_time_sec": round(solve_time, 3),
        "solution": None,
        "verified": None,
        "L_verify": None,
    }

    if status == "SAT":
        model = set(solver.get_model())
        system = []
        for n in mods:
            for a in range(n):
                if var[(n, a)] in model:
                    system.append((n, a))
        distinct_moduli = len({n for n, _ in system})
        assert distinct_moduli == len(system), "AtMostOne violated -- encoding bug"
        ok, L_verify = verify_covering(system)
        result["solution"] = system
        result["verified"] = ok
        result["L_verify"] = L_verify
        print(f"[{name}] solution={system} verified={ok} (L_verify={L_verify})", flush=True)

    solver.delete()
    return result


# ---------- CLI ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["sanity", "erdos"], required=True)
    ap.add_argument("--B", type=int, default=30)
    ap.add_argument(
        "--variant",
        choices=["full", "half"],
        default="full",
        help="full: n=p-1 (the problem). half: n=(p-1)/2 -- heuristic lever, "
        "a SAT solution here does NOT solve the full problem (one congruence "
        "a mod k splits into TWO a mod 2k / a+k mod 2k, but we allow only one "
        "congruence per modulus -- see docstring).",
    )
    ap.add_argument("--subset-lcm-cap", type=int, default=5_000_000)
    ap.add_argument(
        "--timeout",
        type=float,
        default=300.0,
        help="nominal solo-call budget in sec (NOT enforced in-process -- "
        "wrap this script with an OS-level timeout, see module docstring)",
    )
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()
    print(f"nominal timeout budget={args.timeout}s (enforce externally)", flush=True)

    if args.mode == "sanity":
        mods = [2, 3, 4, 6, 12]
        result = run_case("sanity_toy", mods, pool_size=len(mods))
    else:
        full_pool = erdos_pool(args.B, args.variant)
        mods = reduce_pool_by_lcm(full_pool, args.subset_lcm_cap)
        partial = len(mods) < len(full_pool)
        print(f"[B={args.B} variant={args.variant}] full_pool ({len(full_pool)})={full_pool}", flush=True)
        if partial:
            print(
                f"[B={args.B}] subset-lcm-cap={args.subset_lcm_cap} truncated pool "
                f"to {len(mods)}/{len(full_pool)} moduli -- result is PARTIAL, "
                f"not a full bound for this B",
                flush=True,
            )
        result = run_case(f"erdos_B{args.B}_{args.variant}", mods, pool_size=len(full_pool), partial=partial)
        result["B"] = args.B
        result["variant"] = args.variant
        result["subset_lcm_cap"] = args.subset_lcm_cap
        result["full_pool"] = full_pool

    out_path = Path(args.out) if args.out else OUT_DIR / f"result_{result['name']}.json"
    out_path.write_text(json.dumps(result, indent=2))
    print(f"written {out_path}", flush=True)


if __name__ == "__main__":
    main()
