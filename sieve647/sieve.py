"""
Combined segmented sieve for Erdos problems #647 and #385.
One pass over blocks m=1..TARGET, two independent detectors sharing the
tau(m)/spf(m) arrays.

=== Formalization of #647 (Erdos-Selfridge) ===
Exact statement: let tau(n) be the number of divisors of n. Does there exist
n>24 such that
    max_{m<n} (m + tau(m)) <= n+2 ?
Known true for n=24. The bound n+2 is tight: max(tau(n-1)+n-1,
tau(n-2)+n-2) >= n+2 for all n (proved in the literature, not re-checked
here). Erdos called it "extremely doubtful" that infinitely many such n
exist, and offered $44 for a SINGLE example (the problem is open).

Detector: R(n) = max_{m<n} (m+tau(m)) -- a running (non-decreasing) record.
Search for n>24 with R(n) <= n+2. A hit would answer the open problem.

=== Formalization of #385 (Erdos, Eggleton, Selfridge) ===
Exact statement: F(n) = max_{m<n, m composite} (m + p(m)), where p(m) is the
smallest prime factor of m. Questions: (1) is F(n) > n for all sufficiently
large n? (2) does F(n)-n tend to infinity?
Trivially F(n) <= n + sqrt(n) (since p(m) <= sqrt(m) for composite m).

Detector: F(n) is a running record taken ONLY over composite m (primes and
m=1 are excluded from the maximum). Search for n with F(n) <= n -- a
violation of conjecture (1). Small n are ignored as hits (see SMALL_N_385):
the statement is explicitly about "all sufficiently large n", so isolated
violations at small n are expected. In parallel, track the running min of
(F(n)-n) over large n as an indicator for (2).

=== Notes on problem-statement precision ===
- #647's goal is sometimes loosely phrased as "find n disproving the
  statement". More precisely: the question is "does there exist n>24...";
  a found n ANSWERS the question affirmatively (a solution to the open
  problem, Erdos's $44 prize), not a disproof.
- Informal descriptions of #385 sometimes drop "large n"; the original
  statement explicitly restricts the question to "sufficiently large n" --
  implemented here via the SMALL_N_385 threshold.

=== Sharding (--lo/--hi) ===
The running records R(n)/F(n) are global (they depend on the ENTIRE history
m=1..n). An independent shard doesn't know the true record at its left
boundary, so a shard actually sieves starting at (lo-overlap) to "warm up"
the local record, but only logs findings from m=lo onward. Findings with
m-lo < overlap are flagged near_shard_boundary=true -- the local record
there may not have caught up with the true one yet (needs re-checking by
stitching with the neighboring shard).
"""
import os
import sys
import json
import math
import time
import argparse

import numpy as np

WORKDIR = os.path.dirname(os.path.abspath(__file__))
CKPT_PATH = os.path.join(WORKDIR, "checkpoint.json")
FINDINGS_PATH = os.path.join(WORKDIR, "findings.jsonl")

BLOCK = 10_000_000
PROGRESS_EVERY = 5
SMALL_N_385 = 10_000
NEG_INF = np.iinfo(np.int64).min // 2


def sieve_primes_upto(n):
    if n < 2:
        return np.array([], dtype=np.int64)
    is_p = np.ones(n + 1, dtype=bool)
    is_p[:2] = False
    for i in range(2, math.isqrt(n) + 1):
        if is_p[i]:
            is_p[i * i::i] = False
    return np.nonzero(is_p)[0].astype(np.int64)


def block_tau(lo, hi):
    """tau(m) for m in [lo,hi) via the paired-divisor scheme:
    tau(m) = 2*|{d<=sqrt(m): d|m}| - [m is a perfect square].
    counts is uint16 (tau(m)<~11000 for m<1e13, fits with margin);
    the final tau is uint16 too (saves memory on large blocks)."""
    n = hi - lo
    counts = np.zeros(n, dtype=np.uint16)
    D = math.isqrt(hi - 1)
    for d in range(1, D + 1):
        start = max(lo, d * d)
        rem = start % d
        if rem:
            start += d - rem
        if start >= hi:
            continue
        counts[start - lo: hi - lo: d] += 1
    tau = 2 * counts.astype(np.int32)
    s = math.isqrt(lo - 1) + 1
    sq = s * s
    while sq < hi:
        tau[sq - lo] -= 1
        s += 1
        sq = s * s
    return tau.astype(np.uint16)


def block_spf(lo, hi, primes):
    """Smallest prime factor of m for m in [lo,hi). 0 if m is prime or m=1.
    int32 -- values up to sqrt(1e13)~3.16e6, uint16 isn't enough."""
    n = hi - lo
    spf = np.zeros(n, dtype=np.int32)
    for p in primes:
        start = max(lo, int(p))
        rem = start % p
        if rem:
            start += p - rem
        if start >= hi:
            continue
        sub = spf[start - lo: hi - lo: p]
        sub[sub == 0] = p
    return spf


def load_ckpt(path, default_next_lo=1):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {"next_lo": default_next_lo, "max647": 0, "max385": 0, "min_gap385": None}


def save_ckpt(state, path):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(state, f)
    os.replace(tmp, path)


def log_finding(rec, path):
    with open(path, "a") as f:
        f.write(json.dumps(rec) + "\n")


def process_block(lo, hi, primes, max647, max385, min_gap385,
                   report_from=0, overlap=0, verbose=True, findings_path=FINDINGS_PATH):
    """One block: computes tau/spf, updates the running records (in-place, to
    avoid extra int64 copies), logs findings with m>=report_from.
    Returns (max647, max385, min_gap385, n_findings647, n_findings385)."""
    tau = block_tau(lo, hi)
    D_block = math.isqrt(hi - 1)
    idx = np.searchsorted(primes, D_block, side="right")
    spf = block_spf(lo, hi, primes[:idx])

    m = np.arange(lo, hi, dtype=np.int64)
    n_vals = m + 1
    report_mask = m >= report_from

    # --- #647 --- (running647 reuses the excess647 buffer: accumulate/maximum with out=)
    excess647 = m + tau
    np.maximum.accumulate(excess647, out=excess647)
    np.maximum(excess647, max647, out=excess647)
    hit647 = np.where((excess647 <= n_vals + 2) & (n_vals > 24) & report_mask)[0]
    for i in hit647:
        n_i = int(n_vals[i])
        rec = {"problem": 647, "n": n_i, "R_n": int(excess647[i]),
               "near_shard_boundary": bool(int(m[i]) - report_from < overlap)}
        log_finding(rec, findings_path)
        if verbose:
            print("FINDING 647:", rec)
    max647 = int(excess647[-1])

    # --- #385 ---
    composite = (spf > 0) & (spf < m)
    excess385 = np.where(composite, m + spf, NEG_INF)
    np.maximum.accumulate(excess385, out=excess385)
    np.maximum(excess385, max385, out=excess385)
    big = n_vals > SMALL_N_385
    hit385 = np.where((excess385 <= n_vals) & big & report_mask)[0]
    for i in hit385:
        n_i = int(n_vals[i])
        rec = {"problem": 385, "n": n_i, "F_n": int(excess385[i]),
               "near_shard_boundary": bool(int(m[i]) - report_from < overlap)}
        log_finding(rec, findings_path)
        if verbose:
            print("FINDING 385:", rec)
    max385 = int(excess385[-1])
    big_report = big & report_mask
    if big_report.any():
        gap_min = int((excess385[big_report] - n_vals[big_report]).min())
        if min_gap385 is None or gap_min < min_gap385:
            min_gap385 = gap_min

    return max647, max385, min_gap385, len(hit647), len(hit385)


def run(sieve_lo, sieve_hi, report_from, overlap, block, ckpt_path, findings_path,
        reset=False, verbose=True):
    if reset:
        for p in (ckpt_path, findings_path):
            if os.path.exists(p):
                os.remove(p)
    state = load_ckpt(ckpt_path, default_next_lo=sieve_lo)
    lo = state["next_lo"]
    lo0 = lo
    max647, max385, min_gap385 = state["max647"], state["max385"], state["min_gap385"]

    D_target = math.isqrt(max(sieve_hi - 1, 1))
    primes = sieve_primes_upto(D_target)

    t0 = time.time()
    blocks_done = 0
    while lo < sieve_hi:
        hi = min(lo + block, sieve_hi)
        max647, max385, min_gap385, nf647, nf385 = process_block(
            lo, hi, primes, max647, max385, min_gap385,
            report_from=report_from, overlap=overlap, verbose=verbose,
            findings_path=findings_path)
        lo = hi
        blocks_done += 1
        state = {"next_lo": lo, "max647": max647, "max385": max385, "min_gap385": min_gap385}
        save_ckpt(state, ckpt_path)

        if blocks_done % PROGRESS_EVERY == 0 or lo >= sieve_hi:
            elapsed = time.time() - t0
            speed = (lo - lo0) / elapsed if elapsed > 0 else 0
            print(f"[progress] lo={lo:,} elapsed={elapsed:.1f}s speed={speed:,.0f} n/s "
                  f"max647={max647} max385={max385} min_gap385={min_gap385} "
                  f"finds647={nf647} finds385={nf385}", flush=True)

    return state


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("target", nargs="?", type=int, default=None,
                     help="Sequential mode: sieve m=1..target (checkpoint.json/findings.jsonl)")
    ap.add_argument("--block", type=int, default=BLOCK, help="Block size, in numbers")
    ap.add_argument("--reset", action="store_true")
    ap.add_argument("--lo", type=int, default=None, help="Shard: m>=lo is reported (m-range, not n)")
    ap.add_argument("--hi", type=int, default=None, help="Shard: end of the m-range (exclusive)")
    ap.add_argument("--shard-overlap", type=int, default=1_000_000,
                     help="How far before lo to actually start sieving, to warm up the record")
    args = ap.parse_args()

    if args.lo is not None or args.hi is not None:
        if args.lo is None or args.hi is None:
            sys.exit("--lo and --hi are required together")
        sieve_lo = max(1, args.lo - args.shard_overlap)
        ckpt_path = os.path.join(WORKDIR, f"checkpoint_{args.lo}_{args.hi}.json")
        findings_path = os.path.join(WORKDIR, f"findings_{args.lo}_{args.hi}.jsonl")
        run(sieve_lo, args.hi, report_from=args.lo, overlap=args.shard_overlap,
            block=args.block, ckpt_path=ckpt_path, findings_path=findings_path,
            reset=args.reset)
    else:
        tgt = args.target if args.target is not None else 10 ** 9
        run(1, tgt, report_from=0, overlap=0, block=args.block,
            ckpt_path=CKPT_PATH, findings_path=FINDINGS_PATH, reset=args.reset)
