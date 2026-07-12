# Computational bounds for three Erdős problems

Verified computational results for [Erdős problems](https://www.erdosproblems.com)
**#273**, **#385**, and **#647**. All three problems remain open — these are
machine-checked partial bounds, not solutions.

| Problem | Result | Verification |
|---|---|---|
| [#273](https://www.erdosproblems.com/273) | No covering system of ℤ with pairwise distinct moduli of the form p−1 ≤ 50 (p prime ≥ 5) | UNSAT certificate (LRAT) checked by 3 independent checkers incl. the formally verified [cake_lpr](https://github.com/tanyongkiam/cake_lpr) |
| [#385](https://www.erdosproblems.com/385) | F(n) > n for all 10⁴ < n ≤ 1.0011×10¹² except the two known values n = 267672, 267680 | Segmented sieve, deterministic checkpoints, shard cross-checks |
| [#647](https://www.erdosproblems.com/647) | No n > 24 with max_{m<n}(m + τ(m)) ≤ n + 2 up to 1.0011×10¹² | Same sieve pass; known anchor R(24)=26 reproduced |

## #273 — SAT certificate (`sat273/`)

A covering system with all moduli from the full admissible pool
{4, 6, 10, 12, 16, 18, 22, 28, 30, 36, 40, 42, 46} (every p−1 ≤ 50, p ≥ 5 prime;
lcm L = 1,275,120) exists iff a CNF is satisfiable: one variable per congruence
class (n, a), pairwise at-most-one per modulus (distinct moduli), one coverage
clause per residue mod L. The CNF (310 vars, 1,279,875 clauses) is UNSAT.

- `export_dimacs.py` — deterministic CNF generator (~40 lines of logic);
  SHA-256 of the generated CNF is pinned in `erdos273_B50.meta.json`
- `artifacts/` — CNF (1.7 MB xz) and LRAT certificate (11 MB xz)
- See `verify.md` for end-to-end verification (≈ 15 minutes from scratch)

Certificate produced by CaDiCaL 3.0.0 with `--plain` (inprocessing disabled —
without it CaDiCaL emits extension variables that one popular checker silently
accepts; see `verify.md`, "Why --plain").

## #385 / #647 — segmented sieve (`sieve647/`)

One pass computes τ(m) and the smallest prime factor spf(m) block-by-block
(numpy, block 10⁷) and maintains both running records:
R(n) = max_{m<n}(m + τ(m)) for #647 and F(n) = max_{m<n, m composite}(m + p(m))
for #385. Range [1.1×10⁹, 1.0011×10¹²] ran as 8 shards with 10⁶ re-warm overlap
(records are monotone; boundary findings would be flagged — none occurred);
[1, 1.1×10⁹] ran as a single process earlier. ~18 core-hours on Apple M3.

- `sieve.py` — sieve + both detectors + checkpointing (formalizations quoted
  from the problem statements in the module docstring)
- `checkpoints/` — final checkpoints of all shards (`next_lo` == shard end)
- `findings.jsonl` — the only findings: the two known #385 values

For #385 this extends the verified range (previously ~10⁸, see C.K.S.'s
comments on Tao's blog) by ×10⁴ and independently reproduces both known
exceptional values. Per-shard minima of F(n)−n grow monotonically
(23,989 → 868,399), consistent with F(n)−n → ∞.

## License

Apache-2.0
