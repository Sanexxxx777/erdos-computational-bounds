# Computational bounds for three Erdős problems

*На русском: верифицированные вычислительные результаты для трёх открытых проблем Эрдёша (#273, #385, #647) — SAT-решение с проверяемым LRAT-сертификатом и сегментированное решето, независимо проверяемые несколькими чекерами. Полные детали и команды воспроизведения — ниже, на английском.*

Verified computational results for [Erdős problems](https://www.erdosproblems.com)
**#273**, **#385**, and **#647**. All three problems remain open — these are
machine-checked partial bounds, not solutions.

| Problem | Result | Verification |
|---|---|---|
| [#273](https://www.erdosproblems.com/273) | No covering system of ℤ with pairwise distinct moduli of the form p−1 ≤ 57 (p prime ≥ 5) | UNSAT certificates (LRAT, B=50 and B=57) each checked by 3 independent checkers incl. the formally verified [cake_lpr](https://github.com/tanyongkiam/cake_lpr) |
| [#385](https://www.erdosproblems.com/385) | F(n) > n for all 10⁴ < n ≤ 1.0011×10¹² except the two known values n = 267672, 267680 | Segmented sieve, deterministic checkpoints, shard cross-checks |
| [#647](https://www.erdosproblems.com/647) | No n > 24 with max_{m<n}(m + τ(m)) ≤ n + 2 up to 1.0011×10¹² | Same sieve pass; known anchor R(24)=26 reproduced |

## #273 — SAT certificate (`sat273/`)

A covering system with all moduli from the full admissible pool exists iff a
CNF is satisfiable: one variable per congruence class (n, a), pairwise
at-most-one per modulus (distinct moduli), one coverage clause per residue
mod L. Two certified bounds:

- **B=50**: pool {4, …, 46} (13 moduli, L = 1,275,120), CNF 310 vars /
  1,279,875 clauses — UNSAT.
- **B=57** (2026-07-22): pool gains 52 (=4·13), L = 16,576,560, CNF 362 vars /
  16,582,641 clauses — UNSAT.

- `export_dimacs.py` — deterministic CNF generator (~40 lines of logic);
  SHA-256 of each generated CNF is pinned in `erdos273_B{50,57}.meta.json`
- `artifacts/` — B=50 CNF + certificate, B=57 certificate (33 MB xz;
  regenerate the B=57 CNF with `--B 57`, 770 MB raw)
- `proof_B57_verification.json` — checker verdicts, hashes, timings
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

## Related Erdős repositories

This repo is certificate/computation only — SAT + LRAT and a segmented sieve,
no machine learning or theorem proving. Two other repos cover different
Erdős problems with different methods:

- [ProofForge](https://github.com/Sanexxxx777/ProofForge) — Lean 4 / Mathlib
  formal proofs (Erdős #1084, #1052), merged into Google DeepMind's
  `formal-conjectures`.
- [erdos-openevolve](https://github.com/Sanexxxx777/erdos-openevolve) — an
  evolutionary LLM coding pipeline (OpenEvolve/AlphaEvolve) that reproduces a
  numerical SOTA bound for the minimum overlap problem; not a formal proof.

## License

Apache-2.0

---

Built by Aleksandr Shulgin (@Aleksandr_NFA) · [shulgin.is-a.dev](https://shulgin.is-a.dev)
