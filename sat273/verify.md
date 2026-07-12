# Verifying the #273 certificate end-to-end

Nothing here requires trusting us or the solver: regenerate the CNF from the
40-line deterministic generator, compare hashes, and check the LRAT certificate
with independently developed checkers (one of them formally verified).

## 1. Regenerate and pin the CNF

```bash
python3 export_dimacs.py --B 50        # writes erdos273_B50.cnf + .meta.json
shasum -a 256 erdos273_B50.cnf         # must equal cnf_sha256 in meta.json:
# e468d74bf5860a4ee0a58f21ab2f0ac5d6ad7825dc9573f2f4892e73bf6b8606
```

The generator enumerates the *full* pool of admissible moduli (it aborts if an
lcm cap would truncate the pool — a certificate for a partial pool proves
nothing) and writes one at-most-one block per modulus plus one coverage clause
per residue mod L.

## 2. (Optional) re-solve

```bash
brew install cadical                   # 3.0.0
cadical --plain --lrat --no-binary erdos273_B50.cnf proof.lrat
# prints: s UNSATISFIABLE  (~1 min single-thread)
```

Or use `artifacts/proof_B50_plain.lrat.xz` from this repo (`xz -d` it).

## 3. Check the certificate — three independent checkers

```bash
# (a) lrat-check (drat-trim, Heule)
git clone https://github.com/marijnheule/drat-trim && make -C drat-trim
drat-trim/lrat-check erdos273_B50.cnf proof_B50_plain.lrat   # → c VERIFIED

# (b) lrat-trim (Biere)
git clone https://github.com/arminbiere/lrat-trim
cd lrat-trim && ./configure && make && cd ..
lrat-trim/lrat-trim erdos273_B50.cnf proof_B50_plain.lrat    # → s VERIFIED

# (c) cake_lpr — formally verified checker (HOL4/CakeML)
git clone https://github.com/tanyongkiam/cake_lpr
cd cake_lpr && cc basis_ffi.c cake_lpr_arm8.S -o cake_lpr && cd ..
#   (x86-64: use cake_lpr.S instead of cake_lpr_arm8.S)
cake_lpr/cake_lpr erdos273_B50.cnf proof_B50_plain.lrat      # → s VERIFIED UNSAT
```

## Why `--plain`

CaDiCaL 3.0.0 with default inprocessing emits LRAT lines that introduce
extension variables (index > n_vars, empty hints). `lrat-check` silently
accepts such a file, `lrat-trim` correctly rejects it. With `--plain` the
certificate is plain RUP-only LRAT and all three checkers accept it. This is
why the published certificate was produced with `--plain`, and why we recommend
checking any LRAT with more than one checker.

## What exactly is proved

UNSAT of this CNF ⟺ no way to pick at most one residue class per modulus from
{4, 6, 10, 12, 16, 18, 22, 28, 30, 36, 40, 42, 46} covering all residues mod
L = 1,275,120. Since every modulus divides L, covering all residues mod L is
equivalent to covering ℤ. Distinctness of moduli is enforced by the at-most-one
blocks. Hence: **no covering system of ℤ with pairwise distinct moduli, each of
the form p−1 ≤ 50 for a prime p ≥ 5.** Moduli larger than 50 are not
constrained — the general problem #273 remains open.
