"""
Совместное сегментированное решето для задач Erdos #647 и #385.
Один проход по блокам m=1..TARGET, два независимых детектора на общих
tau(m)/spf(m) массивах.

=== Формализация #647 (Erdos-Selfridge) === (см. 647.txt)
Точный текст: пусть tau(n) -- число делителей n. Существует ли n>24 такое что
    max_{m<n} (m + tau(m)) <= n+2 ?
Известно, что для n=24 это верно. Граница n+2 точная: max(tau(n-1)+n-1,
tau(n-2)+n-2) >= n+2 для всех n (доказано в тексте, не нужно проверять).
Эрдёш: "крайне сомнительно", что таких n>24 бесконечно много; предложил $44
за ОДИН такой пример (задача открыта).

Детектор: R(n) = max_{m<n} (m+tau(m)) -- бегущий (неубывающий) рекорд.
Ищем n>24 c R(n) <= n+2. Находка = решение открытой проблемы.

=== Формализация #385 (Erdos, Eggleton, Selfridge) === (см. 385.txt)
Точный текст: F(n) = max_{m<n, m составное} (m + p(m)), где p(m) -- наименьший
простой делитель m. Вопросы: (1) верно ли F(n) > n для всех достаточно
больших n? (2) стремится ли F(n)-n к бесконечности?
Тривиально F(n) <= n + sqrt(n) (т.к. p(m) <= sqrt(m) для составного m).

Детектор: F(n) -- бегущий рекорд ТОЛЬКО по составным m (простые и m=1 не
участвуют в максимуме). Ищем n c F(n) <= n -- нарушение гипотезы (1).
Малые n игнорируются как находки (см. SMALL_N_385) -- текст явно говорит
"for all sufficiently large n", т.е. отдельные нарушения на малых n ожидаемы.
Параллельно трекаем running min (F(n)-n) по большим n -- индикатор для (2).

=== РАСХОЖДЕНИЯ С ТЗ ===
- ТЗ формулирует цель 647 как "n, опровергающее утверждение". Неточно: сам
  вопрос задачи -- "существует ли n>24..."; найденный n ОТВЕЧАЕТ на вопрос
  утвердительно (решение открытой проблемы, приз Эрдёша $44), это не
  "опровержение". Взято по тексту, как требует инструкция.
- ТЗ для 385 не уточняет "большие n"; текст явно ограничивает вопрос
  "sufficiently large n" -- реализовано через порог SMALL_N_385.

=== Шардинг (--lo/--hi) ===
Бегущие рекорды R(n)/F(n) -- глобальные (зависят от ВСЕЙ истории m=1..n).
Независимый шард не знает истинный рекорд на своей левой границе, поэтому:
шард реально решетит с (lo-overlap), чтобы "разогреть" локальный рекорд, но
находки логирует только начиная с m=lo. Находки с m-lo < overlap помечаются
near_shard_boundary=true -- локальный рекорд там мог не успеть "догнать"
истинный (нужна перепроверка склейкой с соседним шардом).
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
    """tau(m) для m в [lo,hi) через парную схему делителей:
    tau(m) = 2*|{d<=sqrt(m): d|m}| - [m -- полный квадрат].
    counts -- uint16 (tau(m)<~11000 для m<1e13, помещается с запасом),
    финальный tau тоже uint16 (экономия памяти на больших блоках)."""
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
    """Наименьший простой делитель m для m в [lo,hi). 0 если m простое или m=1.
    int32 -- значения до sqrt(1e13)~3.16e6, uint16 не хватит."""
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
    """Один блок: считает tau/spf, обновляет бегущие рекорды (in-place, чтобы
    не плодить лишние int64-копии), логирует находки с m>=report_from.
    Возвращает (max647, max385, min_gap385, n_findings647, n_findings385)."""
    tau = block_tau(lo, hi)
    D_block = math.isqrt(hi - 1)
    idx = np.searchsorted(primes, D_block, side="right")
    spf = block_spf(lo, hi, primes[:idx])

    m = np.arange(lo, hi, dtype=np.int64)
    n_vals = m + 1
    report_mask = m >= report_from

    # --- #647 --- (running647 переиспользует буфер excess647: accumulate/maximum с out=)
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
                     help="Последовательный режим: решетить m=1..target (checkpoint.json/findings.jsonl)")
    ap.add_argument("--block", type=int, default=BLOCK, help="Размер блока, чисел")
    ap.add_argument("--reset", action="store_true")
    ap.add_argument("--lo", type=int, default=None, help="Шард: m>=lo репортится (m-диапазон, не n)")
    ap.add_argument("--hi", type=int, default=None, help="Шард: конец m-диапазона (искл.)")
    ap.add_argument("--shard-overlap", type=int, default=1_000_000,
                     help="Насколько раньше lo реально начать решетить для разогрева рекорда")
    args = ap.parse_args()

    if args.lo is not None or args.hi is not None:
        if args.lo is None or args.hi is None:
            sys.exit("--lo и --hi нужны вместе")
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
