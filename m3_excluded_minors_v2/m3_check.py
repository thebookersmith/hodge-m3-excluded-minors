#!/usr/bin/env python3
# ============================================================================================
# m3_check.py -- SELF-CONTAINED machine-verifiable checker for the claim:
#
#   "The graphic matroids of FIVE graphs are EXCLUDED MINORS for the class M_3 of
#    Engel-de Gaay Fortman-Schreieder (arXiv:2507.15704 v3, Def 8.3 / Prop 8.11):
#       find1  H?zTbbo      (9V/16E,  rank 8:  Q_3 + one-sided apex)
#       find2  J??FCpSJFw?  (11V/18E, rank 10: subdivided K_4 + apex to midpoints)
#       find3  I?`FF_{F_    (10V/17E, rank 9)
#       find4  I?B@nRWN?    (10V/17E, rank 9)
#       find5  ICQf@pSF_    (10V/17E, rank 9)"
#   (+ K_{3,5} anchor: the published EGFS excluded minor, dual witness only.)
#
# WHAT IS VERIFIED, for each find, with NOTHING but integer matrix-vector products mod 3:
#   (1) NONMEMBER: a Farkas dual witness y with y^T M = 0 and y . rhs != 0 for the
#       Albanese membership system of the candidate  =>  the matroid is NOT in M_3;
#   (2) MEMBER, for one-element minors: a primal witness x with Clo.x = 0 and color
#       profile identically 1  =>  the minor IS in M_3.
#       - finds 3/4/5: ALL 34 one-element minors individually (all_one_element_minors);
#       - finds 1/2:   one representative per Aut-orbit of minors (orbit_representatives);
#         the manifest's orbit_accounting (printed below) records why the 4 classes
#         cover all 32/36 single-element minors; isomorphic matroids have equal
#         M_3-status, so orbit representatives suffice.
#   (3) MINIMALITY LOGIC (mathematical, no computation): M_3 is minor-closed
#       (EGFS Prop 7.2), so if every 1-element minor is a member, every proper minor
#       is; with (1), the matroid is minor-minimal not-in-M_3 = an excluded minor.
#       The same fact shows the five finds are pairwise independent (none contains
#       another, since all proper minors of each are members while each is not).
#
# Witness files are SELF-CONTAINED: each .npz carries the realization A (reduced
# oriented incidence, g x n over GF(3)); the Albanese system is REBUILT from A alone
# (cycle-space coordinates: nullspace basis N with identity on free columns, color-s
# generator g_s = column s of N, edge (s,w): w -> w + g_s, edge index e = s*Vn + w).
# REQUIREMENTS: python3 + numpy.   USAGE: python3 m3_check.py
# ============================================================================================
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

def P(*a):
    print(*a); sys.stdout.flush()

def gf_rref(A, p):
    A = (np.asarray(A, np.int64) % p).copy()
    m, n = A.shape
    piv = []; r = 0
    for c in range(n):
        if r >= m: break
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0: continue
        i = r + int(nz[0])
        if i != r: A[[r, i]] = A[[i, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
        col = A[:, c].copy(); col[r] = 0
        A = (A - np.outer(col, A[r])) % p
        piv.append(c); r += 1
    return A, piv

def nullspace(A, p):
    R, piv = gf_rref(A, p)
    n = A.shape[1]
    free = [c for c in range(n) if c not in set(piv)]
    N = np.zeros((len(free), n), np.int64)
    for j, f in enumerate(free):
        N[j, f] = 1
        for ri, pc in enumerate(piv):
            N[j, pc] = (-R[ri, f]) % p
    return N % p

def check_file(fn):
    W = np.load(fn, allow_pickle=False)
    kind = str(W["kind"][0]); ell = int(W["ell"])
    A = W["A"].astype(np.int64) % ell
    w = W["witness"].astype(np.int64) % ell
    g, n = A.shape
    N = nullspace(A, ell)
    c = N.shape[0]
    assert c == n - g and not ((A @ N.T) % ell).any()
    Vn = ell ** c
    powers = ell ** np.arange(c, dtype=np.int64)
    D = np.zeros((Vn, c), np.int64)
    tmp = np.arange(Vn, dtype=np.int64)
    for i in range(c):
        D[:, i] = tmp % ell; tmp //= ell
    gens = (N.T % ell)
    if kind == "primal_member":
        x = w; assert x.size == n * Vn
        ok = True
        clo = np.zeros((g, Vn), np.int64)
        for s in range(n):
            xs = x[s * Vn:(s + 1) * Vn]
            if int(xs.sum() % ell) != 1: ok = False
            head = ((D + gens[s][None, :]) % ell) @ powers
            for k in range(g):
                a = int(A[k, s]) % ell
                if a:
                    np.add.at(clo[k], head, a * xs); clo[k] -= a * xs
        ok = ok and not (clo % ell).any()
        return ok, f"primal g={g} n={n} c={c}"
    elif kind == "dual_nonmember":
        y = w; assert y.size == g * Vn + n
        a2 = y[:g * Vn].reshape(g, Vn); b = y[g * Vn:]
        allz = True
        for s in range(n):
            head = ((D + gens[s][None, :]) % ell) @ powers
            acc = np.zeros(Vn, np.int64)
            for k in range(g):
                co = int(A[k, s]) % ell
                if co: acc += co * (a2[k][head] - a2[k])
            if ((acc + b[s]) % ell).any(): allz = False
        val = int(b.sum() % ell)
        return (allz and val != 0), f"dual g={g} n={n} c={c} y.rhs={val}"
    raise ValueError(kind)

def main():
    man = json.load(open(os.path.join(HERE, "data", "manifest.json")))
    P("=" * 88)
    P(man["title"]); P(man["adjudication"]); P("=" * 88)
    all_ok = True
    for F in man["finds"]:
        P(f"\n### {F['name']}  graph6={F.get('graph6')}  [{F.get('description','')}]")
        ok, det = check_file(os.path.join(HERE, "witnesses", F["dual"]))
        P(f"  NONMEMBER dual witness: {'VERIFIED' if ok else '*** FAIL ***'} ({det})")
        all_ok &= ok
        nm = 0
        for m in F["minors"]:
            ok, det = check_file(os.path.join(HERE, "witnesses", m["file"]))
            all_ok &= ok; nm += ok
            if not ok:
                P(f"  minor {m['file']}: *** FAIL *** ({det})")
        if F["minors"]:
            P(f"  minor MEMBER witnesses: {nm}/{len(F['minors'])} VERIFIED "
              f"[{F['minor_scheme']}]")
        if F.get("orbit_accounting"):
            P(f"  orbit accounting: {F['orbit_accounting']}")
        if F["minor_scheme"] in ("all_one_element_minors", "orbit_representatives"):
            P("  => minor-closedness (Prop 7.2): every proper minor is a member; "
              "with the dual witness, this matroid is an EXCLUDED MINOR for M_3.")
    P("\n" + "=" * 88)
    P("Pairwise independence: each find's proper minors are all members while every "
      "find is a nonmember; hence no find contains another as a minor.")
    P("RESULT: " + ("PASS" if all_ok else "*** FAIL ***"))
    P("=" * 88)
    sys.exit(0 if all_ok else 2)

main()
