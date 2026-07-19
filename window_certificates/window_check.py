#!/usr/bin/env python3
# ============================================================================================
# window_check.py -- checker for socle-window certificates of nonmembership in the
# EGFS classes M_ell (arXiv:2507.15704v3, Def. 8.3), for two systems too large to
# construct directly:
#
#     K_{3,9} at ell=5:  corank 16, full system has 5^16 = 152,587,890,625 vertices
#     K_8     at ell=3:  corank 21, full system has 3^21 =  10,460,353,203 vertices
#
# Requirements: Python 3 and numpy.  Usage: python3 window_check.py
# The checker needs no other file from the repository: it rebuilds all data from the
# edge lists below and reconstructs each realization from first principles.
#
# --------------------------------------------------------------------------------------------
# What a certificate is, and why the window identity proves nonmembership
# (the correspondence between window-ring identities and full-system dual certificates)
#
# The Albanese membership system for a graph G with reduced oriented incidence A (g x n
# over F_ell) lives on the vertex group V = F_ell^c (c = corank): writing N for a basis
# of the cycle space (A N^T = 0, with identity on the free columns), each edge s has a
# generator g_s = column s of N^T, and the (5.3)-graph has an edge w -> w + g_s of color
# s for every w in V.  G is in M_ell iff a "profile-1" 1-chain exists (Def 8.3).  By the
# Farkas alternative, G is NOT in M_ell iff there exist functions a_1..a_g : V -> F_ell
# and b in F_ell^n with, for EVERY edge s and EVERY w in V:
#     sum_k A[k,s] * ( a_k(w + g_s) - a_k(w) )  +  b_s  =  0,        (*)
# and sum_s b_s != 0.
#
# DICTIONARY.  Encode a function a : V -> F_ell as the group-algebra element
# a^ := sum_w a(w) x^w  in  R = F_ell[x_1..x_c]/(x_i^ell - 1) = F_ell[t_1..t_c]/(t_i^ell),
# where t_i = x_i - 1 (the second equality is the mod-ell binomial identity
# (x-1)^ell = x^ell - 1).  Then:
#   . the shifted function w |-> a(w + g)  corresponds to  x^{-g} * a^ ;
#   . the CONSTANT function 1 corresponds to sum_w x^w = prod_i (1 + x_i + .. + x_i^{ell-1})
#     = prod_i t_i^{ell-1} =: t^socle  (the socle monomial), because
#     (x^ell - 1)/(x - 1) = t^{ell-1} in F_ell[t]/(t^ell);
#   . therefore the per-edge identity (*) -- "this function of w is the constant -b_s" --
#     is precisely the ring identity
#     sum_k A[k,s] * (x^{-g_s} - 1) * a_k^  +  b_s * t^socle  =  0   in R.        (**)
#
# WINDOW.  If every a_k^ is supported on t-monomials of degree >= (ell-1)c - D (the
# "socle window", depth D), then so is every term of (**): multiplication by
# (x^{-g} - 1) = prod_i (1+t_i)^{(-g)_i mod ell} - 1 never LOWERS any t-exponent.  Hence
# (**) can be checked entirely inside the window -- a finite computation of size
# (number of window monomials), independent of ell^c.  A verified (**) for all s, plus
# sum_s b_s != 0, yields a full-system dual certificate via the correspondence read in
# reverse, so the graph is not in M_ell.  No step of this argument requires evaluating
# any function on the ell^c vertices.
#
# COORDINATES USED IN THE FILES.  Window monomials are stored as DEFICIENCY tuples
# delta = (ell-1, .., ell-1) - alpha (so the socle is delta = 0), listed in the array
# `mons`; a_k^ is stored as its coefficient row awin[k] over `mons`.  Multiplication by
# (x^gamma - 1) in deficiency coordinates sends the monomial delta to
#     sum_{0 < d <= delta (componentwise)}  prod_i C(gamma_i, d_i)  *  (delta - d),
# because (1+t_i)^{gamma_i} t_i^{alpha_i} = sum_{d_i} C(gamma_i, d_i) t_i^{alpha_i+d_i}
# and exponents above ell-1 vanish (t_i^ell = 0) -- i.e. negative deficiencies are dropped.
#
# For each claim the checker (1) rebuilds the incidence matrix from the edge list
# embedded below, (2) row-reduces it to the realization A by plain Gauss-Jordan
# (leftmost-pivot RREF), (3) requires it to equal the realization stored in the witness
# file, (4) recomputes the cycle-space generators, and (5) verifies (**) monomial by
# monomial together with sum_s b_s != 0.  The only data taken from the witness file are
# the certificate coefficients being checked.
# ============================================================================================
import os, sys
import numpy as np
from math import comb
from itertools import product as iproduct

HERE = os.path.dirname(os.path.abspath(__file__))

CLAIMS = [
    {"name": "K_{3,9} is not in M_5   (the p=5 nonmembership underlying Engel-Schreieder "
             "Thm 1.7; certificate depth 5)",
     "file": "witnesses/K39_ell5_windowD5_dual.npz",
     "ell": 5,
     "edges": [(i, 3 + j) for i in range(3) for j in range(9)],
     "nverts": 12},
    {"name": "K_8 is not in M_3   (depth-3 window certificate at corank 21)",
     "file": "witnesses/K8_ell3_windowD3_dual.npz",
     "ell": 3,
     "edges": [(i, j) for i in range(8) for j in range(i + 1, 8)],
     "nverts": 8},
]

def P(*a):
    print(*a); sys.stdout.flush()

# ------------------------------------------------------------------ exact GF(p) RREF
def gf_rref(M, p):
    M = (np.asarray(M, np.int64) % p).copy()
    m, n = M.shape
    piv = []
    r = 0
    for c in range(n):
        if r >= m:
            break
        nz = np.nonzero(M[r:, c])[0]
        if nz.size == 0:
            continue
        i = r + int(nz[0])
        if i != r:
            M[[r, i]] = M[[i, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        col = M[:, c].copy(); col[r] = 0
        M = (M - np.outer(col, M[r])) % p
        piv.append(c); r += 1
    return M, piv

def realization_from_edges(nverts, edges, p):
    """Oriented incidence (+1 at u, -1 at v per edge (u,v)), then keep the row set
    selected by RREF of the transpose (leftmost pivots) -- the producing convention."""
    n = len(edges)
    A0 = np.zeros((nverts, n), np.int64)
    for j, (u, v) in enumerate(edges):
        A0[u, j] += 1
        A0[v, j] -= 1
    _, prow = gf_rref(A0.T % p, p)          # pivot columns of A0^T = independent rows of A0
    return A0[sorted(prow), :] % p

def cycle_generators(A, p):
    """Nullspace basis with identity on free columns; generator g_s = column s of N^T."""
    R, piv = gf_rref(A, p)
    n = A.shape[1]
    free = [c for c in range(n) if c not in set(piv)]
    N = np.zeros((len(free), n), np.int64)
    for j, f in enumerate(free):
        N[j, f] = 1
        for ri, pc in enumerate(piv):
            N[j, pc] = (-R[ri, f]) % p
    return (N % p).T % p                     # n x c : row s = g_s

# ------------------------------------------------------------------ window arithmetic
def mult_by_xgamma_minus_1(coefs, gamma, p, c):
    """(x^gamma - 1) * (window element) in deficiency coordinates.
    coefs: dict delta-tuple -> coefficient.  See the header derivation."""
    out = {}
    for delta, v in coefs.items():
        if v % p == 0:
            continue
        for d in iproduct(*[range(x + 1) for x in delta]):
            if not any(d):
                continue                      # the d = 0 term cancels against the "-1"
            co = v
            for i in range(c):
                co = (co * comb(int(gamma[i]), d[i])) % p
                if co == 0:
                    break
            if co:
                tgt = tuple(delta[i] - d[i] for i in range(c))
                out[tgt] = (out.get(tgt, 0) + co) % p
    return out

def check_claim(cl):
    P("-" * 88)
    P("CLAIM:", cl["name"])
    ell = cl["ell"]
    W = np.load(os.path.join(HERE, cl["file"]), allow_pickle=False)
    assert str(W["kind"][0]) == "dual_nonmember_window" and int(W["ell"]) == ell
    A_stored = W["A"].astype(np.int64) % ell
    awin = W["awin"].astype(np.int64) % ell
    b = W["b"].astype(np.int64) % ell
    mons = [tuple(int(x) for x in row) for row in W["mons"]]
    # provenance: rebuild the realization from the edge list; must equal the stored one
    A = realization_from_edges(cl["nverts"], cl["edges"], ell)
    same = A.shape == A_stored.shape and not ((A - A_stored) % ell).any()
    P(f"  realization rebuilt from edge list ({cl['nverts']}V/{len(cl['edges'])}E): "
      f"g={A.shape[0]} n={A.shape[1]}; matches witness file: {'YES' if same else 'NO'}")
    if not same:
        return False
    g, n = A.shape
    c = n - g
    gens = cycle_generators(A, ell)
    P(f"  corank c={c}: full system would have {ell}^{c} = {ell**c:,} vertices "
      f"(not constructed); window monomials: {len(mons)}")
    socle = tuple([0] * c)
    ok_all = True
    for s in range(n):
        gamma = tuple(int((-int(x)) % ell) for x in gens[s])
        acc = {}
        for k in range(g):
            co = int(A[k, s]) % ell
            if co == 0:
                continue
            coefs = {mons[j]: int(awin[k, j]) for j in range(len(mons))
                     if awin[k, j] % ell}
            part = mult_by_xgamma_minus_1(coefs, gamma, ell, c)
            for key, v in part.items():
                acc[key] = (acc.get(key, 0) + co * v) % ell
        acc[socle] = (acc.get(socle, 0) + int(b[s])) % ell
        residue = {k: v for k, v in acc.items() if v % ell}
        if residue:
            P(f"  edge {s}: identity (**) FAILS, residue at {list(residue.items())[:3]}")
            ok_all = False
    sb = int(b.sum() % ell)
    P(f"  ring identity (**) for all {n} edges: {'PASS' if ok_all else 'FAIL'}")
    P(f"  sum_s b_s = {sb} (must be nonzero): {'PASS' if sb != 0 else 'FAIL'}")
    verdict = ok_all and sb != 0
    P(f"  => {'verified: not in M_%d (window dual certificate)' % ell if verdict else 'CERTIFICATE INVALID'}")
    return verdict

def main():
    ok = True
    for cl in CLAIMS:
        ok &= check_claim(cl)
    P("=" * 88)
    P("RESULT: " + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 2)

main()
