#!/usr/bin/env python3
"""Verify the membership certificates at ell = 5 for the three M_3 excluded
minors for which the accompanying note reports ell = 5 computations:

    M(K_{3,5})  (8V/15E,  rank 7)
    M(G_1)      H?zTbbo      (9V/16E,  rank 8)
    M(G_2)      J??FCpSJFw?  (11V/18E, rank 10)

Each is shown to lie in the class M_5 of Engel-de Gaay Fortman-Schreieder
(arXiv:2507.15704v3, Def. 8.3). For each graph the script first identifies the
stored realization with the named graph (graph6 or edge list embedded here;
explicit isomorphism verified edge by edge), then reconstructs the full
5^8 = 390,625-vertex Albanese membership system from the realization and verifies
that the stored primal certificate solves the closedness equations with color
profile identically 1. Same design as m3_check.py, with ell = 5.

Requirements: Python 3 and numpy.   Usage: python3 l5_check.py   (a few seconds).
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

EMBED = {
    "K35": {"edges": [(i, 3 + j) for i in range(3) for j in range(5)], "nV": 8, "nE": 15,
            "file": "K35_ell5_member.npz", "label": "M(K_{3,5})"},
    "G1": {"graph6": "H?zTbbo", "nE": 16,
           "file": "G1_ell5_member.npz", "label": "M(G_1)"},
    "G2": {"graph6": "J??FCpSJFw?", "nE": 18,
           "file": "G2_ell5_member.npz", "label": "M(G_2)"},
}

def P(*a):
    print(*a); sys.stdout.flush()

# ---------------------------------------------------------------- graph layer (pure python)
def graph6_decode(s):
    b = [ord(c) - 63 for c in s]
    assert all(0 <= x < 64 for x in b), "graph6: byte out of range"
    n = b[0]
    assert 0 <= n < 63, "graph6: only short-form n supported here"
    bits = []
    for x in b[1:]:
        bits += [(x >> (5 - i)) & 1 for i in range(6)]
    need = n * (n - 1) // 2
    assert len(bits) >= need and not any(bits[need:]), "graph6: bad length/padding"
    E, k = [], 0
    for j in range(1, n):
        for i in range(j):
            if bits[k]:
                E.append((i, j))
            k += 1
    return n, sorted(E)

def _adj(n, E):
    a = [set() for _ in range(n)]
    for u, v in E:
        a[u].add(v); a[v].add(u)
    return a

def isomorphisms(n, E1, E2):
    if len(E1) != len(E2):
        return []
    a1, a2 = _adj(n, E1), _adj(n, E2)
    d1 = [len(x) for x in a1]; d2 = [len(x) for x in a2]
    if sorted(d1) != sorted(d2):
        return []
    order = sorted(range(n), key=lambda u: (-d1[u], u))
    out, phi, used = [], [-1] * n, [False] * n
    def bt(i):
        if i == n:
            if sorted(tuple(sorted((phi[u], phi[v]))) for u, v in E1) == E2:
                out.append(phi.copy())
            return bool(out)
        u = order[i]
        for v in range(n):
            if used[v] or d2[v] != d1[u]:
                continue
            ok = True
            for w in order[:i]:
                if (w in a1[u]) != (phi[w] in a2[v]):
                    ok = False; break
            if ok:
                phi[u] = v; used[v] = True
                if bt(i + 1):
                    return True
                phi[u] = -1; used[v] = False
        return False
    bt(0)
    return out

def support_graph(A, p):
    g, n = A.shape
    edges = []
    for j in range(n):
        rows = np.nonzero(A[:, j] % p)[0]
        if rows.size == 2:
            u, v = int(rows[0]), int(rows[1])
            if (int(A[u, j]) + int(A[v, j])) % p != 0:
                return None
            edges.append((u, v))
        elif rows.size == 1:
            edges.append((int(rows[0]), g))
        else:
            return None
    edges = [tuple(sorted(e)) for e in edges]
    if len(set(edges)) != len(edges):
        return None
    return g + 1, edges

def verify_binding(spec, A, p):
    sg = support_graph(A, p)
    if sg is None:
        return False, "A is not a reduced oriented incidence matrix"
    nv, E_sup = sg
    if "graph6" in spec:
        n_t, E_t = graph6_decode(spec["graph6"])
    else:
        n_t, E_t = spec["nV"], sorted(tuple(sorted(e)) for e in spec["edges"])
    if nv != n_t or len(E_sup) != spec["nE"] or len(E_t) != spec["nE"]:
        return False, "vertex/edge count mismatch"
    if not isomorphisms(nv, sorted(E_sup), E_t):
        return False, "no isomorphism support(A) ~ named graph"
    return True, f"isomorphism verified ({nv}V/{len(E_sup)}E)"

# ---------------------------------------------------------------- GF(p) linear algebra
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

# ---------------------------------------------------------------- primal verification
def check_member(fn):
    W = np.load(fn, allow_pickle=False)
    kind = str(W["kind"][0]); ell = int(W["ell"])
    assert kind == "primal_member" and ell == 5
    A = W["A"].astype(np.int64) % ell
    x = W["witness"].astype(np.int64) % ell
    g, n = A.shape
    N = nullspace(A, ell)
    c = N.shape[0]
    assert c == n - g and not ((A @ N.T) % ell).any()
    Vn = ell ** c
    assert x.size == n * Vn
    powers = ell ** np.arange(c, dtype=np.int64)
    D = np.zeros((Vn, c), np.int64)
    tmp = np.arange(Vn, dtype=np.int64)
    for i in range(c):
        D[:, i] = tmp % ell; tmp //= ell
    gens = (N.T % ell)
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
    return ok, A, f"g={g} n={n} c={c} Vn={Vn}"

def main():
    P("=" * 88)
    P("Membership certificates at l = 5: M(K_{3,5}), M(G_1), M(G_2) in M_5 "
      "(EGFS arXiv:2507.15704v3, Def. 8.3)")
    P("=" * 88)
    all_ok = True
    for key in ("K35", "G1", "G2"):
        spec = EMBED[key]
        P(f"\n{spec['label']}" + (f"  (graph6 {spec['graph6']})" if 'graph6' in spec else ""))
        ok, A, det = check_member(os.path.join(HERE, spec["file"]))
        bok, bdet = verify_binding(spec, A, 5)
        P(f"  graph identification: {'verified' if bok else 'FAILED'} ({bdet})")
        st = 'VERIFIED: MEMBER of M_5' if ok else 'FAILED'
        P(f"  full-system membership certificate ({det}): {st}")
        all_ok &= ok and bok
    P("\n" + "=" * 88)
    P("RESULT: " + ("PASS" if all_ok else "FAIL"))
    P("=" * 88)
    sys.exit(0 if all_ok else 2)

main()
