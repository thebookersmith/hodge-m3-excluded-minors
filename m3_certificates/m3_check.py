#!/usr/bin/env python3
"""Verify the certificates for five excluded minors of the class M_3 of
Engel-de Gaay Fortman-Schreieder (arXiv:2507.15704v3, Def. 8.3 / Prop. 8.11):

    G_1  H?zTbbo      (9V/16E,  rank 8)
    G_2  J??FCpSJFw?  (11V/18E, rank 10)
    G_3  I?`FF_{F_    (10V/17E, rank 9)
    G_4  I?B@nRWN?    (10V/17E, rank 9)
    G_5  ICQf@pSF_    (10V/17E, rank 9)

together with M(K_{3,5}) as a reference case. For each graph the script:

  1. identifies the stored realization with the stated graph -- the graph6
     strings are embedded below and decoded here, and an explicit vertex
     bijection to the support graph of the realization is verified edge by edge;
  2. verifies a dual certificate of nonmembership in M_3 (a left-kernel vector y
     with y^T M = 0 and y . rhs != 0 for the Albanese membership system);
  3. verifies primal membership certificates covering every single-element
     deletion and contraction -- for G_3, G_4, G_5 all 34 individually; for
     G_1, G_2 by automorphism-orbit representatives, with Aut(G) and its edge
     orbits computed here. Each stored minor certificate is matched to a
     realization reconstructed from the parent, by equal row space or by an
     explicit graph isomorphism.

Since M_3 is minor-closed (EGFS Prop. 7.2), membership of every single-element
minor gives membership of every proper minor; with the nonmembership certificate,
each graph is an excluded minor, and the five are pairwise incomparable. The
membership systems are reconstructed from the realizations alone; the
verification does not use the search or solver that found the certificates.

Requirements: Python 3 and numpy.   Usage: python3 m3_check.py
"""
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# Embedded graph identities -- the ground truth this checker certifies against.
# (K_{3,5} is given as an explicit edge list: parts {0,1,2} and {3..7}.)
EMBED = {
    "find1": {"graph6": "H?zTbbo", "nE": 16},
    "find2": {"graph6": "J??FCpSJFw?", "nE": 18},
    "find3": {"graph6": "I?`FF_{F_", "nE": 17},
    "find4": {"graph6": "I?B@nRWN?", "nE": 17},
    "find5": {"graph6": "ICQf@pSF_", "nE": 17},
    "K35_anchor": {"edges": [(i, 3 + j) for i in range(3) for j in range(5)],
                   "nV": 8, "nE": 15},
}

DISPLAY = {"find1": "G_1", "find2": "G_2", "find3": "G_3", "find4": "G_4",
           "find5": "G_5", "K35_anchor": "M(K_{3,5})"}

def P(*a):
    print(*a); sys.stdout.flush()

# ---------------------------------------------------------------- graph layer (pure python)
def graph6_decode(s):
    """Decode a graph6 string to (n, sorted edge list). From-scratch, format per McKay."""
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

def isomorphisms(n, E1, E2, find_all=False):
    """Explicit backtracking graph isomorphism: returns list of vertex bijections phi
    (phi[u] in G2 for u in G1) with phi(E1) = E2. Small n only; exhaustive and exact."""
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
            # full independent re-verification of the candidate bijection
            if sorted(tuple(sorted((phi[u], phi[v]))) for u, v in E1) == E2:
                out.append(phi.copy())
            return bool(out) and not find_all
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
    """Interpret A (g x n over GF(p)) as a reduced oriented incidence matrix: each column
    must have exactly two nonzeros a, -a (edge between those rows) or one nonzero (edge
    from that row to the deleted root vertex, index g). Returns (nV, edge list, col->edge)
    or None if A is not of this shape or the graph is not simple."""
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

def verify_binding(name, A, p):
    """A is a realization of the graphic matroid of the graph named `name` (embedded)."""
    spec = EMBED[name]
    sg = support_graph(A, p)
    if sg is None:
        return False, "A is not a reduced oriented incidence matrix"
    nv, E_sup = sg
    if "graph6" in spec:
        n_t, E_t = graph6_decode(spec["graph6"])
    else:
        n_t, E_t = spec["nV"], sorted(tuple(sorted(e)) for e in spec["edges"])
    if nv != n_t or len(E_sup) != spec["nE"] or len(E_t) != spec["nE"]:
        return False, f"vertex/edge count mismatch ({nv}V/{len(E_sup)}E vs {n_t}V/{len(E_t)}E)"
    phis = isomorphisms(nv, sorted(E_sup), E_t)
    if not phis:
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

def gf_rank(A, p):
    return len(gf_rref(A, p)[1])

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

def same_rowspace(B1, B2, p):
    if B1.shape[1] != B2.shape[1]:
        return False
    r1, r2 = gf_rank(B1, p), gf_rank(B2, p)
    return r1 == r2 == gf_rank(np.vstack([B1 % p, B2 % p]), p)

def minor_graph(nv, edges, e, op):
    """Combinatorial one-element minor of the support graph (labels compacted).
    Returns (nV, sorted edge list), or None if the minor is not simple."""
    rest = [edges[j] for j in range(len(edges)) if j != e]
    if op == "delete":
        out, n2 = rest, nv
    else:
        u, v = edges[e]
        def m(w): return u if w == v else (w - 1 if w > v else w)
        out = [tuple(sorted((m(a), m(b)))) for a, b in rest]
        n2 = nv - 1
    if len(set(out)) != len(out) or any(a == b for a, b in out):
        return None
    return n2, sorted(out)

def minor_realization(A, e, op, p):
    """Derive the realization of the one-element minor (delete/contract edge = column e)
    from the parent realization A, inside this checker."""
    A = np.asarray(A, np.int64) % p
    if op == "delete":
        return np.delete(A, e, axis=1)
    col = A[:, e].copy()
    u = int(np.nonzero(col)[0][0])
    factor = (col * pow(int(col[u]), p - 2, p)) % p
    B = (A - np.outer(factor, A[u, :])) % p
    return np.delete(np.delete(B, u, axis=0), e, axis=1)

# ---------------------------------------------------------------- witness verification
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
        return ok, A, f"primal g={g} n={n} c={c}"
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
        return (allz and val != 0), A, f"dual g={g} n={n} c={c} y.rhs={val}"
    raise ValueError(kind)

# ---------------------------------------------------------------- per-find driver
def edge_orbits(nv, E, ncols, col_edges):
    """Compute the orbits of Aut(G) on edges (as column indices), from scratch."""
    autos = isomorphisms(nv, E, E, find_all=True)
    idx = {e: j for j, e in enumerate(col_edges)}
    seen, orbits = set(), []
    for j in range(ncols):
        if j in seen: continue
        orb = set()
        for phi in autos:
            u, v = col_edges[j]
            orb.add(idx[tuple(sorted((phi[u], phi[v])))])
        orbits.append(sorted(orb)); seen |= orb
    return len(autos), orbits

def run_find(F, ell=3):
    name = F["name"]
    label = DISPLAY.get(name, name)
    g6 = F.get("graph6")
    P(f"\n{label}" + (f"  (graph6 {g6})" if g6 else ""))
    all_ok = True

    # manifest graph6 must agree with the embedded ground truth
    if "graph6" in EMBED[name] and F.get("graph6") != EMBED[name]["graph6"]:
        P("  FAILED: manifest graph6 disagrees with the graph6 embedded in this checker")
        return False

    ok, A_parent, det = check_file(os.path.join(HERE, "witnesses", F["dual"]))
    bok, bdet = verify_binding(name, A_parent, ell)
    P(f"  graph identification: {'verified' if bok else 'FAILED'} ({bdet})")
    P(f"  nonmembership certificate: {'verified' if ok else 'FAILED'} ({det})")
    all_ok &= ok and bok
    if not bok:
        P("  (skipping minor coverage: identification failed)")
        return False
    if not F["minors"]:
        return all_ok

    # verify every stored minor witness once; keep its realization + support graph
    stored = []
    for m in F["minors"]:
        ok, A_m, det = check_file(os.path.join(HERE, "witnesses", m["file"]))
        all_ok &= ok
        if not ok:
            P(f"  minor certificate {m['file']}: FAILED ({det})")
        stored.append((m["file"], A_m, ok, support_graph(A_m, ell)))

    g, n = A_parent.shape
    sg = support_graph(A_parent, ell)
    nv, E_sup = sg
    col_edges = [tuple(sorted(e)) for e in E_sup]

    def covered(op, e):
        # tier 1: identical GF(3) row space against the internally derived realization
        target = minor_realization(A_parent, e, op, ell)
        for (_, A_m, ok, _sg) in stored:
            if ok and same_rowspace(target, A_m, ell):
                return True
        # tier 2: the stored file is a reduced incidence of a graph explicitly
        # isomorphic to the derived minor graph (labeling-independent, still exact)
        tg = minor_graph(nv, col_edges, e, op)
        if tg is None:
            return False
        tnv, tE = tg
        for (_, A_m, ok, sgm) in stored:
            if not ok or sgm is None:
                continue
            mnv, mE = sgm
            if mnv == tnv and len(mE) == len(tE) and isomorphisms(tnv, sorted(mE), tE):
                return True
        return False

    if F["minor_scheme"] == "all_one_element_minors":
        miss = [(op, e) for op in ("delete", "contract") for e in range(n)
                if not covered(op, e)]
        if miss:
            P(f"  FAILED: single-element minors not covered by a verified certificate: {miss}")
            all_ok = False
        else:
            P(f"  minor coverage: all {2*n} single-element minors verified.")
    elif F["minor_scheme"] == "orbit_representatives":
        naut, orbits = edge_orbits(nv, sorted(E_sup), n, col_edges)
        bad = []
        for orb in orbits:
            for op in ("delete", "contract"):
                if not any(covered(op, e) for e in orb):
                    bad.append((op, f"orbit{orbits.index(orb)}"))
        if bad:
            P(f"  FAILED: orbit classes without a verified representative: {bad}")
            all_ok = False
        else:
            P(f"  minor coverage: all {2*n} single-element minors verified, via the "
              f"{naut} automorphisms of the graph (computed here).")
    if all_ok:
        P(f"  conclusion: {label} is an excluded minor of M_3.")
    return all_ok

def main():
    man = json.load(open(os.path.join(HERE, "data", "manifest.json")))
    P("=" * 88)
    P(man["title"])
    P("=" * 88)
    all_ok = True
    for F in man["finds"]:
        all_ok &= run_find(F)
    P("\n" + "=" * 88)
    P("Each of G_1..G_5 is a nonmember of M_3 whose single-element minors are all "
      "members; hence none contains another as a minor.")
    P("RESULT: " + ("PASS" if all_ok else "FAIL"))
    P("=" * 88)
    sys.exit(0 if all_ok else 2)

main()
