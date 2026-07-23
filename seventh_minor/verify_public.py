#!/usr/bin/env python3
# verify_public.py -- run the PUBLIC checker's primal verification (m3_check.check_file)
# and graph identification (support_graph + isomorphisms) on the seventh-minor witnesses.
# m3_check.py auto-runs its manifest on import, so we exec a truncated copy (everything
# before `def main(`) and call the functions directly -- identical source, no side effects.
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
M3 = os.path.abspath(os.path.join(HERE, "..", "m3_certificates", "m3_check.py"))

_src = open(M3).read()
_src = _src[:_src.index("def main(")]
CK = {"__file__": M3}
exec(_src, CK)
check_file = CK["check_file"]
support_graph = CK["support_graph"]
isomorphisms = CK["isomorphisms"]
graph6_decode = CK["graph6_decode"]

# K_{1,2,2,2} and its two deletion-class representatives (verified: 2 iso-classes)
K1222 = [(0,2),(0,3),(0,4),(0,5),(0,6),(1,2),(1,3),(1,4),(1,5),(1,6),
         (2,4),(2,5),(2,6),(3,4),(3,5),(3,6),(4,6),(5,6)]
def _del(E, e): return [x for j, x in enumerate(E) if j != e]

TARGETS = {
    "K6_member":       {"nV": 6, "edges": [(i, j) for i in range(6) for j in range(i+1, 6)],
                        "desc": "K6 (complete graph on 6 vertices)"},
    "K1222_del_oct":   {"nV": 7, "edges": _del(K1222, 10),
                        "desc": "K_{1,2,2,2} minus octahedron edge (2,4); degseq (6,5,5,5,5,4,4)"},
    "K1222_del_apx":   {"nV": 7, "edges": _del(K1222, 4),
                        "desc": "K_{1,2,2,2} minus apex edge (0,6); degseq (5,5,5,5,5,5,4)"},
}

def P(f, *a):
    s = " ".join(str(x) for x in a)
    print(s); f.write(s + "\n"); f.flush()

def run(f, name, ell=3):
    spec = TARGETS[name]
    fn = os.path.join(HERE, "witnesses", name + "_primal.npz")
    P(f, f"\n{name}  ({spec['desc']})")
    P(f, f"  file: witnesses/{name}_primal.npz")
    ok, A, det = check_file(fn)
    P(f, f"  primal membership certificate: {'VERIFIED' if ok else 'FAILED'}  ({det})")
    # graph identification: support(A) ~ named simple graph
    sg = support_graph(A % ell, ell)
    bok = False
    if sg is None:
        P(f, "  graph identification: FAILED (A is not a reduced oriented incidence matrix)")
    else:
        nv, E_sup = sg
        n_t, E_t = spec["nV"], sorted(tuple(sorted(e)) for e in spec["edges"])
        if nv != n_t or len(E_sup) != len(E_t):
            P(f, f"  graph identification: FAILED (size {nv}V/{len(E_sup)}E vs {n_t}V/{len(E_t)}E)")
        elif isomorphisms(nv, sorted(E_sup), E_t):
            P(f, f"  graph identification: verified (support(A) ~ target, {nv}V/{len(E_sup)}E)")
            bok = True
        else:
            P(f, "  graph identification: FAILED (no isomorphism)")
    return ok and bok

def main():
    out = os.path.join(HERE, "K1222_CHECK_OUTPUT.txt")
    with open(out, "a", encoding="utf-8") as f:
        P(f, "=" * 88)
        P(f, "PUBLIC-CHECKER primal verification (m3_check.check_file) -- seventh excluded minor")
        P(f, "K_{1,2,2,2} closure: K6 membership + two deletion-class witnesses")
        import datetime
        P(f, "run:", datetime.datetime.now().isoformat(timespec="seconds"))
        P(f, "=" * 88)
        allok = True
        names = sys.argv[1:] or ["K6_member", "K1222_del_oct", "K1222_del_apx"]
        for nm in names:
            try:
                allok &= run(f, nm)
            except FileNotFoundError:
                P(f, f"\n{nm}: witness file not yet present (skipped)")
                allok = False
        P(f, "\n" + "=" * 88)
        P(f, "RESULT: " + ("PASS -- all present witnesses verified by the public checker"
                           if allok else "INCOMPLETE / FAIL"))
        P(f, "=" * 88)

if __name__ == "__main__":
    main()
