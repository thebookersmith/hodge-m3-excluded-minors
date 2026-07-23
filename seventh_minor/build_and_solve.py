# build_and_solve.py -- full-system GF(3) primal MEMBER witnesses for the seventh
# excluded-minor closure (K_{1,2,2,2}).  Targets:
#   K6      = complete graph on 6 vertices          (member; covers all contractions)
#   del_oct = K_{1,2,2,2} minus one octahedron edge (7V/17E, deletion class A)
#   del_apx = K_{1,2,2,2} minus one apex edge        (7V/17E, deletion class B)
# Route: mfx_core system -> full sparse M x = b over GF(3) -> SpaSM `solve` (feed M^T,
# b as 1-row SMS) -> exact scipy matvec verification -> save npz (kind='primal_member',
# ell, A=int8, witness=int8) in the public find*-file schema.
import sys, os, time, subprocess
sys.stdout.reconfigure(encoding="utf-8")
CODE = "/mnt/c/Users/Booker/Desktop/hodge/matroid/code"
sys.path.insert(0, CODE)
import numpy as np
import scipy.sparse as sp
from mfx_core import (graph_incidence, full_rank_rows, gf_nullspace, build_full,
                      rhs_full, write_sms, check_primal, system_data)

SOLVE = "/root/spasm_build/build/tools/solve"
SMSDIR = "/root/mfx_lab/sms"
OUT = "/mnt/c/Users/Booker/Desktop/hodge/m3_repo_staging/seventh_minor_staging/witnesses"
ELL = 3

# ---- graphs -----------------------------------------------------------------
def K6():
    return 6, [(i, j) for i in range(6) for j in range(i + 1, 6)]

# K_{1,2,2,2}: parts {0,1},{2,3},{4,5},{6}; graph6 F]~vw
K1222 = [(0,2),(0,3),(0,4),(0,5),(0,6),(1,2),(1,3),(1,4),(1,5),(1,6),
         (2,4),(2,5),(2,6),(3,4),(3,5),(3,6),(4,6),(5,6)]

def del_edge(edges, e):
    return 7, [x for j, x in enumerate(edges) if j != e]

def degseq(nv, E):
    d = [0]*nv
    for u, v in E:
        d[u]+=1; d[v]+=1
    return tuple(sorted(d, reverse=True))

def p(*a):
    print(*a); sys.stdout.flush()

def read_sms_row(fn, n, ell):
    x = np.zeros(n, np.int64)
    with open(fn) as f:
        f.readline()
        for line in f:
            parts = line.split()
            if len(parts) != 3 or parts[0] == "0":
                continue
            i, j, v = int(parts[0]), int(parts[1]), int(parts[2])
            assert i == 1
            x[j - 1] = v % ell
    return x

def solve_and_save(name, nv, edges, ell, threads=8):
    S = system_data(nv, edges, ell)
    p(f"\n[{name}] nv={nv} n={S['n']} g={S['g']} c={S['c']} Vn={S['Vn']} "
      f"degseq={degseq(nv, edges)}")
    t0 = time.time()
    Clo, Col = build_full(S)
    M = sp.vstack([Clo, Col]).tocsr(); M.data %= ell; M.eliminate_zeros(); M.sort_indices()
    b = rhs_full(S).astype(np.int64)
    p(f"[{name}] system M {M.shape} nnz={M.nnz} built ({time.time()-t0:.1f}s)")
    os.makedirs(SMSDIR, exist_ok=True)
    base = f"{SMSDIR}/_solve_{name}"
    t0 = time.time()
    write_sms(base + "_MT.sms", M.T.tocsr(), ell)
    bco = sp.csr_matrix((b % ell).reshape(1, -1))
    write_sms(base + "_b.sms", bco, ell)
    p(f"[{name}] SMS written ({time.time()-t0:.1f}s); launching SpaSM solve...")
    t0 = time.time()
    r = subprocess.run([SOLVE, "-p", str(ell), "-m", base + "_MT.sms",
                        "--rhs", base + "_b.sms", "-o", base + "_x.sms"],
                       capture_output=True, text=True, timeout=36000,
                       env=dict(os.environ, OMP_NUM_THREADS=str(threads)))
    dt = time.time() - t0
    p(f"[{name}] SpaSM solve rc={r.returncode} ({dt:.1f}s)")
    if r.stderr.strip():
        p(f"[{name}] solve stderr tail: {r.stderr.strip()[-300:]}")
    if not os.path.exists(base + "_x.sms"):
        p(f"[{name}] SOLVE FAILED: {(r.stdout + r.stderr)[-400:]}")
        return False
    x = read_sms_row(base + "_x.sms", M.shape[1], ell)
    # exact verification (scipy matvec) BEFORE any save
    resid = (M.dot(x) - b) % ell
    ok_mv = not resid.any()
    ok_cp, info = check_primal(S, Clo, Col, x)
    p(f"[{name}] verify: matvec_zero={ok_mv}  check_primal={ok_cp}  {info}")
    for suf in ("_MT.sms", "_b.sms", "_x.sms"):
        try: os.unlink(base + suf)
        except OSError: pass
    if not (ok_mv and ok_cp):
        p(f"[{name}] VERIFICATION FAILED -- not saving")
        return False
    fn = os.path.join(OUT, f"{name}_primal.npz")
    np.savez_compressed(fn, kind=np.array(["primal_member"]), ell=ell,
                        witness=x.astype(np.int8), A=S["A"].astype(np.int8))
    p(f"[{name}] SAVED {fn}  (walltime solve {dt:.1f}s)")
    return True

def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    threads = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    jobs = []
    if which in ("all", "K6"):
        nv, E = K6(); jobs.append(("K6_member", nv, E))
    if which in ("all", "del_oct"):
        nv, E = del_edge(K1222, 10); jobs.append(("K1222_del_oct", nv, E))  # edge (2,4)
    if which in ("all", "del_apx"):
        nv, E = del_edge(K1222, 4); jobs.append(("K1222_del_apx", nv, E))   # edge (0,6)
    okall = True
    for nm, nv, E in jobs:
        okall &= solve_and_save(nm, nv, E, ELL, threads)
    p("\nRESULT:", "ALL OK" if okall else "SOME FAILED")

if __name__ == "__main__":
    main()
