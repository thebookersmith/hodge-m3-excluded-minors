# Membership certificates at ℓ = 5

This directory contains explicit certificates showing that M(K₃,₅), M(G₁), and
M(G₂) lie in M₅. These are the three M₃ excluded minors for which the accompanying
note reports computations at ℓ = 5.

| matroid | graph | (rank, edges) | full system | witness |
|---|---|---|---|---|
| **M(K₃,₅) ∈ M₅** | K₃,₅ | (7, 15) | 5⁸ = 390,625 vertices | `K35_ell5_member.npz` |
| **M(G₁) ∈ M₅** | ``H?zTbbo`` | (8, 16) | 390,625 vertices | `G1_ell5_member.npz` |
| **M(G₂) ∈ M₅** | ``J??FCpSJFw?`` | (10, 18) | 390,625 vertices | `G2_ell5_member.npz` |

## Verify

    python3 l5_check.py        # requires Python 3 and NumPy

A successful run ends with `RESULT: PASS` (archived reference run: `CHECK_OUTPUT.txt`).
For each matroid the checker first identifies the stored realization with the named
graph. It then reconstructs the full 5⁸-vertex membership system and verifies that
the stored primal certificate satisfies the closedness equations with color profile
identically 1.

The rank-9 excluded minors G₃, G₄, G₅ have not been tested at ℓ = 5 (as noted in the
accompanying note). Companion directories: `../m3_certificates/` and
`../window_certificates/`.
