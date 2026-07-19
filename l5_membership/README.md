# l5_membership — full-system membership certificates at ℓ = 5

Machine-checkable certificates that the three known excluded minors of M₃ with
recorded ℓ = 5 data all lie in the class M₅ of Engel–de Gaay Fortman–Schreieder
(arXiv:2507.15704v3, Def 8.3):

| matroid | graph | (rank, edges) | full system | witness |
|---|---|---|---|---|
| **M(K₃,₅) ∈ M₅** | K₃,₅ | (7, 15) | 5⁸ = 390,625 vertices | `K35_ell5_member.npz` |
| **M(G₁) ∈ M₅** | `H?zTbbo` | (8, 16) | 390,625 vertices | `G1_ell5_member.npz` |
| **M(G₂) ∈ M₅** | `J??FCpSJFw?` | (10, 18) | 390,625 vertices | `G2_ell5_member.npz` |

## Verify

    python3 l5_check.py        # python3 + numpy only

Expected final line: `RESULT: PASS` (archived run: `CHECK_OUTPUT.txt`). For each
system the checker verifies (a) the realization in the witness file is a reduced
oriented incidence matrix of the named graph — the graph identities are embedded in
the checker and an explicit isomorphism is found and re-verified edge-by-edge — and
(b) the primal witness solves the full 5⁸-vertex Albanese membership system, rebuilt
from the realization alone, with color profile identically 1, reporting
`VERIFIED: MEMBER of M_5`.

The rank-9 excluded minors G₃, G₄, G₅ have not been tested at ℓ = 5 (as stated in
the accompanying note). Companion bundles: `../m3_excluded_minors_v2/` (the ℓ = 3
excluded-minor certificates) and `../window_certificates/` (non-membership beyond
materializable scale).
