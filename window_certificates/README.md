# window_certificates — socle-window Farkas certificates beyond materializable scale

Machine-checkable certificates of NON-membership in the EGFS classes M_ℓ
(arXiv:2507.15704v3, Def 8.3) for two systems whose full Albanese graphs cannot be
written down:

| claim | corank | full-system size | witness |
|---|---|---|---|
| **K₃,₉ ∉ M₅** (Engel–Schreieder Thm 1.7 instance, p=5) | 16 | 5¹⁶ = 152,587,890,625 vertices | depth-5 window dual (`K39_ell5_windowD5_dual.npz`) |
| **K₈ ∉ M₃** | 21 | 3²¹ = 10,460,353,203 vertices | depth-3 window dual (`K8_ell3_windowD3_dual.npz`) |

## Verify

    python3 window_check.py        # python3 + numpy only; no repo context needed

Expected final line: `RESULT: PASS` (archived run: `CHECK_OUTPUT.txt`).

`window_check.py` is fully self-contained: it embeds the edge lists, rebuilds each
realization by plain RREF and **requires** it to match the one stored in the witness
file, recomputes the cycle-space generators, and verifies the window ring identity

  Σₖ A[k,s]·(x^{−g_s} − 1)·âₖ + b_s·t^socle = 0  for every edge s,  Σ_s b_s ≠ 0

monomial-by-monomial in F_ℓ[t₁..t_c]/(tᵢ^ℓ). The header comments derive the
function↔element dictionary showing this finite identity IS a Farkas dual of the full
membership system (shift ↔ x^{−g}-multiplication; constant-1 ↔ the socle monomial;
window stability of the multiplication operator) — so non-membership follows without
ever materializing ℓ^c vertices.

Depth context (adjudicated theory): no Farkas dual of depth ≤ ℓ−1 exists for any graph
(so depth 5 is minimal at ℓ=5, depth 3 at ℓ=3); for K₃,₉ the depths 2,3,4 additionally
carry explicit GPU-verified infeasibility certificates (program archive,
`mfx_witnesses\` / `gpu_lab`). Companion full-system bundle (five M₃ excluded minors,
checkable directly): `..\m3_excluded_minors_v2\`.
