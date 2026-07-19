# window_certificates — socle-window Farkas certificates beyond materializable scale

Machine-checkable certificates of NON-membership in the EGFS classes M_ℓ
(arXiv:2507.15704v3, Def 8.3) for two systems whose full Albanese graphs cannot be
written down:

| claim | corank | full-system size | witness |
|---|---|---|---|
| **K₃,₉ ∉ M₅** | 16 | 5¹⁶ = 152,587,890,625 vertices | depth-5 window dual (`K39_ell5_windowD5_dual.npz`) |
| **K₈ ∉ M₃** | 21 | 3²¹ = 10,460,353,203 vertices | depth-3 window dual (`K8_ell3_windowD3_dual.npz`) |

The K₃,₉ certificate gives the p = 5 non-membership underlying Theorem 1.7 of
Engel–Schreieder (arXiv:2606.31894); the theorem's full nonzero-profile form then
follows by their Cor. 2.10, K₃,₉ being biconnected. (K₈ ∉ M₃ also follows from
M(K₇) ∉ M₃ [arXiv:2512.04902, Thm 1.2] together with minor-closedness; it is
included here as a scale demonstration of the window method.)

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

Depth context: for K₃,₉ at ℓ = 5, the depths 2, 3, and 4 additionally carry explicit
infeasibility certificates (program archive; not part of this bundle), so depth 5 is
minimal for this witness. A general lower bound on witness depth is deferred to a
companion note.

Companion bundles: `../m3_excluded_minors_v2/` (five M₃ excluded minors, full-system
certificates) and `../l5_membership/` (ℓ = 5 membership witnesses).
