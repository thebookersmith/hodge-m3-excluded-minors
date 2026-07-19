# Window certificates for nonmembership in M_ℓ

This directory contains explicit nonmembership certificates for two Albanese
membership systems that are too large to construct directly:

| claim | corank | full system | witness |
|---|---:|---|---|
| **M(K₃,₉) ∉ M₅** | 16 | 5¹⁶ = 152,587,890,625 vertices | `K39_ell5_windowD5_dual.npz` |
| **M(K₈) ∉ M₃** | 21 | 3²¹ = 10,460,353,203 vertices | `K8_ell3_windowD3_dual.npz` |

The K₃,₉ certificate gives the p = 5 nonmembership underlying Theorem 1.7 of
Engel–Schreieder (arXiv:2606.31894); the theorem's nonzero-profile form then follows
by their Corollary 2.10, K₃,₉ being biconnected. (M(K₈) ∉ M₃ also follows from
M(K₇) ∉ M₃ [arXiv:2512.04902, Thm 1.2] together with minor-closedness; it is included
here to demonstrate the method.)

## Verify

    python3 window_check.py        # requires Python 3 and NumPy

A successful run ends with `RESULT: PASS` (archived reference run: `CHECK_OUTPUT.txt`).

The checker embeds the edge lists, rebuilds each realization by row reduction and
requires it to match the one stored in the certificate, recomputes the cycle-space
generators, and verifies, monomial by monomial in F_ℓ[t₁..t_c]/(tᵢ^ℓ), the identity

  Σₖ A[k,s]·(x^{−g_s} − 1)·âₖ + b_s·t^socle = 0   for every edge s,   Σ_s b_s ≠ 0.

The comments at the top of `window_check.py` explain the correspondence between this
truncated-ring identity and a dual certificate for the full membership system, so
that nonmembership follows without constructing the ℓ^c vertices.

Separate computations, not included here, rule out certificates of depths 2, 3, and 4
for K₃,₉; the depth-5 certificate shown is therefore minimal among the window
certificates tested. Companion directories: `../m3_certificates/` and
`../l5_membership/`.
