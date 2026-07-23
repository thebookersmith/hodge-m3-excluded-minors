# Nonmembership certificate for M(K₁,₂,₂,₂)

This directory contains a certificate that the graphic matroid M(K₁,₂,₂,₂) is a
nonmember of the class M₃ of Engel–de Gaay Fortman–Schreieder
(arXiv:2507.15704v3, Def. 8.3).

K₁,₂,₂,₂ is the complete four-partite graph with parts of sizes 1, 2, 2, 2 —
equivalently, K₇ with a perfect matching on six of its vertices deleted (the seventh
vertex remaining adjacent to all others). It has 7 vertices, 18 edges, and matroid
rank 6.

| object | graph6 | rank | corank | certificate |
|---|---|---:|---:|---|
| M(K₁,₂,₂,₂) | ``F]~vw`` | 6 | 12 | dual nonmembership certificate |

## What is certified

The witness `K1222_dual.npz` is a dual certificate of nonmembership: a left-kernel
vector of the Albanese membership system over F₃. It is given as explicit functions
a₁,…,a₆ on the 3¹² = 531,441 vertices of the system, together with a vector b on the
18 edges, satisfying for every edge s and every vertex w

    φ_s(w + g_s) − φ_s(w) = −b_s,        φ_s := Σ_k A[k,s]·a_k,

with Σ_s b_s = 1 ≠ 0. By the Farkas alternative this rules out membership:

**M(K₁,₂,₂,₂) is not in M₃.**

That is all this directory certifies. Minor-minimality is **not** established here:
showing that K₁,₂,₂,₂ is an excluded minor of M₃ would require, in addition, primal
membership certificates for each of its 36 single-element minors (18 deletions and 18
contractions), and those are not yet included.

## Verify

The witness uses the same `.npz` schema and dual-verification convention as the
certificates in `../m3_certificates/` (kind `dual_nonmember`; keys `A`, `witness`,
`ell`). A checker rebuilds the realization from `A`, recomputes the cycle-space
generators, and verifies the edge identity above by exact arithmetic modulo 3 over all
3¹² vertices, together with Σ_s b_s ≠ 0. Archived reference run: `K1222_CHECK_OUTPUT.txt`.

## Provenance

The dual was found in a computational search (2026-07-22); the certificate stands on
its own — its verification uses only the stored realization and witness, and does not
depend on the search that produced it.

Companion directories: `../m3_certificates/`, `../window_certificates/`,
`../l5_membership/`.
