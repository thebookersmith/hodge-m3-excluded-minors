# Machine-checkable certificates: five new excluded minors for the class M₃ of regular matroids

This repository accompanies the note **"New excluded minors for the class M₃ of regular
matroids"** (Booker Smith, July 2026; `newminor_note.pdf` in this repository). The class
M₃ was introduced by Engel, de Gaay Fortman, and Schreieder in
[arXiv:2507.15704](https://arxiv.org/abs/2507.15704), where its excluded-minor
characterization is posed as Problem 8.8.

The note exhibits five graphic matroids — of ranks 8, 9, 9, 9, and 10 — that are excluded
minors for M₃ (and, combined with [arXiv:2606.31894](https://arxiv.org/abs/2606.31894),
Cor. 2.10, for the larger class M̃₃). Every claim is carried by explicit certificates
verifiable with the self-contained checkers below: no rank engine, solver, or other
software of ours sits in the trust path.

## Verify everything

Requirements: Python 3 and numpy. Nothing else.

```
cd m3_excluded_minors_v2
python3 m3_check.py
```

Expected final line: `RESULT: PASS` (≈1 second, well under 1 GB of memory). This rebuilds
each membership system from the embedded realizations and verifies, for each of the five
graphs: that the realization is a reduced oriented incidence matrix of the named graph
(graph6 strings embedded in the checker; explicit isomorphism verified edge-by-edge), a
left-kernel certificate of non-membership, and explicit solutions certifying membership
of every single-element deletion and contraction (minimality; coverage and Aut-orbits
verified inside the checker), plus the K₃,₅ anchor of [EGFS, Prop. 8.11].

```
cd window_certificates
python3 window_check.py
```

Expected final line: `RESULT: PASS`. This verifies two further certificates at scales
where the membership system cannot be materialized, via ring identities in
F_ℓ[t₁..t_c]/(tᵢ^ℓ): **M(K₈) ∉ M₃** (corank 21; the full system would have 3²¹ vertices)
and **M(K₃,₉) ∉ M₅** (corank 16; 5¹⁶ vertices) — the non-membership underlying the
p = 5 case of [arXiv:2606.31894], Theorem 1.7.

```
cd l5_membership
python3 l5_check.py
```

Expected final line: `RESULT: PASS`. This verifies the ℓ = 5 **membership** witnesses of
the note's §4: M(K₃,₅), M(G₁), M(G₂) ∈ M₅, each by an explicit solution of the full
5⁸-vertex membership system (`VERIFIED: MEMBER of M_5`).

## Contents

- `newminor_note.pdf` — the 7-page note (statements, certificates table, verification details).
- `m3_excluded_minors_v2/` — the five-find bundle: realizations, dual witnesses, per-minor
  membership witnesses, checker (`m3_check.py`), and the archived reference run
  (`CHECK_OUTPUT.txt`).
- `window_certificates/` — the beyond-scale certificates: window witnesses for K₈ (ℓ=3)
  and K₃,₉ (ℓ=5), checker (`window_check.py`), reference run, and a README with the
  depth-minimality context.
- `l5_membership/` — the ℓ = 5 membership witnesses (K₃,₅, G₁, G₂ ∈ M₅) with their own
  self-contained checker (`l5_check.py`) and reference run.

## License

Verification code and certificate data: MIT (see `LICENSE`). The note
`newminor_note.pdf` is © 2026 Booker Smith (distributed here for reference).

## Contact

Booker Smith — bookers7@gmail.com
