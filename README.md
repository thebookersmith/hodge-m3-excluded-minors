# Five excluded minors for M₃: certificates and verification code

This repository accompanies the note **"New excluded minors for the class M₃ of
regular matroids"** by Booker Smith (`newminor_note.pdf`).

Engel, de Gaay Fortman, and Schreieder introduced the minor-closed classes M_ℓ in
*Matroids and the integral Hodge conjecture for abelian varieties*
([arXiv:2507.15704](https://arxiv.org/abs/2507.15704)) and posed the excluded-minor
characterization of M_ℓ as Problem 8.8.

The note gives five graphic excluded minors for M₃, of ranks 8, 9, 9, 9, and 10.
Together with the previously known example M(K₃,₅) (rank 7), this gives excluded
minors in every rank from 7 through 10. Corollary 2.10 of the recent
Engel–Schreieder preprint ([arXiv:2606.31894](https://arxiv.org/abs/2606.31894))
also shows that the five new examples are excluded minors for the larger
nonzero-profile class M̃₃.

Each computational claim comes with an explicit linear-algebra certificate. The
scripts below reconstruct the relevant membership systems from the stored
realizations and verify the certificates directly; verification does not depend on
the search or solver used to find the examples.

## Verification

Requirements: Python 3 and NumPy.

### Five excluded minors for M₃

```
cd m3_certificates
python3 m3_check.py
```

A successful run ends with `RESULT: PASS` (about one second, well under 1 GB of
memory). For each of the five graphs the checker verifies the identification of the
stored realization with the stated graph, a dual certificate of nonmembership in
M₃, and primal membership certificates for every single-element deletion and
contraction (individually, or through automorphism-orbit representatives). It also
verifies the K₃,₅ reference case.

### Nonmembership certificates beyond constructible scale

```
cd window_certificates
python3 window_check.py
```

A successful run ends with `RESULT: PASS`. This verifies **M(K₈) ∉ M₃** (whose full
membership system would have 3²¹ vertices) and **M(K₃,₉) ∉ M₅** (5¹⁶ vertices)
through identities in a truncated polynomial ring, without constructing the full
Albanese graphs. The K₃,₉ certificate is the p = 5 case underlying Theorem 1.7 of
[arXiv:2606.31894](https://arxiv.org/abs/2606.31894).

### Membership certificates at ℓ = 5

```
cd l5_membership
python3 l5_check.py
```

A successful run ends with `RESULT: PASS`. This verifies M(K₃,₅), M(G₁), M(G₂) ∈ M₅
by explicit solutions of the full 5⁸-vertex membership systems.

## Contents

- `newminor_note.pdf` — the accompanying note.
- `m3_certificates/` — certificates for the five excluded minors of M₃.
- `window_certificates/` — nonmembership certificates for systems too large to
  construct explicitly.
- `l5_membership/` — membership certificates at ℓ = 5.

## License

The verification code and certificate data are available under the MIT License
(see `LICENSE`). The note `newminor_note.pdf` is © 2026 Booker Smith and is not
included under that license.

## Contact

Booker Smith — bookers7@gmail.com
