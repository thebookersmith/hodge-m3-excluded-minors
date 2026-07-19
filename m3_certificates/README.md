# Certificates for five excluded minors of M₃

This directory contains certificates that five graphic matroids are excluded minors
for the class M₃ of Engel–de Gaay Fortman–Schreieder (arXiv:2507.15704v3, Def. 8.3,
Prop. 8.11), together with a certificate for the previously known example M(K₃,₅),
included as a reference case.

| example | graph6 | rank | certificates |
|---|---|---:|---|
| G₁ | ``H?zTbbo`` | 8 | dual nonmembership certificate; primal certificates for minor-orbit representatives |
| G₂ | ``J??FCpSJFw?`` | 10 | dual nonmembership certificate; primal certificates for minor-orbit representatives |
| G₃ | ``I?`FF_{F_`` | 9 | dual nonmembership certificate; primal certificates for all 34 single-element minors |
| G₄ | ``I?B@nRWN?`` | 9 | dual nonmembership certificate; primal certificates for all 34 single-element minors |
| G₅ | ``ICQf@pSF_`` | 9 | dual nonmembership certificate; primal certificates for all 34 single-element minors |
| M(K₃,₅) | — | 7 | dual nonmembership certificate; previously established in EGFS |

## Verify

    python3 m3_check.py        # requires Python 3 and NumPy

A successful run ends with `RESULT: PASS` (archived reference run: `CHECK_OUTPUT.txt`).

The checker reconstructs each Albanese membership system from the stored realization
and verifies the corresponding certificate by exact arithmetic modulo 3. It also
verifies the graph identifications and confirms that the stored minor certificates
cover every single-element deletion and contraction; for G₁ and G₂ it computes the
automorphism group and its edge orbits directly.

Layout: `data/manifest.json`, `witnesses/*.npz`, `m3_check.py`, `CHECK_OUTPUT.txt`.
Companion directories: `../window_certificates/` and `../l5_membership/`.
