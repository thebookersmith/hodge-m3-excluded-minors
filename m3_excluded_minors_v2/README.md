# m3_excluded_minors_v2 — unified certificate bundle

Machine-checkable certificates that **five graphic matroids are excluded minors for the
class M₃** of Engel–de Gaay Fortman–Schreieder (arXiv:2507.15704v3, Def 8.3, Prop 8.11),
plus the published K₃,₅ excluded minor as a pipeline anchor:

| find | graph6 | rank | scheme |
|---|---|---|---|
| find1 | `H?zTbbo` | 8 | dual + orbit-representative minor primals (S₄ accounting) |
| find2 | `J??FCpSJFw?` | 10 | dual + orbit-representative minor primals |
| find3 | `I?`FF_{F_` | 9 | dual + ALL 34 one-element minor primals |
| find4 | `I?B@nRWN?` | 9 | dual + ALL 34 one-element minor primals |
| find5 | `ICQf@pSF_` | 9 | dual + ALL 34 one-element minor primals |
| K35_anchor | K₃,₅ | 7 | dual only (published result; anchor) |

## Verify

    python3 m3_check.py        # requires python3 + numpy, nothing else

Expected final line: `RESULT: PASS` (full archived run: `CHECK_OUTPUT.txt`).

Every witness file in `witnesses/` is self-contained (`kind`, `ell`, `witness`, and the
integer realization `A`); the checker rebuilds each Albanese membership system from `A`
alone and verifies the witness with integer matrix–vector products mod 3 — no rank
engine, no external libraries, no trust in the software that *found* the witnesses.
Graph identity and minor coverage are verified inside the checker as well: the graph6
strings are embedded in `m3_check.py` itself, each parent realization is verified to be
a reduced oriented incidence matrix of the named graph by an explicit isomorphism, each
minor realization is matched against a realization *derived in the checker* from the
parent (identical row space, or explicit isomorphism of the minor graphs), and coverage
of **all** single-element deletions and contractions is confirmed — with the
automorphism groups and edge orbits of G₁/G₂ computed from scratch in the checker, not
quoted from the manifest.

Layout: `data/manifest.json` (claims + orbit accounting), `witnesses/*.npz`,
`m3_check.py`, `CHECK_OUTPUT.txt` (archived reference run).

Companion bundles: `../window_certificates/` (claims beyond materializable scale:
K₃,₉ ∉ M₅ at 5¹⁶, K₈ ∉ M₃ at 3²¹, via the socle-window ring identity) and
`../l5_membership/` (ℓ = 5 full-system membership witnesses).
