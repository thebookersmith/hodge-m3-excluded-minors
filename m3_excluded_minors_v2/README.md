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

Layout: `data/manifest.json` (claims + orbit accounting), `witnesses/*.npz`,
`m3_check.py`, `CHECK_OUTPUT.txt`. Companion summary: `..\fivefinds_summary.md`.
Predecessor bundle (finds #1–#2 only, xh-format witnesses): `..\newminor_certificate\`
— its witnesses were converted to the unified format here, with every conversion
re-verified at build time (`code\mfx_bundle_build.py`).

Companion bundle for claims BEYOND materializable scale (K₃,₉ ∉ M₅ at 5¹⁶, K₈ ∉ M₃ at
3²¹, verified through the socle-window ring identity with its own self-contained
checker): `..\window_certificates\`.
