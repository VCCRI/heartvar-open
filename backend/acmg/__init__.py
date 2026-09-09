"""ACMG/AMP scoring engine.

This package holds the deterministic ACMG criterion scorer extracted from
``backend.app``. The modules are leaf-ordered and import only from each other
(plus clients / prompt / stdlib) — never from ``backend.app`` — so the engine
is importable standalone and free of an app-level import cycle.

  - ``constants``  : the ACMG constant layer (single source of truth, loaded
                     from ``static/acmg_constants.json`` + ``backend/data/*``).
  - ``tiers``      : point summation + tier classification.
  - ``hard_coded`` : ``compute_hard_coded_criteria`` and its helper closure
                     (the 16 deterministic criteria + ClinVar PP5/BP6 + the
                     cross-criterion exclusion pass + merge + family history),
                     so 18 of the 21 codes in ``hard_coded_criteria_codes``.
                     PM1, PP1 and BS4 are the other three and are derived by
                     ``no_ai`` — see ``compute_hard_coded_criteria``'s docstring.
  - ``no_ai``      : the no-key (AI-free) supplementary inference layer.

``backend.app`` re-exports every public/consumed symbol so existing importers
(`from backend.app import X`) keep working unchanged.
"""
