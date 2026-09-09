"""Pydantic request/response models for the HeartVar API.

Extracted verbatim from ``backend/app.py`` as the first step of splitting the
FastAPI god-module into a package. This module must NOT import from
``backend.app`` (it depends only on stdlib + pydantic) to avoid an import cycle.
"""

import re
from typing import Literal

from pydantic import BaseModel, Field, field_validator


_GENE_SYMBOL_RE = re.compile(r"^[A-Za-z0-9.\-]+$")


class CurationRequest(BaseModel):
    gene: str = Field("", max_length=64)
    hgvs_c: str = Field(..., max_length=256)
    hpo: str = Field("", max_length=2000)

    @field_validator("gene")
    @classmethod
    def _validate_gene_symbol(cls, v: str) -> str:
        v = (v or "").strip()
        if v and not _GENE_SYMBOL_RE.match(v):
            raise ValueError("gene must contain only letters, digits, '.' and '-'")
        return v
    inheritance: str = Field("", max_length=128)
    family: str = Field("", max_length=5000)
    segregation_context: str = Field("", max_length=5000)
    zygosity: str = Field("", max_length=32)
    inheritance_input: str = Field("", max_length=32)
    proband_sex: str = Field("", max_length=32)
    trio_status: str = Field("", max_length=32)
    denovo_status: str = Field("", max_length=48)
    seg_affected_carriers: int = Field(0, ge=0, le=100000)
    seg_affected_noncarriers: int = Field(0, ge=0, le=100000)
    seg_meioses: int = Field(0, ge=0, le=100000)
    in_trans_pathogenic: str = Field("", max_length=8)
    denovo_confirmed_count: int = Field(0, ge=0, le=1000)
    denovo_unconfirmed_count: int = Field(0, ge=0, le=1000)
    alt_cause_present: str = Field("", max_length=8)
    alt_cause_detail: str = Field("", max_length=1000)
    amino_acid: str = Field("", max_length=64)
    genome_build: str | None = Field(None, max_length=16)
    ai_mode: str = Field("none", max_length=16)
    enable_erepo: bool = True


class ChatTurn(BaseModel):
    """One prior exchange in the "Ask about this variant" chat history."""

    role: str = Field("", max_length=8)
    text: str = Field("", max_length=4000)


class ChatRequest(BaseModel):
    """Request body for the server-side chatbot (POST /api/chat).

    The browser sends the evidence ``context`` it already has loaded (the
    rendered variant/criteria/classification summary) plus the user's
    ``question`` and a little recent ``history``. The server owns the scoping
    system prompt (see ``claude.chat_reply``); nothing here can override the
    "only this variant / refuse off-topic" guardrail. Caps mirror the curation
    request — pure abuse protection, well above any real chat turn.
    """

    context: str = Field("", max_length=20000)
    question: str = Field(..., min_length=1, max_length=2000)
    history: list[ChatTurn] = Field(default_factory=list, max_length=20)
    mode: Literal["chat", "report"] = "chat"


class CriterionOverride(BaseModel):
    """One criterion as the Criteria tab posts it back for re-scoring.

    Mirrors the shape the browser was handed, so the client can send its edited
    list straight back without reshaping it.

    ``curator_override`` marks a line the human changed. It is what lets a
    deliberately-enabled PP5/BP6 score points while the engine's own PP5/BP6
    stay at zero (SVI 2018) — see ``acmg.tiers._points_for``. Nothing sets it
    server-side.
    """

    code: str = Field("", max_length=64)
    status: str = Field("not_met", max_length=32)
    criteria_strength: str | None = Field(None, max_length=64)
    direction: str = Field("", max_length=16)
    curator_override: bool = False


class InterpretRequest(BaseModel):
    """Request body for POST /api/curate/interpret.

    One field: the token POST /api/curate/stream handed back with an
    evidence-only result. Everything the AI half of the pipeline needs — the
    gathered evidence, the resolved gene/HGVS, the deterministic criteria, the
    original request — stays SERVER-side under that token (see
    backend/interpret_cache.py). The client cannot supply evidence, so it cannot
    shape the prompt or the ACMG score beyond the query it already ran.
    """

    token: str = Field(..., min_length=16, max_length=128)


class RescoreRequest(BaseModel):
    """Request body for POST /api/acmg/rescore.

    ``gene``/``inheritance_input`` are not decoration: the cross-criterion
    exclusion pass uses them for the CSpec applicability guard on the recessive
    PS4 rule, so omitting them would let an override reach a combination the
    automatic call would have demoted.

    ``max_length`` on the list is abuse protection in the same spirit as the
    caps on CurationRequest — the real ACMG set is ~28 criteria.
    """

    gene: str = Field("", max_length=64)
    inheritance_input: str = Field("", max_length=32)
    criteria: list[CriterionOverride] = Field(..., max_length=200)
