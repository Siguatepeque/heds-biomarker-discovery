"""Annotate ABC-model candidate genes with weak plausibility context (no ranking effect).

Reference symbols (KLK15, ACKR3, SLC39A13, MIA3, TNXB) are literature research
anchors, not ground truth, not causal claims, and not diagnostic criteria: they
mark where recent hEDS-adjacent reports cluster, and any locus-to-gene link is
provisional. OTHER_EDS_GENES in eds_genes.py is a curated
annotation set for other EDS subtypes used only for differential-diagnosis
context; absence from it (None) means not in the curated set, not biological
exclusion, and flagged genes are never removed for carrying a subtype flag.

Each annotation is heuristic, unfrozen/live, and unscored when its source is
unavailable: a failed lookup returns None (unknown), never a biological
False/0. STRING v12.0 here is the default interaction_partners endpoint
(functional associations including text-mining at score >= 700), i.e.
non-independent literature-adjacent context, not independent validation and
not a physical-only network. The ClinVar check is a limited keyword match over
at most 5 esummary records, i.e. keyword hits only, not established variant
evidence. GTEx v8 bulk medians are proxy-tissue context only.

Row metadata records annotation_retrieved_at_utc (run timestamp),
string_version, gtex_dataset, and annotation_mode='live_unfrozen'; the
timestamp marks when live sources were queried, not raw-snapshot
reproducibility, since unfrozen sources can drift.
"""
import math
import time
from datetime import datetime, timezone

import requests

from eds_genes import OTHER_EDS_GENES

HEDS_REFERENCE_GENES = {"KLK15", "ACKR3", "SLC39A13", "MIA3", "TNXB"}
# Full differential context for STRING: every curated other-EDS symbol, not a
# classical/vascular subset. Derived from eds_genes.py so the flag and the
# annotation reference cannot drift apart.

ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
STRING_BASE = "https://version-12-0.string-db.org/api"
GTEX_BASE = "https://gtexportal.org/api/v2"
GTEX_TISSUES = [
    "Cells_Cultured_fibroblasts",  # collagen/ECM-producing cell type - best available proxy
    "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Artery_Aorta",
]
# gtex_v10 is listed as a valid dataset ID but the medianGeneExpression endpoint returns
# no rows for it (as of 2026) - only gtex_v8 actually has bulk-tissue expression data here.
GTEX_DATASET = "gtex_v8"
REQUEST_DELAY = 1.0
CLINVAR_RELEVANT_TERMS = {"ehlers-danlos", "connective tissue", "joint", "collagen", "hypermobility"}


def _json_body(resp):
    """Parsed JSON body, or None when the body is missing/malformed.

    Malformed, null, or wrong-shaped bodies are lookup failures (unknown), never
    negative biological evidence.
    """
    try:
        body = resp.json()
    except (ValueError, AttributeError, TypeError):
        return None
    return body if isinstance(body, (dict, list)) else None


def gene_symbol(entrez_gene_id):
    try:
        resp = requests.get(
            f"{ENTREZ_BASE}/esummary.fcgi",
            params={"db": "gene", "id": entrez_gene_id, "retmode": "json"}, timeout=30,
        )
    except requests.RequestException:
        return None
    if resp.status_code != 200:
        return None
    body = _json_body(resp)
    if not isinstance(body, dict):
        return None
    result = body.get("result")
    if not isinstance(result, dict):
        return None
    doc = result.get(str(entrez_gene_id), {})
    if not isinstance(doc, dict):
        return None
    name = doc.get("name")
    return name if isinstance(name, str) and name.strip() else None


def string_connected_to_reference(symbol):
    """True/False on a successful lookup, None when the lookup itself failed."""
    all_refs = HEDS_REFERENCE_GENES | {entry["symbol"] for entry in OTHER_EDS_GENES.values() if isinstance(entry, dict) and isinstance(entry.get("symbol"), str)}
    try:
        resp = requests.get(
            f"{STRING_BASE}/json/interaction_partners",
            params={"identifiers": symbol, "species": 9606, "required_score": 700,
                    "caller_identity": "heds-biomarker-discovery"},
            timeout=30,
        )
    except requests.RequestException:
        return None
    if resp.status_code != 200:
        return None
    partners = _json_body(resp)
    if not isinstance(partners, list):
        return None
    names = set()
    for row in partners:
        # A malformed row (or one with no usable partner name) makes the whole
        # response unknown: silently skipping it could hide a reference hit.
        if not isinstance(row, dict):
            return None
        name = row.get("preferredName_B")
        if not isinstance(name, str) or not name.strip():
            return None
        names.add(name)
    return bool(names & all_refs)


def gtex_expression(symbol):
    """Max measured median TPM across proxy tissues, or None when unknown.

    Only a returned measurement of median 0 establishes zero expression; an empty
    expression result (gene unseen, no rows) is NO DATA (None), not evidence of
    zero expression.
    """
    try:
        resp = requests.get(f"{GTEX_BASE}/reference/gene", params={"geneId": symbol}, timeout=30)
    except requests.RequestException:
        return None
    if resp.status_code != 200:
        return None
    body = _json_body(resp)
    if not isinstance(body, dict):
        return None
    data = body.get("data")
    if not isinstance(data, list) or not data or not isinstance(data[0], dict):
        return None
    gencode_id = data[0].get("gencodeId")
    if not gencode_id:
        return None
    try:
        resp = requests.get(
            f"{GTEX_BASE}/expression/medianGeneExpression",
            params={"gencodeId": gencode_id, "datasetId": GTEX_DATASET, "tissueSiteDetailId": GTEX_TISSUES},
            timeout=30,
        )
    except requests.RequestException:
        return None
    if resp.status_code != 200:
        return None
    expr = _json_body(resp)
    if not isinstance(expr, dict):
        return None
    rows = expr.get("data", [])
    if not isinstance(rows, list):
        return None
    values = []
    for row in rows:
        # Any unusable row makes the result unknown rather than silently
        # narrowing the claimed max over the remaining rows. Only finite,
        # non-negative numbers are valid TPM medians (bools excluded).
        if not isinstance(row, dict):
            return None
        median = row.get("median")
        if (isinstance(median, bool) or not isinstance(median, (int, float))
                or not math.isfinite(median) or median < 0):
            return None
        values.append(median)
    return max(values) if values else None


def clinvar_relevant(symbol):
    """True/False on a successful lookup, None when the lookup itself failed. A
    limited keyword match over up to 5 esummary records, not a variant review."""
    try:
        resp = requests.get(
            f"{ENTREZ_BASE}/esearch.fcgi",
            params={"db": "clinvar", "term": f"{symbol}[gene]", "retmode": "json", "retmax": 5}, timeout=30,
        )
    except requests.RequestException:
        return None
    if resp.status_code != 200:
        return None
    body = _json_body(resp)
    if not isinstance(body, dict):
        return None
    esearch = body.get("esearchresult")
    if not isinstance(esearch, dict):
        return None
    if "errorlist" in esearch or "ERROR" in body:
        return None
    ids = esearch.get("idlist", None)
    if not isinstance(ids, list):
        # A missing/non-list idlist means no actual search result - unknown,
        # not a negative. Only a present-but-empty idlist is a genuine False.
        return None
    if not ids:
        return False
    try:
        resp = requests.get(
            f"{ENTREZ_BASE}/esummary.fcgi", params={"db": "clinvar", "id": ",".join(str(i) for i in ids), "retmode": "json"}, timeout=30,
        )
    except requests.RequestException:
        return None
    if resp.status_code != 200:
        return None
    summary = _json_body(resp)
    if not isinstance(summary, dict):
        return None
    result = summary.get("result")
    # Keyword matching needs actual records: a missing/null/non-dict/empty result
    # is unknown, and matching against an error string must never read as biology.
    if not isinstance(result, dict) or not result:
        return None
    records = [result.get(str(identifier)) for identifier in ids]
    if any(not isinstance(record, dict) or not record or "error" in record
           for record in records):
        return None
    text = str(records).lower()
    return any(term in text for term in CLINVAR_RELEVANT_TERMS)


def validate(long_list_df, top_n=20):
    """Annotate each long-list candidate with STRING/GTEx/ClinVar context.

    Annotations never affect rank: rows are sorted by descending
    adamic_adar_score with a deterministic node_id tie-break. Unknown lookups
    (None) are marked unknown (annotation only) in the rationale. Returns the
    top_n rows. Unresolved symbols and other-EDS-subtype flags are kept, never
    dropped: other_eds_subtype is a curation annotation only (None means not
    in the curated set, not biological exclusion).
    """
    retrieved_at = datetime.now(timezone.utc).isoformat()
    rows = []
    for _, row in long_list_df.iterrows():
        node_id = row["node_id"]
        subtype = OTHER_EDS_GENES.get(node_id, {}).get("subtype")
        entrez_id = node_id.split(":", 1)[1]
        symbol = gene_symbol(entrez_id)
        time.sleep(REQUEST_DELAY)
        if not symbol:
            # Never drop a candidate: keep the row with a null symbol and all
            # sources unknown, and skip protein API queries without a symbol.
            name = row["name"] if "name" in row and row["name"] else node_id
            rows.append({
                "gene_symbol": None, "node_id": node_id,
                "adamic_adar_score": row["adamic_adar_score"],
                "string_connected": None, "gtex_max_median_tpm": None,
                "clinvar_relevant": None, "other_eds_subtype": subtype,
                "annotation_retrieved_at_utc": retrieved_at,
                "string_version": "12.0", "gtex_dataset": GTEX_DATASET,
                "annotation_mode": "live_unfrozen",
                "rationale": (
                    f"{node_id} ({name}): gene symbol unresolved, "
                    f"literature bridge score {row['adamic_adar_score']:.2f}; "
                    "STRING/GTEx/ClinVar unknown (annotation only, does not affect rank); "
                    f"other-EDS subtype flag (curation annotation only, null means not in curated set): {subtype}."
                ),
            })
            continue

        string_hit = string_connected_to_reference(symbol)
        time.sleep(REQUEST_DELAY)
        expression = gtex_expression(symbol)
        time.sleep(REQUEST_DELAY)
        clinvar_hit = clinvar_relevant(symbol)
        time.sleep(REQUEST_DELAY)

        string_txt = str(string_hit) if string_hit is not None else "unknown (annotation only, does not affect rank)"
        expression_txt = f"{expression:.1f}" if expression is not None else "unknown (annotation only, does not affect rank)"
        clinvar_txt = str(clinvar_hit) if clinvar_hit is not None else "unknown (annotation only, does not affect rank)"
        rationale = (
            f"{symbol}: literature bridge score {row['adamic_adar_score']:.2f}; "
            "STRING v12.0 functional (incl. text-mining) context, not independent validation "
            f"- connected to a reference anchor set: {string_txt}; "
            f"connective-tissue expression (max median TPM across proxies): {expression_txt}; "
            "ClinVar keyword hits among <=5 records (not established variant evidence) "
            f": {clinvar_txt}; "
            f"other-EDS subtype flag (curation annotation only, null means not in curated set): {subtype}."
        )
        rows.append({
            "gene_symbol": symbol, "node_id": node_id,
            "adamic_adar_score": row["adamic_adar_score"],
            "string_connected": string_hit, "gtex_max_median_tpm": expression,
            "clinvar_relevant": clinvar_hit, "other_eds_subtype": subtype,
            "annotation_retrieved_at_utc": retrieved_at,
            "string_version": "12.0", "gtex_dataset": GTEX_DATASET,
            "annotation_mode": "live_unfrozen",
            "rationale": rationale,
        })
    rows.sort(key=lambda r: (-r["adamic_adar_score"], r["node_id"]))
    return rows[:top_n]
