"""Annotate ABC-model candidate genes with plausibility signals from public sources.

Reference gene sets reflect the current (2024-2025) hEDS-specific research landscape,
not the classical/vascular-EDS collagen genes those OTHER subtypes are defined by:
- KLK15 (kallikrein cluster): recurrent variant in WES of 200 hEDS patients, reproduces
  connective-tissue defects in mouse knock-ins (Norris Lab, PMC11213194, 2024).
- ACKR3, SLC39A13: reported genome-wide-significant loci from the first hEDS GWAS
  meta-analysis (Petrucci-Nelson et al., medRxiv 2025.09.19.25336146, preprint -
  still requiring replication), provisionally suggesting a neuroimmune/stromal
  model rather than a single collagen gene.
- MIA3: separately proposed 2024-2025 candidate, still unresolved.
- TNXB: long-studied partial-deficiency lead; real but explains only ~1% of cases, and
  serum tenascin-X failed as a screening test (PMC6820911) - kept as a weak legacy
  plausibility anchor, not a strong reference.

Each annotation is heuristic and unscored when its source is unavailable: a failed
lookup returns None (unknown), never a biological False/0. STRING here is the default
interaction_partners endpoint (functional + physical associations at score >= 700),
not a physical-only network and not an independent validation. The ClinVar check is
a limited keyword match over the first few esummary records, not a variant-level
review.
"""
import math
import time

import requests

HEDS_REFERENCE_GENES = {"KLK15", "ACKR3", "SLC39A13", "MIA3", "TNXB"}
# Genes that define OTHER EDS subtypes (classical, vascular) - relevant for differential
# diagnosis, not as hEDS-specific ground truth, so kept in a separate bucket.
DIFFERENTIAL_DIAGNOSIS_GENES = {"COL5A1", "COL5A2", "COL1A1", "COL1A2", "COL3A1"}

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
    all_refs = HEDS_REFERENCE_GENES | DIFFERENTIAL_DIAGNOSIS_GENES
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
    """For each long-list candidate, resolve its gene symbol and annotate with
    STRING/GTEx/ClinVar plausibility signals. Unknown lookups (None) contribute
    nothing to the heuristic composite and are marked unscored/unknown in the
    rationale. Returns the top_n rows by composite score, sorted descending."""
    rows = []
    for _, row in long_list_df.iterrows():
        entrez_id = row["node_id"].split(":", 1)[1]
        symbol = gene_symbol(entrez_id)
        time.sleep(REQUEST_DELAY)
        if not symbol:
            # Never drop a candidate: keep the row with a null symbol and all
            # sources unscored, and skip protein API queries without a symbol.
            name = row["name"] if "name" in row and row["name"] else row["node_id"]
            rows.append({
                "gene_symbol": None, "node_id": row["node_id"],
                "adamic_adar_score": row["adamic_adar_score"],
                "string_connected": None, "gtex_max_median_tpm": None,
                "clinvar_relevant": None, "composite_score": row["adamic_adar_score"],
                "rationale": (
                    f"{row['node_id']} ({name}): gene symbol unresolved, "
                    f"literature bridge score {row['adamic_adar_score']:.2f}; "
                    "STRING/GTEx/ClinVar unscored (unknown)."
                ),
            })
            continue

        string_hit = string_connected_to_reference(symbol)
        time.sleep(REQUEST_DELAY)
        expression = gtex_expression(symbol)
        time.sleep(REQUEST_DELAY)
        clinvar_hit = clinvar_relevant(symbol)
        time.sleep(REQUEST_DELAY)

        composite = (
            row["adamic_adar_score"]
            + (1 if string_hit is True else 0)
            + (min(expression / 10.0, 1.0) if expression is not None else 0)
            + (1 if clinvar_hit is True else 0)
        )
        string_txt = str(string_hit) if string_hit is not None else "unknown (source unscored)"
        expression_txt = f"{expression:.1f}" if expression is not None else "unknown (source unscored)"
        clinvar_txt = str(clinvar_hit) if clinvar_hit is not None else "unknown (source unscored)"
        rationale = (
            f"{symbol}: literature bridge score {row['adamic_adar_score']:.2f}; "
            f"STRING-connected to a known hEDS/EDS-family gene: {string_txt}; "
            f"connective-tissue expression (max median TPM across proxies): {expression_txt}; "
            f"ClinVar has a connective-tissue-relevant record: {clinvar_txt}."
        )
        rows.append({
            "gene_symbol": symbol, "node_id": row["node_id"],
            "adamic_adar_score": row["adamic_adar_score"],
            "string_connected": string_hit, "gtex_max_median_tpm": expression,
            "clinvar_relevant": clinvar_hit, "composite_score": composite,
            "rationale": rationale,
        })
    rows.sort(key=lambda r: r["composite_score"], reverse=True)
    return rows[:top_n]
