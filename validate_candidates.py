"""Cross-validate ABC-model candidate genes against independent public evidence.

Reference gene sets reflect the current (2024-2025) hEDS-specific research landscape,
not the classical/vascular-EDS collagen genes those OTHER subtypes are defined by:
- KLK15 (kallikrein cluster): recurrent variant in WES of 200 hEDS patients, reproduces
  connective-tissue defects in mouse knock-ins (Norris Lab, PMC11213194, 2024).
- ACKR3, SLC39A13: genome-wide-significant loci from the first hEDS GWAS meta-analysis
  (Petrucci-Nelson et al., medRxiv 2025.09.19.25336146), pointing to a neuroimmune/
  stromal model rather than a single collagen gene.
- MIA3: separately proposed 2024-2025 candidate, still unresolved.
- TNXB: long-studied partial-deficiency lead; real but explains only ~1% of cases, and
  serum tenascin-X failed as a screening test (PMC6820911) - kept as a weak legacy
  positive control, not a strong reference.
"""
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


def gene_symbol(entrez_gene_id):
    resp = requests.get(
        f"{ENTREZ_BASE}/esummary.fcgi",
        params={"db": "gene", "id": entrez_gene_id, "retmode": "json"}, timeout=30,
    )
    if resp.status_code != 200:
        return None
    doc = resp.json().get("result", {}).get(str(entrez_gene_id), {})
    return doc.get("name") or None


def string_connected_to_reference(symbol):
    all_refs = HEDS_REFERENCE_GENES | DIFFERENTIAL_DIAGNOSIS_GENES
    resp = requests.get(
        f"{STRING_BASE}/json/interaction_partners",
        params={"identifiers": symbol, "species": 9606, "required_score": 700,
                "caller_identity": "heds-biomarker-discovery"},
        timeout=30,
    )
    if resp.status_code != 200:
        return False
    partners = {row.get("preferredName_B") for row in resp.json()}
    return bool(partners & all_refs)


def gtex_expression(symbol):
    resp = requests.get(f"{GTEX_BASE}/reference/gene", params={"geneId": symbol}, timeout=30)
    if resp.status_code != 200 or not resp.json().get("data"):
        return 0.0
    gencode_id = resp.json()["data"][0]["gencodeId"]
    resp = requests.get(
        f"{GTEX_BASE}/expression/medianGeneExpression",
        params={"gencodeId": gencode_id, "datasetId": GTEX_DATASET, "tissueSiteDetailId": GTEX_TISSUES},
        timeout=30,
    )
    if resp.status_code != 200:
        return 0.0
    values = [row["median"] for row in resp.json().get("data", [])]
    return max(values) if values else 0.0


def clinvar_relevant(symbol):
    resp = requests.get(
        f"{ENTREZ_BASE}/esearch.fcgi",
        params={"db": "clinvar", "term": f"{symbol}[gene]", "retmode": "json", "retmax": 5}, timeout=30,
    )
    ids = resp.json().get("esearchresult", {}).get("idlist", []) if resp.status_code == 200 else []
    if not ids:
        return False
    resp = requests.get(
        f"{ENTREZ_BASE}/esummary.fcgi", params={"db": "clinvar", "id": ",".join(ids), "retmode": "json"}, timeout=30,
    )
    if resp.status_code != 200:
        return False
    text = str(resp.json().get("result", {})).lower()
    return any(term in text for term in CLINVAR_RELEVANT_TERMS)


def validate(long_list_df, top_n=20):
    """For each long-list candidate, resolve its gene symbol and check STRING/GTEx/ClinVar
    evidence. Returns the top_n rows by composite score, sorted descending."""
    rows = []
    for _, row in long_list_df.iterrows():
        entrez_id = row["node_id"].split(":", 1)[1]
        symbol = gene_symbol(entrez_id)
        time.sleep(REQUEST_DELAY)
        if not symbol:
            continue

        string_hit = string_connected_to_reference(symbol)
        time.sleep(REQUEST_DELAY)
        expression = gtex_expression(symbol)
        time.sleep(REQUEST_DELAY)
        clinvar_hit = clinvar_relevant(symbol)
        time.sleep(REQUEST_DELAY)

        composite = (
            row["adamic_adar_score"]
            + (1 if string_hit else 0)
            + min(expression / 10.0, 1.0)
            + (1 if clinvar_hit else 0)
        )
        rationale = (
            f"{symbol}: literature bridge score {row['adamic_adar_score']:.2f}; "
            f"STRING-connected to a known hEDS/EDS-family gene: {string_hit}; "
            f"connective-tissue expression (max median TPM across proxies): {expression:.1f}; "
            f"ClinVar has a connective-tissue-relevant record: {clinvar_hit}."
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
