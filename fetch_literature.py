"""Fetch hEDS/HSD-related PubMed literature and PubTator3 entity annotations.

Scoped to hEDS/HSD terminology specifically (not bare "Ehlers-Danlos syndrome")
to reduce genes from other EDS subtypes (classical, vascular) in the corpus.
Filtering by query terms cannot fully exclude them: shared terminology and
co-mentions still let other-subtype genes in.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path

import requests

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBTATOR_EXPORT_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson"
CACHE_DIR = Path(__file__).parent / "cache"
BATCH_SIZE = 100
REQUEST_DELAY = 1.0  # seconds - NCBI/PubTator courtesy pacing

QUERY = (
    '("hypermobile Ehlers-Danlos"[tiab] OR "Ehlers-Danlos syndrome, hypermobility type"[tiab] '
    'OR "EDS-HT"[tiab] OR "Ehlers-Danlos syndrome type III"[tiab] OR "Ehlers-Danlos syndrome type 3"[tiab] '
    'OR "joint hypermobility syndrome"[tiab] OR "hypermobility spectrum disorder"[tiab] '
    'OR "hypermobility spectrum disorders"[tiab] OR "generalized joint hypermobility"[tiab])'
)


def search_pmids(retmax=5000):
    params = {
        "db": "pubmed", "term": QUERY, "retmode": "json", "retmax": retmax,
        "tool": "heds-biomarker-discovery", "email": "darlinton.or21@gmail.com",
    }
    resp = requests.get(ESEARCH_URL, params=params, timeout=30)
    resp.raise_for_status()
    return resp.json()["esearchresult"]["idlist"]


def _batch_cache_path(batch):
    key = hashlib.sha1(",".join(batch).encode()).hexdigest()[:16]
    return CACHE_DIR / f"pubtator_{key}.json"


def _pubtator_get(pmids):
    """One export call for a list of PMIDs. Returns the 'PubTator3' document list,
    or None if PubTator3 rejects the whole request (e.g. an unindexed/too-recent PMID
    in the batch - it 400s the entire batch rather than skipping the bad one)."""
    resp = requests.get(PUBTATOR_EXPORT_URL, params={"pmids": ",".join(pmids)}, timeout=60)
    if resp.status_code != 200:
        return None
    return resp.json().get("PubTator3", [])


def fetch_documents(pmids, skip_cache=False):
    """Pull PubTator3 entity-annotated documents for a list of PMIDs, batched and cached."""
    CACHE_DIR.mkdir(exist_ok=True)
    documents = []
    for i in range(0, len(pmids), BATCH_SIZE):
        batch = pmids[i:i + BATCH_SIZE]
        cache_file = _batch_cache_path(batch)
        if not skip_cache and cache_file.exists():
            documents.extend(json.loads(cache_file.read_text(encoding="utf-8")))
            continue

        docs = _pubtator_get(batch)
        if docs is None:
            # Whole batch rejected - fall back to one-at-a-time so a single bad PMID
            # doesn't lose the other ~99 good ones in the batch.
            docs = []
            for pmid in batch:
                single = _pubtator_get([pmid])
                if single:
                    docs.extend(single)
                time.sleep(REQUEST_DELAY)
        else:
            time.sleep(REQUEST_DELAY)

        cache_file.write_text(json.dumps(docs), encoding="utf-8")
        documents.extend(docs)
    return documents


def fetch_corpus(retmax=5000, skip_cache=False):
    """End-to-end: search PubMed, pull annotated documents. Returns BioC document dicts."""
    pmids = search_pmids(retmax=retmax)
    return fetch_documents(pmids, skip_cache=skip_cache)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retmax", type=int, default=5000)
    parser.add_argument("--skip-cache", action="store_true")
    parser.add_argument("--limit", type=int, default=None, help="Only fetch first N PMIDs (smoke test)")
    args = parser.parse_args()

    pmids = search_pmids(retmax=args.retmax)
    print(f"Found {len(pmids)} PMIDs for hEDS/HSD query")
    if args.limit:
        pmids = pmids[:args.limit]
    docs = fetch_documents(pmids, skip_cache=args.skip_cache)
    print(f"Fetched {len(docs)} annotated documents")
