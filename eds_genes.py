"""Frozen differential-diagnosis annotations, not hEDS ground truth.

GeneReviews NBK1279, Table 1 (chapter updated 2024-02-22):
https://www.ncbi.nlm.nih.gov/books/NBK1279/
Human Entrez IDs identify genes independently of PubTator display-name aliases.
Includes brittle cornea syndrome as listed in that differential table. Membership
does not rule out a role in hEDS; absence means only 'not in this curated set'.
This present-day set is a retrospective baseline, not a historical knowledge set.
"""

OTHER_EDS_GENES = {
    "Gene:9509": {"symbol": "ADAMTS2", "subtype": "dermatosparaxis EDS"},
    "Gene:165": {"symbol": "AEBP1", "subtype": "classical-like EDS type 2"},
    "Gene:126792": {"symbol": "B3GALT6", "subtype": "spondylodysplastic EDS"},
    "Gene:11285": {"symbol": "B4GALT7", "subtype": "spondylodysplastic EDS"},
    "Gene:715": {"symbol": "C1R", "subtype": "periodontal EDS"},
    "Gene:716": {"symbol": "C1S", "subtype": "periodontal EDS"},
    "Gene:113189": {"symbol": "CHST14", "subtype": "musculocontractural EDS"},
    "Gene:1277": {"symbol": "COL1A1", "subtype": "classical/vascular/arthrochalasia EDS; COL1 overlap"},
    "Gene:1278": {"symbol": "COL1A2", "subtype": "arthrochalasia/cardiac-valvular EDS; COL1 overlap"},
    "Gene:1281": {"symbol": "COL3A1", "subtype": "vascular EDS"},
    "Gene:1289": {"symbol": "COL5A1", "subtype": "classical EDS"},
    "Gene:1290": {"symbol": "COL5A2", "subtype": "classical EDS"},
    "Gene:1303": {"symbol": "COL12A1", "subtype": "myopathic EDS"},
    "Gene:29940": {"symbol": "DSE", "subtype": "musculocontractural EDS"},
    "Gene:55033": {"symbol": "FKBP14", "subtype": "kyphoscoliotic EDS"},
    "Gene:5351": {"symbol": "PLOD1", "subtype": "kyphoscoliotic EDS"},
    "Gene:11107": {"symbol": "PRDM5", "subtype": "brittle cornea syndrome"},
    "Gene:91252": {"symbol": "SLC39A13", "subtype": "spondylodysplastic EDS"},
    "Gene:7148": {"symbol": "TNXB", "subtype": "classical-like EDS"},
    "Gene:84627": {"symbol": "ZNF469", "subtype": "brittle cornea syndrome"},
}
