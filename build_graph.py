"""Build a co-occurrence knowledge graph from PubTator3-annotated hEDS/HSD literature."""
import math

import networkx as nx

# Only biologically meaningful entity types. "Species" is deliberately dropped: PubTator
# tags nearly every clinical abstract with a "patients"/human Species annotation, which
# would otherwise become a degenerate hub node connecting almost every document.
KEPT_TYPES = {"Gene", "Disease", "Chemical"}


def doc_identity(doc):
    """Resolved identity of one BioC document: explicit ``pmid`` or ``id``.

    Numeric and string forms of the same identifier collapse (``123`` == ``"123"``).
    Returns None for anonymous documents (no usable ``pmid``/``id``).
    """
    for field in ("pmid", "id"):
        value = doc.get(field)
        if value is None:
            continue
        text = str(value).strip()
        if text:
            return text
    return None


def dedupe_documents(documents):
    """Collapse repeat documents to one each (first occurrence wins).

    Counts each PMID (or BioC ``id`` fallback) only once so a cached article stored
    twice cannot double edge weights or seed counts. Anonymous documents without a
    usable ``pmid``/``id`` (e.g. synthetic test fixtures) are keyed by object identity
    so distinct ones never collapse together.
    """
    seen = set()
    unique = []
    for doc in documents:
        ident = doc_identity(doc)
        key = ("doc", ident) if ident is not None else ("obj", id(doc))
        if key in seen:
            continue
        seen.add(key)
        unique.append(doc)
    return unique


def _node_id(infons):
    identifier = infons.get("identifier")
    entity_type = infons.get("type")
    if entity_type not in KEPT_TYPES or not identifier or identifier == "-":
        return None
    return f"{entity_type}:{identifier}"


def document_entities(document):
    """Extract {node_id: display_name} for one PubTator3 BioC document."""
    entities = {}
    for passage in document.get("passages", []):
        for ann in passage.get("annotations", []):
            node_id = _node_id(ann.get("infons", {}))
            if node_id:
                entities.setdefault(node_id, ann.get("text", node_id))
    return entities


def _pairs(items):
    for i, a in enumerate(items):
        for b in items[i + 1:]:
            yield a, b


def build_graph(documents):
    """Co-occurrence graph: nodes = bioconcepts, edges weighted by shared-abstract count,
    with a PMI attribute (log(P(a,b) / (P(a)*P(b)))) estimated from document frequency."""
    documents = dedupe_documents(documents)
    graph = nx.Graph()
    doc_entity_sets = [e for e in (document_entities(d) for d in documents) if e]
    doc_count = len(doc_entity_sets)
    node_doc_count = {}

    for entities in doc_entity_sets:
        for node_id, name in entities.items():
            graph.add_node(node_id, name=name, type=node_id.split(":", 1)[0])
            node_doc_count[node_id] = node_doc_count.get(node_id, 0) + 1
        for a, b in _pairs(list(entities)):
            if graph.has_edge(a, b):
                graph[a][b]["weight"] += 1
            else:
                graph.add_edge(a, b, weight=1)

    for a, b, data in graph.edges(data=True):
        p_a = node_doc_count[a] / doc_count
        p_b = node_doc_count[b] / doc_count
        p_ab = data["weight"] / doc_count
        data["pmi"] = math.log(p_ab / (p_a * p_b))

    return graph


def save_graph(graph, path):
    nx.write_graphml(graph, path)


def load_graph(path):
    return nx.read_graphml(path)
