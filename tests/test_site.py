"""Static guards for docs/index.html. Stdlib only: python tests/test_site.py."""
import sys
from html.parser import HTMLParser
from pathlib import Path

PAGE = Path(__file__).resolve().parent.parent / "docs" / "index.html"
errors = []


class P(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids, self.dupes, self.anchors = [], set(), []
        self.h1 = self.main = self.lang = 0
        self.tables = []  # (has_caption, col_scopes, row_scopes)
        self.cur = None
        self.ext = []  # external runtime assets

    def handle_starttag(self, t, a):
        a = dict(a)
        if i := a.get("id"):
            if i in self.ids:
                self.dupes.add(i)
            else:
                self.ids.append(i)
        if t == "a" and (h := a.get("href", "")).startswith("#"):
            self.anchors.append(h[1:])
        if t == "h1":
            self.h1 += 1
        if t == "main":
            self.main += 1
        if t == "html":
            self.lang = a.get("lang", "")
        if t == "table":
            self.cur = [False, 0, 0]
            self.tables.append(self.cur)
        if t == "caption" and self.cur is not None:
            self.cur[0] = True
        if t == "th" and self.cur is not None:
            if a.get("scope") == "col":
                self.cur[1] += 1
            elif a.get("scope") == "row":
                self.cur[2] += 1
            else:
                errors.append("table header without row/col scope")
        if t in ("script", "img"):
            self.ext.append(t)
        if t == "link" and "stylesheet" in a.get("rel", "").split():
            self.ext.append("link-css")

    def handle_endtag(self, tag):
        if tag == "table":
            self.cur = None


p = P()
p.feed(PAGE.read_text(encoding="utf-8"))
if p.dupes:
    errors.append(f"duplicate ids: {sorted(p.dupes)}")
if missing := [a for a in p.anchors if a and a not in p.ids]:
    errors.append(f"anchors without target: {missing}")
if p.h1 != 1:
    errors.append(f"h1 count={p.h1}")
if p.main != 1:
    errors.append(f"main count={p.main}")
if p.lang != "en":
    errors.append(f"lang={p.lang!r}")
for i, (cap, col, row) in enumerate(p.tables):
    if not cap or not col or not row:
        errors.append(f"table{i}: caption={cap} colscopes={col} rowscopes={row}")
if p.ext:
    errors.append(f"runtime assets: {p.ext}")
print(f"ids={len(p.ids)} anchors={len(p.anchors)} tables={len(p.tables)}")
if errors:
    print("FAIL")
    for error in errors:
        print(" -", error)
    sys.exit(1)
print("PASS static site checks")
