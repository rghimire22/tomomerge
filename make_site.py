#!/usr/bin/env python3
"""
Build the publishable web version of the tool into docs/, for GitHub Pages.

    python make_site.py                       # placeholder URL
    python make_site.py --url https://riwaj.github.io/tomomerge/

docs/index.html is generated from merge_tomo.html, not maintained beside it, so
the published page cannot drift from the tested one.  The script proves that by
comparing the two <script> blocks byte-for-byte after writing, and refuses to
finish if they differ.

Why docs/ : GitHub Pages can serve a repository's docs/ folder directly, so one
repository holds the code, the tests and the paper while also serving the app.
"""
import argparse
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "merge_tomo.html")
DOCS = os.path.join(HERE, "docs")

AUTHOR = "Riwaj Ghimire"
AFFIL = "University of Houston"
TITLE = "TomoMerge — merge overlapping seismic tomography models"
DESC = ("Free browser tool that merges two overlapping seismic tomography maps "
        "(phase velocity, group velocity or any lon/lat/value grid) with a "
        "Gaussian distance-weighted blend and an explicit criterion for the "
        "hand-over width. Runs entirely in your browser; files are never "
        "uploaded.")
KEYWORDS = ("seismic tomography, surface wave tomography, phase velocity, model "
            "merging, mosaicking, ambient noise, subarray, GMT, xyz grid, "
            "checkerboard resolution test")


def scripts(text):
    """
    Executable script blocks only.

    The published page carries a JSON-LD block for search engines, which is data
    and not code.  Comparing raw <script> tags would count it and report a
    difference that does not exist, so the type attribute is checked: no type, or
    a JavaScript type, means code.
    """
    out = []
    for attrs, body in re.findall(r"<script([^>]*)>([\s\S]*?)</script>", text):
        m = re.search(r'type\s*=\s*["\']([^"\']+)', attrs)
        if m and "javascript" not in m.group(1).lower() and m.group(1) != "module":
            continue
        out.append(body)
    return out


def build(url):
    src = open(SRC, encoding="utf-8").read()
    canonical = url.rstrip("/") + "/"

    head = f'''<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{TITLE}</title>
<meta name="description" content="{DESC}">
<meta name="keywords" content="{KEYWORDS}">
<meta name="author" content="{AUTHOR}, {AFFIL}">
<link rel="canonical" href="{canonical}">
<meta property="og:type" content="website">
<meta property="og:title" content="{TITLE}">
<meta property="og:description" content="{DESC}">
<meta property="og:url" content="{canonical}">
<meta name="twitter:card" content="summary">
<meta name="twitter:title" content="{TITLE}">
<meta name="twitter:description" content="{DESC}">
<link rel="icon" href="data:image/svg+xml,\
%3Csvg xmlns=&#39;http://www.w3.org/2000/svg&#39; viewBox=&#39;0 0 40 40&#39;%3E\
%3Crect x=&#39;2&#39; y=&#39;9&#39; width=&#39;21&#39; height=&#39;21&#39; rx=&#39;4&#39; fill=&#39;%232f6fd0&#39;/%3E\
%3Crect x=&#39;17&#39; y=&#39;9&#39; width=&#39;21&#39; height=&#39;21&#39; rx=&#39;4&#39; fill=&#39;%23d9730d&#39; opacity=&#39;.8&#39;/%3E\
%3C/svg%3E">
<script type="application/ld+json">
{{"@context":"https://schema.org","@type":"SoftwareApplication",
 "name":"TomoMerge","applicationCategory":"ScienceApplication",
 "operatingSystem":"Any browser","url":"{canonical}",
 "description":"{DESC}",
 "author":{{"@type":"Person","name":"{AUTHOR}","affiliation":"{AFFIL}"}},
 "offers":{{"@type":"Offer","price":"0","priceCurrency":"USD"}}}}
</script>'''

    # replace the source page's minimal head tags with the published set
    out = src.replace('<meta name="viewport" content="width=device-width,initial-scale=1">\n'
                      '<title>Tomographic Merge</title>', head, 1)
    if out == src:
        # the source head changed; fall back to inserting after the charset line
        out = src.replace('<meta charset="utf-8">',
                          '<meta charset="utf-8">\n' + head, 1)
        out = re.sub(r"<title>Tomographic Merge</title>\n?", "", out, count=1)

    os.makedirs(DOCS, exist_ok=True)
    with open(os.path.join(DOCS, "index.html"), "w", encoding="utf-8") as fh:
        fh.write(out)

    # Pages must serve the files as they are, not run them through Jekyll
    open(os.path.join(DOCS, ".nojekyll"), "w").close()

    with open(os.path.join(DOCS, "robots.txt"), "w", encoding="utf-8") as fh:
        fh.write(f"User-agent: *\nAllow: /\nSitemap: {canonical}sitemap.xml\n")
    with open(os.path.join(DOCS, "sitemap.xml"), "w", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n'
                 '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n'
                 f'  <url><loc>{canonical}</loc><changefreq>monthly</changefreq>'
                 '<priority>1.0</priority></url>\n</urlset>\n')

    for extra in ("palette_chart.png",):
        p = os.path.join(HERE, extra)
        if os.path.isfile(p):
            shutil.copyfile(p, os.path.join(DOCS, extra))

    # the published page must be the tested page
    a, b = scripts(src), scripts(out)
    if len(a) != 2 or a != b:
        sys.exit("ERROR: the published page's code differs from merge_tomo.html "
                 f"({len(a)} vs {len(b)} script blocks, identical={a == b}). "
                 "Not shipping a page that has not been tested.")

    size = os.path.getsize(os.path.join(DOCS, "index.html"))/1024
    print(f"docs/index.html written ({size:.0f} KB)")
    print("  code verified byte-identical to merge_tomo.html")
    print(f"  canonical URL: {canonical}")
    if "example" in url or "your-" in url:
        print("\n  ! placeholder URL. Re-run with --url once you know your Pages "
              "address,\n    so the canonical link and sitemap are right.")
    print("\nNext: see PUBLISH.md")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default="https://your-username.github.io/tomomerge/",
                    help="the GitHub Pages URL this will be served from")
    build(ap.parse_args().url)
