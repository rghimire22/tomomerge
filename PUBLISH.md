# Putting TomoMerge online

Five minutes to a real `https://` link you can open anywhere, send to a
collaborator, and cite in the paper.

Everything served is already built into `docs/`, generated from
`merge_tomo.html` by `make_site.py`. That script also verifies the published
page runs **byte-identical code** to the tested one, and refuses to write
otherwise — so the live tool is always the tool the test suite passed.

---

## 1. Make a GitHub repository

Sign in at [github.com](https://github.com) → **New repository**.

- Name: `tomomerge` (this becomes part of your URL)
- **Public** — GitHub Pages needs this on a free account
- Do **not** add a README or .gitignore; the folder already has files

## 2. Push the folder

From a terminal in `Merging/`:

```bash
git init -b main
git add .
git commit -m "TomoMerge: merge overlapping tomography models"
git remote add origin https://github.com/YOUR-USERNAME/tomomerge.git
git push -u origin main
```

If `git` asks for a password, use a **personal access token**, not your account
password: GitHub → Settings → Developer settings → Personal access tokens →
Fine-grained → generate one with *Contents: read and write* on this repository.

## 3. Turn on Pages

Repository → **Settings** → **Pages** → under *Build and deployment*:

- Source: **Deploy from a branch**
- Branch: **main**, folder: **/docs**
- **Save**

Wait about a minute. The URL appears at the top of that page:

```
https://YOUR-USERNAME.github.io/tomomerge/
```

## 4. Fix the canonical URL

The first build used a placeholder. Now that you know the real address:

```bash
python make_site.py --url https://YOUR-USERNAME.github.io/tomomerge/
git add docs && git commit -m "Set canonical URL" && git push
```

That corrects the `<link rel="canonical">` and `sitemap.xml`, which is what
search engines use to decide which address is the real one.

---

## Being findable in a search engine

Publishing does not make you searchable — indexing does, and it takes days to
weeks. To help it along:

1. **Google Search Console** — [search.google.com/search-console](https://search.google.com/search-console),
   add `https://YOUR-USERNAME.github.io/tomomerge/` as a URL-prefix property,
   verify by the HTML-tag method (paste the tag into `make_site.py`'s `head`
   block and rebuild), then *URL Inspection* → **Request indexing**.
2. **Bing Webmaster Tools** — same idea, and it feeds DuckDuckGo.
3. **Get one real link to it.** Search engines weight inbound links heavily. The
   paper's Data and Code Availability section, your university profile page, or
   an AGU/SSA abstract will do more than any meta tag.

The page already carries a description, keywords, Open Graph tags for link
previews, and JSON-LD structured data marking it as a `SoftwareApplication`.

**Search terms it should eventually answer:** *merge overlapping tomography
models*, *combine two phase velocity maps*, *tomography mosaic seam*,
*subarray inversion merge*.

---

## What is and is not public

**Public:** the tool, its source (it is one HTML file — anyone can read it), and
whatever else you push to the repository.

**Not public, ever:** your velocity data. The page has no server and makes no
network requests; files you drop are read by the browser and stay on your
machine. That is worth stating on the page, and it is — in the footer.

**Think before pushing** `vel02.fg.sa13(*)`, `vel02lht*`, `vel03*`, and
`paper/`. If the data or manuscript are not ready to be public, exclude them:

```bash
cat > .gitignore <<'EOF'
vel0*
merged.xyz
merged_diag.txt
merged_qc.png
paper/
tests/checkerboard/
tests/resolution/
tests/fixtures/
EOF
git rm -r --cached . && git add . && git commit -m "Exclude data and manuscript"
```

The tool, the tests and the documentation are enough for the link to be useful
and for the method to be reproducible.

---

## Citing it in the paper

Once the URL is live, the Data and Code Availability section can read:

> The implementation is available at `https://YOUR-USERNAME.github.io/tomomerge/`
> and in the repository `https://github.com/YOUR-USERNAME/tomomerge`.

For a permanent, citable version with a DOI — which reviewers increasingly expect
and which a moving `main` branch cannot provide — connect the repository to
[Zenodo](https://zenodo.org), then cut a GitHub **release**. Zenodo archives that
exact snapshot and mints a DOI for it.

---

## Alternatives, if GitHub is inconvenient

| Option | Notes |
|---|---|
| **Netlify Drop** — [app.netlify.com/drop](https://app.netlify.com/drop) | Drag the `docs/` folder onto the page. A URL in about ten seconds, no account needed to start. Easiest possible route. |
| **Cloudflare Pages** | Similar, generous free tier, custom domains. |
| **University web space** | Many departments give you `people.uh.edu/~user/`. Copy `docs/index.html` there. Often the most durable address. |
| **Just the file** | `merge_tomo.html` works offline by double-click and can be emailed. No URL, but no hosting either. |

## Running it locally over HTTP

Not needed — the page works fine from `file://` — but if you want a local URL:

```bash
cd docs
python -m http.server 8000
```

then open `http://localhost:8000/`. Visible only on your own machine.
