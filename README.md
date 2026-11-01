# Learn You a Tropical Geometry for Great Good! — website

The tutorial site for **minplus.xyz**, built with
[Zola](https://www.getzola.org/) (Rust static site generator) and deployed to
Cloudflare Pages. Idiom follows `cycle-five/cyclefive-site`.

## Stack

- **Zola 0.22.1** — static site generator (Tera templates + built-in libsass).
- No external theme. Custom templates in `templates/`, styles in `sass/`.
- No JS framework, no JS. The build is `zola build` plus
  [serve-the-source](https://github.com/cycle-five/serve-the-source), which
  writes each page's Markdown source beside its HTML for agents
  (`npm run build` runs both).
- Content: `content/chapters/` is a section sorted by `weight` (fixed
  reading order); each chapter is a page bundle (`chapter-NN/index.md`)
  with its figures in `figures/` and its full program beside it as a
  downloadable `.py`.
- Figures ship as **SVG** (rendered by the same matplotlib scripts that
  produced the study's PNGs), keeping the whole repo UTF-8 text.

## Local development

Requires Zola **0.22.1**.

```sh
zola serve     # live-reload dev server at http://127.0.0.1:1111
zola build     # HTML only, into ./public
npm run build  # what Pages runs: zola build, then the Markdown for agents
zola check     # validate internal links & markup
```

## Deploy — Cloudflare Pages

The Pages project (`minplus-site`, production branch `master`) is a
**direct-upload** project: deploys are pushes of the prebuilt `public/`,
not Cloudflare builds. That keeps the exact local build (Zola +
serve-the-source) as what ships.

From this repo, with the Cloudflare connector helper
(`~/workspace/custom-connectors/custom.cloudflare/bin/cf.py`):

```sh
npm run build   # zola build + serve-the-source into ./public
python3 ~/workspace/custom-connectors/custom.cloudflare/bin/cf.py \
  pages-deploy <account-id> minplus-site public master \
  functions/_middleware.js
```

The helper follows wrangler's upload protocol (upload-token,
check-missing, assets/upload, upsert-hashes, then a manifest-only
deployment) and adapts `functions/_middleware.js` into an advanced-mode
`_worker.js` for the deployment, so the Accept negotiation works the
same as it does in a git-connected build.

**The build must be `npm run build`, not `zola build`.** Plain
`zola build` deploys a site that looks fine, but it writes no `index.md`
files, so the middleware has nothing to serve and agents asking for
`text/markdown` silently get HTML. Check after a deploy:
`curl -sI -H 'Accept: text/markdown' https://minplus.xyz/chapters/chapter-01/`
should say `content-type: text/markdown`.

If the project is ever switched to git-connected builds instead, the
settings are: framework preset none, build command `npm run build`,
output dir `public`, and a `ZOLA_VERSION = 0.22.1` environment variable
for Production **and** Preview (the build image's default Zola is older).

## Marginalia (margin notes + print editions)

Chapters carry a handwriting-friendly layer, built by `tools/marginalia.py`
(the last step of `npm run build`):

- **Anchors.** Every top-level paragraph of a chapter is numbered within its
  section and stamped with an id (`s1-3-p2`) and a label (`§1.3 ¶2`).
- **Print edition.** Each chapter gets a `print/` page: Letter paper, a wide
  right margin for handwriting, anchor labels visible in the left margin.
  Print it, scribble the marginal notes you would make anyway as a reader,
  scan/photograph the pages.
- **Digitising.** Transcriptions go into the chapter bundle's
  `annotations.json`:

  ```json
  {"chapter": 1, "notes": [
    {"anchor": "s1-3-p2", "text": "the note, transcribed",
     "ink": "ink/note-01.png"}]}
  ```

  `ink` is an optional cropped image of the actual handwriting, kept beside
  the transcription (the build injects notes as real HTML — selectable,
  copyable, searchable — behind a CSS-only "author's margin notes" toggle;
  no JS). Ink images are binary, so they deploy with the site but are not
  part of the text-only GitHub snapshot push.

House rule (2026-11-01): chapter prose never dumps the full program inline.
The "whole chapter, as a program" section points at the downloadable `.py`
in the chapter bundle; the assertions in that file are the chapter's claims.

## Domains

- **`minplus.xyz`** (apex) and **`www.minplus.xyz`** — attached to the
  Pages project; the zone's apex and `www` CNAMEs point at
  `minplus-site.pages.dev` (proxied). Live once the zone is active
  (registrar nameservers delegated to Cloudflare).
- **`minpl.us`** — planned 301 redirect to the canonical `minplus.xyz`,
  either as its own Cloudflare zone with a redirect rule or as registrar
  URL forwarding, so SEO and links consolidate on one host.

## License

MIT (see LICENSE). Chapter content: the tutorial text of *Learn You a
Tropical Geometry for Great Good!*, converted from the tropical geometry
study (Days 0–20); every chapter's code is asserted in-line.
