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

| Setting          | Value           |
| ---------------- | --------------- |
| Framework preset | None            |
| Build command    | `npm run build` |
| Build output dir | `public`        |
| Root directory   | `/`             |

**The build command must be `npm run build`, not `zola build`.** Plain
`zola build` deploys a site that looks fine, but it writes no `index.md`
files, so `functions/_middleware.js` has nothing to serve and agents asking
for `text/markdown` silently get HTML. Check after a deploy:
`curl -sI -H 'Accept: text/markdown' https://minplus.xyz/chapters/chapter-01/`
should say `content-type: text/markdown`.

**Pin the Zola version** — add an environment variable (Settings →
Environment variables, Production **and** Preview):

```
ZOLA_VERSION = 0.22.1
```

## Domains

- **`minplus.xyz`** (apex) — the canonical domain.
- **`minpl.us`** — redirect to the canonical `minplus.xyz` (301), set up
  as a custom domain / redirect rule in Cloudflare, so SEO and links
  consolidate on one host.

## License

MIT (see LICENSE). Chapter content: the tutorial text of *Learn You a
Tropical Geometry for Great Good!*, converted from the tropical geometry
study (Days 0–20); every chapter's code is asserted in-line.
