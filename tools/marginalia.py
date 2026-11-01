#!/usr/bin/env python3
"""Marginalia post-processing for minplus.xyz (run after zola build +
serve-the-source, as the last step of `npm run build`).

For every chapter page in public/chapters/chapter-NN/:

1. ANCHORS. Number the chapter's top-level paragraphs within each section
   and stamp each with an id (s1-3-p2) plus a hidden human label span
   (§1.3 ¶2). The labels become visible in the print edition, so a note
   scribbled in a paper margin can be addressed exactly.

2. MARGIN NOTES. Read content/chapters/chapter-NN/annotations.json
   ({"notes": [{"anchor": "s1-3-p2", "text": "...", "ink": "ink/x.png"}]}).
   Each note is injected after its paragraph as an <aside class="margin-note">
   — real HTML text: selectable, copyable, searchable. If a chapter has any
   notes, a CSS-only toggle (checkbox hack, no JS) is added to show/hide them.

3. PRINT EDITION. Write public/chapters/chapter-NN/print/index.html: the same
   page with a print stylesheet (Letter, wide right margin for handwriting,
   anchor labels visible in the left margin, chrome hidden).

The .md mirrors from serve-the-source are produced before this step and are
not touched: the raw Markdown stays the clean chapter text.
"""
import json
import re
import sys
from copy import deepcopy
from pathlib import Path

from lxml import etree
from lxml import html as LH

ROOT = Path(__file__).resolve().parent.parent
PUBLIC = ROOT / "public"
CONTENT = ROOT / "content"

PRINT_CSS = """
@page { size: Letter; margin: 20mm 58mm 20mm 24mm; }
body { font-size: 11pt; }
.site-header, .site-footer, .post-footer, .skip-link, .margin-toggle,
.margin-toggle-label { display: none !important; }
.container { max-width: none; padding: 0; }
.prose { max-width: none; }
.para-anchor { display: block !important; float: left; width: 20mm;
  margin-left: -22mm; text-align: right; padding-right: 2mm;
  font-family: "JetBrains Mono", monospace; font-size: 7.5pt; color: #999; }
.margin-note { display: block !important; border-left: 2px solid #6d4aff;
  padding-left: 3mm; font-style: italic; }
.print-note { font-family: "JetBrains Mono", monospace; font-size: 8.5pt;
  color: #666; border: 1px dashed #bbb; padding: 3mm; margin-bottom: 8mm; }
img { max-width: 100%; }
"""


def section_of(h2_text):
    """Return (id_fragment, display_label) for a section heading."""
    text = " ".join(h2_text.split())
    m = re.match(r"(\d+)\.(\d+)\b", text)
    if m:
        return f"s{m.group(1)}-{m.group(2)}", f"§{m.group(1)}.{m.group(2)}"
    slug = re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")[:24] or "front"
    return f"s-{slug}", text.split("—")[0].strip()[:28]


def process(chapter_dir: Path, content_dir: Path):
    page_path = chapter_dir / "index.html"
    # Force UTF-8 both ways: lxml's encoding sniffing misreads Zola's
    # output as Latin-1 and double-encodes every non-ASCII character
    # (mojibake: ⊕ -> âŠ•). Learned the hard way, 2026-11-01.
    parser = LH.HTMLParser(encoding="utf-8")
    tree = LH.parse(str(page_path), parser)
    root = tree.getroot()
    prose = None
    for div in root.iter("div"):
        classes = (div.get("class") or "").split()
        if "prose" in classes:
            prose = div
            break
    if prose is None:
        print(f"  !! no .prose container in {page_path}", file=sys.stderr)
        return

    # --- anchors ------------------------------------------------------------
    anchors = {}  # pid -> element
    sec_id, sec_label = "s-front", "Intro"
    pcount = 0
    for child in list(prose):
        if child.tag == "h2":
            sec_id, sec_label = section_of(child.text_content())
            pcount = 0
        elif child.tag == "p":
            if not child.text_content().strip():
                continue
            pcount += 1
            pid = f"{sec_id}-p{pcount}"
            child.set("id", pid)
            span = LH.fromstring(
                f'<span class="para-anchor" aria-hidden="true"></span>')
            span.text = f"{sec_label} ¶{pcount}"
            child.insert(0, span)
            anchors[pid] = child

    # --- margin notes ---------------------------------------------------------
    notes = []
    ann_file = content_dir / "annotations.json"
    if ann_file.exists():
        notes = json.loads(ann_file.read_text()).get("notes", [])
    placed = 0
    for note in notes:
        target = anchors.get(note["anchor"])
        if target is None:
            print(f"  !! note anchor {note['anchor']} not found in "
                  f"{chapter_dir.name} — skipped", file=sys.stderr)
            continue
        aside = LH.Element("aside")
        aside.set("class", "margin-note")
        aside.set("data-anchor", note["anchor"])
        label = target.find("span").text if target.find("span") is not None else note["anchor"]
        a_span = etree.SubElement(aside, "span")
        a_span.set("class", "margin-note__anchor")
        a_span.text = label
        text_p = etree.SubElement(aside, "p")
        text_p.set("class", "margin-note__text")
        text_p.text = note["text"]
        if note.get("ink"):
            img = etree.SubElement(aside, "img")
            img.set("class", "margin-note__ink")
            img.set("src", note["ink"])
            img.set("alt", f"Handwritten margin note at {label}")
            img.set("loading", "lazy")
        target.addnext(aside)
        placed += 1

    if placed:
        toggle = LH.fromstring(
            '<input type="checkbox" class="margin-toggle" '
            'id="margin-toggle">')
        label_el = LH.fromstring(
            '<label class="margin-toggle-label" for="margin-toggle">'
            "✎ The author's margin notes</label>")
        prose.insert(0, label_el)
        prose.insert(0, toggle)

    tree.write(str(page_path), method="html", encoding="utf-8")

    # --- print edition --------------------------------------------------------
    print_tree = deepcopy(tree)
    proot = print_tree.getroot()
    head = proot.find("head")
    style = etree.SubElement(head, "style")
    style.text = PRINT_CSS
    pprose = None
    for div in proot.iter("div"):
        if "prose" in (div.get("class") or "").split():
            pprose = div
            break
    note = LH.fromstring(
        '<p class="print-note">Print edition for handwritten margin notes. '
        "Anchors in the left margin (§1.3 ¶2 style) match the web edition — "
        "notes digitised from this paper appear at the same anchors online. "
        "minplus.xyz</p>")
    header = pprose.find("header")
    if header is not None:
        header.addnext(note)
    else:
        pprose.insert(0, note)
    out_dir = chapter_dir / "print"
    out_dir.mkdir(exist_ok=True)
    print_tree.write(str(out_dir / "index.html"), method="html", encoding="utf-8")
    print(f"  {chapter_dir.name}: {len(anchors)} anchors, {placed} notes, "
          f"print edition written")


def main():
    chapters = sorted((PUBLIC / "chapters").glob("chapter-*"))
    if not chapters:
        print("no chapter output found — run zola build first", file=sys.stderr)
        return 1
    for ch in chapters:
        process(ch, CONTENT / "chapters" / ch.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
