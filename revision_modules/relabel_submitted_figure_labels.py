#!/usr/bin/env python3
"""Replace the 'Strong uORF' / 'Weak uORF' wording inside two submitted figure PDFs by the uIF-class wording of the revised
manuscript (Larry, 2026-11-01; document tooling, no analysis). The old text objects are removed (PDF redaction restricted to
text; the drawings are untouched) and the new labels are written in Arial at the same size and baseline, centred where the
old label was centred (left-aligned where it was left-aligned). The output stays a vector PDF.
  Figure 1 (submitted Fig1_strong_vs_weak_uORFs_v20260817.pdf): panel A 'Strong uORF (High uIF)' -> 'High-uIF uORF' and
    'Weak uORF (Low uIF)' -> 'Low-uIF uORF'; the titles of D-G 'Strong/Weak uORF in <gene>' -> 'High-uIF/Low-uIF TuORF in <gene>'
    -> output_revised_figs_20261028/composite_work/Fig1_submitted_relabelled.pdf (compose_main_figures.R crops A and D-G from it).
  Figure 7 (submitted Fig7_Model_v4.pdf), Option A of Figure7_model_options_2026-11-01.md (Larry, 2026-11-02): the rows are
    defined by uORF inclusion instead of uIF, and the outcome is stated per mRNA. Left labels (9.5 pt, centred in the column):
    'High %uORF+ (>80%) / most TSSs upstream / e.g. CPuORFs', 'Low %uORF+ (<30%) / few TSSs upstream', 'Undetected uORFs /
    no TSS upstream'; right labels (10 pt, old baselines): 'mORF repressed / in most mRNAs', 'mORF repressed / in few mRNAs',
    'mORF / not repressed'. The captions above the mRNAs ('Most / Few / No mRNAs contain uORFs') and the drawing are unchanged.
    -> output_revised_figs_20261028/Figure7_revised_model.pdf (the journal file Wu_Fig7.pdf). A variant with a grey
    'observed uIF: high / low / none' line under each left label is written to composite_work/Figure7_model_optionA_uIF_note.pdf
    for comparison (not in the package). The V12 word swap ('High-uIF uORFs', 'Low-uIF uORFs', 'Repressed mORF') is kept in
    output_revised_figs_20261028/_before_V13_20261102/.
usage: python3 revision_modules/relabel_submitted_figure_labels.py        (run from Weak_uORF_submission/; needs PyMuPDF)"""
import os, sys, warnings
warnings.simplefilter('ignore')
import pymupdf
SUB = '/Users/wu/Desktop/weak_uORF_revision/Submit_Genome Research/Figures 20260826'
OUT = 'output_revised_figs_20261028'; WK = os.path.join(OUT, 'composite_work'); os.makedirs(WK, exist_ok=True)
REG = '/System/Library/Fonts/Supplemental/Arial.ttf'; ITA = '/System/Library/Fonts/Supplemental/Arial Italic.ttf'
F = {False: pymupdf.Font(fontfile=REG), True: pymupdf.Font(fontfile=ITA)}

def spans(page):
    for b in page.get_text('dict')['blocks']:
        if b['type'] == 0:
            for l in b['lines']:
                for s in l['spans']: yield s

def find(page, text):
    hits = [s for s in spans(page) if s['text'].strip() == text]
    assert len(hits) == 1, (text, [s['text'] for s in hits]); return hits[0]

def visible_extent(s):
    """x range of the glyphs without the leading / trailing spaces of the span."""
    sp = F[False].text_length(' ', fontsize=s['size']); t = s['text']
    return s['bbox'][0] + sp * (len(t) - len(t.lstrip(' '))), s['bbox'][2] - sp * (len(t) - len(t.rstrip(' ')))

def relabel(page, jobs, tag):
    """jobs: (old span texts of one label line, new parts [(text, italic)], 'left' | 'centre')."""
    all_before = [(s['text'], round(s['origin'][0], 1), round(s['origin'][1], 1), s['size'], s['font']) for s in spans(page)]
    rects, plan, expected = [], [], set()
    for olds, new, align in jobs:
        ss = [find(page, o) for o in olds]; size = ss[0]['size']; base = ss[0]['origin'][1]
        x0 = min(visible_extent(s)[0] for s in ss); x1 = max(visible_extent(s)[1] for s in ss)
        for s in ss:                                      # a rectangle inside the line box, so the lines above and below are not touched
            rects.append(pymupdf.Rect(s['bbox'][0] + 0.3, s['bbox'][1] + 1.0, s['bbox'][2] - 0.3, s['bbox'][3] - 1.5)); expected.add(s['text'])
        w = sum(F[it].text_length(t, fontsize=size) for t, it in new)
        plan.append((x0 if align == 'left' else (x0 + x1) / 2 - w / 2, base, size, new))
    for r in rects: page.add_redact_annot(r)
    page.apply_redactions(images=0, graphics=0, text=0)   # text only
    after = set((s['text'], round(s['origin'][0], 1), round(s['origin'][1], 1), s['size'], s['font']) for s in spans(page))
    removed = [s for s in all_before if s not in after]
    missing = [o for olds, _, _ in jobs for o in olds if any(s['text'].strip() == o for s in spans(page))]   # checked before the new text goes in (gene names recur)
    assert not missing, ('old text still present', missing)
    collateral = [s for s in removed if s[0] not in expected]
    for t, x, y, size, font in collateral:                # a neighbouring span caught by a rectangle: put it back unchanged
        page.insert_text(pymupdf.Point(x, y), t, fontsize=size, fontname='ArialI' if 'Italic' in font else 'ArialR', fontfile=ITA if 'Italic' in font else REG, color=0)
    for x, base, size, new in plan:
        for t, it in new:
            page.insert_text(pymupdf.Point(x, base), t, fontsize=size, fontname='ArialI' if it else 'ArialR', fontfile=ITA if it else REG, color=0)
            x += F[it].text_length(t, fontsize=size)
    print(f'{tag}: removed {len(removed)} spans ({len(expected)} expected; {len(collateral)} neighbouring spans restored: {[c[0] for c in collateral]})')

# ---- Figure 1 -------------------------------------------------------------------------------------------------------------
doc = pymupdf.open(os.path.join(SUB, 'Fig1_strong_vs_weak_uORFs_v20260817.pdf')); pg = doc[0]
relabel(pg, [(['Strong uORF (High uIF)'], [('High-uIF uORF', False)], 'left'),
             (['Weak uORF (Low uIF)'], [('Low-uIF uORF', False)], 'left'),
             (['Strong uORF in AT1G22882', 'SUN3'], [('High-uIF TuORF in AT1G22882 ', False), ('SUN3', True)], 'centre'),
             (['Strong uORF in AT2G46690', 'SAUR32'], [('High-uIF TuORF in AT2G46690 ', False), ('SAUR32', True)], 'centre'),
             (['Weak uORF in AT2G43010', 'PIF4'], [('Low-uIF TuORF in AT2G43010 ', False), ('PIF4', True)], 'centre'),
             (['Weak uORF in AT3G54220', 'SCR'], [('Low-uIF TuORF in AT3G54220 ', False), ('SCR', True)], 'centre')], 'Figure 1')
f1 = os.path.join(WK, 'Fig1_submitted_relabelled.pdf'); doc.save(f1, garbage=3, deflate=True); doc.close(); print('saved', f1)

# ---- Figure 7 (Option A: rows by uORF inclusion, outcome per mRNA) -------------------------------------------------------------
GREY = (0x58 / 255, 0x59 / 255, 0x5b / 255)   # the grey of the 'uORF  mORF' labels in the drawing
BLACK = (0, 0, 0)

def find_near(page, text, y=None):
    """the one span whose stripped text is `text` (and whose baseline is within 3 pt of y when given; 'mORF' occurs three times)."""
    hits = [s for s in spans(page) if s['text'].strip() == text and (y is None or abs(s['origin'][1] - y) < 3)]
    assert len(hits) == 1, (text, y, [(s['text'], s['origin']) for s in hits]); return hits[0]

def redraw(page, remove, blocks, tag):
    """remove: [(span text, baseline y or None)] -> those text objects are removed (text-only redaction of a rectangle inside the
    line box). blocks: [(centre x, first baseline, [(text, size, leading before this line, colour)])] -> each line centred at
    centre x, written in Arial."""
    all_before = [(s['text'], round(s['origin'][0], 1), round(s['origin'][1], 1), s['size'], s['font']) for s in spans(page)]
    rects, expected = [], set()
    for t, y in remove:
        s = find_near(page, t, y); expected.add(s['text'])
        rects.append(pymupdf.Rect(s['bbox'][0] + 0.3, s['bbox'][1] + 1.0, s['bbox'][2] - 0.3, s['bbox'][3] - 1.5))
    for r in rects: page.add_redact_annot(r)
    page.apply_redactions(images=0, graphics=0, text=0)   # text only; the drawing is untouched
    after = set((s['text'], round(s['origin'][0], 1), round(s['origin'][1], 1), s['size'], s['font']) for s in spans(page))
    removed = [s for s in all_before if s not in after]
    missing = [t for t, y in remove if any(s['text'].strip() == t and (y is None or abs(s['origin'][1] - y) < 3) for s in spans(page))]
    assert not missing, ('old text still present', missing)
    collateral = [s for s in removed if s[0] not in expected]
    for t, x, y, size, font in collateral:                # a neighbouring span caught by a rectangle: put it back unchanged
        page.insert_text(pymupdf.Point(x, y), t, fontsize=size, fontname='ArialI' if 'Italic' in font else 'ArialR', fontfile=ITA if 'Italic' in font else REG, color=0)
    for cx, b, lines in blocks:
        for i, (t, size, lead, col) in enumerate(lines):
            if i: b += lead
            w = F[False].text_length(t, fontsize=size)
            page.insert_text(pymupdf.Point(cx - w / 2, b), t, fontsize=size, fontname='ArialR', fontfile=REG, color=col)
            assert cx - w / 2 >= 4 and cx + w / 2 <= page.rect.width - 4, ('line leaves the page', t)
    print(f'{tag}: removed {len(removed)} spans ({len(expected)} expected; {len(collateral)} neighbouring spans restored: {[c[0] for c in collateral]})')

def first_baseline(centre_y, lines):
    """first baseline of a block whose visual centre (x-height centre of the middle of the block) is at centre_y."""
    return centre_y - sum(l[2] for l in lines[1:]) / 2 + 0.35 * lines[0][1]

L, LS = 9.5, 11.5       # left labels: size and leading (the column between the page edge and the gene drawing is about 100 pt wide)
R, RS = 10.0, 12.0      # right labels: the size and leading of the submitted labels (old baselines 65.8 / 160.0 / 256.1)
LEFT_X, RIGHT_X = 56.0, 418.0                                   # column centres (gene drawing starts at x = 111.6; A(n) ends at 370.8)
ROW_Y = (59.0, 159.0, 250.0)                                    # visual centres of the three rows (TSS label to the uORF/mORF labels)
REMOVE = [('Strong uORFs', None), ('or CPuORFs', None), ('Weak uORFs', None), ('Undetected uORFs', None),
          ('Strongly repressed', None), ('mORF', 77.8), ('Less repressed', None), ('mORF', 172.0), ('Non-repressed', None), ('mORF', 268.1)]
LEFT = [[('High %uORF+ (>80%)', L, 0, BLACK), ('most TSSs upstream', L, LS, BLACK), ('e.g. CPuORFs', L, LS, BLACK)],
        [('Low %uORF+ (<30%)', L, 0, BLACK), ('few TSSs upstream', L, LS, BLACK)],
        [('Undetected uORFs', L, 0, BLACK), ('no TSS upstream', L, LS, BLACK)]]
UIF_NOTE = [('observed uIF: high', 8.5, 10.5, GREY), ('observed uIF: low', 8.5, 10.5, GREY), ('observed uIF: none', 8.5, 10.5, GREY)]
RIGHT = [(RIGHT_X, 65.8, [('mORF repressed', R, 0, BLACK), ('in most mRNAs', R, RS, BLACK)]),
         (RIGHT_X, 160.0, [('mORF repressed', R, 0, BLACK), ('in few mRNAs', R, RS, BLACK)]),
         (RIGHT_X, 256.1, [('mORF', R, 0, BLACK), ('not repressed', R, RS, BLACK)])]
f7 = os.path.join(OUT, 'Figure7_revised_model.pdf'); f7b = os.path.join(WK, 'Figure7_model_optionA_uIF_note.pdf')
for out, note in ((f7, False), (f7b, True)):
    doc = pymupdf.open(os.path.join(SUB, 'Fig7_Model_v4.pdf')); pg = doc[0]
    left = [lines + ([UIF_NOTE[i]] if note else []) for i, lines in enumerate(LEFT)]
    blocks = [(LEFT_X, first_baseline(ROW_Y[i], lines), lines) for i, lines in enumerate(left)] + RIGHT
    redraw(pg, REMOVE, blocks, 'Figure 7' + (' (uIF-note variant)' if note else ''))
    doc.save(out, garbage=3, deflate=True); doc.close(); print('saved', out)
for f in (f1, f7, f7b):
    d = pymupdf.open(f); words = d[0].get_text('words'); d.close()
    left = [w for w in words if 'Strong' in w[4] or 'Weak' in w[4]]
    if f == f1:   # the key of the submitted panels B/C (x > 370 pt, y < 110 pt) keeps the old wording; those panels are regenerated and never cropped
        assert all(w[0] > 370 and w[3] < 110 for w in left), (f, 'old wording inside a cropped region', left)
        print(f'{f}: old wording left only in the unused panels B/C key: {[w[4] for w in left]}')
    else:
        assert not left, (f, 'strong/weak wording left', left)
        assert (f == f7b) == any('uIF' in w[4] for w in words), (f, 'uIF wording')
        print(f'{f}: no Strong/Weak wording; labels: {[w[4] for w in words if w[0] < 110 or w[0] > 370]}')
