#!/usr/bin/env python3
"""compose_main_figures_vector.py — vector (Illustrator-editable) composites of Figures 1-4, with the same geometry as
compose_main_figures.R (which rasterises at 300 dpi). Panels are the compact vector PDFs of the figure driver / ladder script /
A19, and the schematic and browser panels are clipped regions of the submitted figure PDFs (clipped through the page media box).
Document tooling only (2026-11-07, Larry: 'Make the pdf too so I can modify in the illustrator').
usage: python3 revision_modules/compose_main_figures_vector.py DRIVER_OUT_DIR      (run from Weak_uORF_submission/)"""
import sys, os
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import RectangleObject
DRV = sys.argv[1] if len(sys.argv) > 1 else '_to_delete/driver_compact_20261028'
OUT = 'output_revised_figs_20261028'; WK = f'{OUT}/composite_work'; VEC = f'{OUT}/vector'; os.makedirs(VEC, exist_ok=True)
SUB = '/Users/wu/Desktop/weak_uORF_revision/Submit_Genome Research/Figures 20260826'
PT = 72.0
import warnings; warnings.filterwarnings('ignore', category=FutureWarning)
import fitz   # PyMuPDF, document tooling (V23, 2026-11-08)
_CLEAN = {}
def clean_crop(src_pdf, crop):
    """V23 (2026-11-08): a copy of the source page that keeps only what lies inside the crop (x, y, w, h in points from the top-left).
    Text, vector graphics and image pixels outside the crop are removed by redaction, so the composite carries no hidden content under
    the clip (before, the whole submitted page sat under every crop: old labels stayed searchable and Figure 1 was 42 MB)."""
    key = (src_pdf, tuple(crop))
    if key in _CLEAN: return _CLEAN[key]
    doc = fitz.open(src_pdf); pg = doc[0]; W, H = pg.rect.width, pg.rect.height
    cx, cy, cw, ch = crop; r = fitz.Rect(cx, cy, cx + cw, cy + ch)
    for rr in (fitz.Rect(0, 0, W, r.y0), fitz.Rect(0, r.y1, W, H), fitz.Rect(0, r.y0, r.x0, r.y1), fitz.Rect(r.x1, r.y0, W, r.y1)):
        if rr.width > 0.5 and rr.height > 0.5: pg.add_redact_annot(rr, fill=False)
    pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_PIXELS, graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED, text=fitz.PDF_REDACT_TEXT_REMOVE)
    out = f'{WK}/_clean_{len(_CLEAN)}_{os.path.basename(src_pdf)}'; doc.save(out, garbage=4, deflate=True); doc.close()
    _CLEAN[key] = out; return out

class Page:
    def __init__(self, file, H):
        self.w = PdfWriter(); self.H = H; self.p = self.w.add_blank_page(width=6.5 * PT, height=H * PT); self.file = file
    def place(self, src_pdf, x, y, w, crop=None):
        """Place a PDF page (or a crop of it: x, y, w, h in points from the page's top-left) at (x, y) inches from the top-left, w inches wide; returns the placed height in inches."""
        if crop: src_pdf = clean_crop(src_pdf, crop)   # V23: no hidden content outside the crop
        pg = PdfReader(src_pdf).pages[0]; pw, ph = float(pg.mediabox.width), float(pg.mediabox.height)
        if crop:
            cx, cy, cw, ch = crop; x0, y0 = cx, ph - cy - ch
            pg.mediabox = RectangleObject([x0, y0, x0 + cw, y0 + ch]); pg.cropbox = pg.mediabox
        else:
            x0, y0, cw, ch = float(pg.mediabox.left), float(pg.mediabox.bottom), pw, ph
        s = w * PT / cw; h_in = w * ch / cw
        t = Transformation().translate(-x0, -y0).scale(s, s).translate(x * PT, (self.H - y) * PT - h_in * PT)
        self.p.merge_transformed_page(pg, t, over=True)
        return h_in
    def white(self, x, y, w, h):
        """A white rectangle (x, y, w, h in inches from the top-left), drawn as a tiny white page scaled to size."""
        from pypdf import PdfWriter as W2
        tmp = W2(); bp = tmp.add_blank_page(width=10, height=10)
        from pypdf.generic import ContentStream, NameObject, DecodedStreamObject
        cs = DecodedStreamObject(); cs.set_data(b'1 1 1 rg 0 0 10 10 re f'); bp[NameObject('/Contents')] = tmp._add_object(cs)
        tmpfile = f'{WK}/_white_tmp.pdf'; tmp.write(tmpfile)
        pg = PdfReader(tmpfile).pages[0]
        t = Transformation().scale(w * PT / 10, h * PT / 10).translate(x * PT, (self.H - y - h) * PT)
        self.p.merge_transformed_page(pg, t, over=True)
    def save(self):
        self.w.write(self.file); print('saved', self.file)

## Figure 1: A schematic (submitted, relabelled), B uIF distribution (driver), C-F browser rows (submitted rows D-G)
f1 = f'{WK}/Fig1_submitted_relabelled_v16.pdf'
P = Page(f'{VEC}/Figure1_vector_composite.pdf', 8.3)
P.place(f1, 0.00, 0.05, 1.96, crop=(0, 0, 141, 134))
P.place(f'{DRV}/Figure1_compact_B_dist.pdf', 2.00, 0.00, 3.40)
y = 2.05; s = 5.9 / 468
for crop in [(0, 134, 468, 129), (0, 273, 468, 119), (0, 404, 468, 117), (0, 532, 468, 118)]:
    h = P.place(f1, 0.3, y, 5.9, crop=crop)
    if crop[1] == 134: P.white(0.3 + 100 * s, y, 320 * s, 3 * s)   # blank the descenders of the submitted B/C axis titles at the top of row C (as the raster composite does)
    y += h + 0.02
P.save(); print('  Figure 1 bottom at %.2f in' % y)

## Figure 2: A schematic (submitted), B-F compact panels, G table (A19)
f2 = f'{SUB}/Fig2_uORF+_v20260414.pdf'
H_BF, A_POS = 4.4, (0.08, 0.05, 1.85)   # V30 (2026-11-09): as in compose_main_figures.R (B-F grid height; x, y, width of the schematic A)
# V54 (2026-11-23, Larry: "Reletter"): the table (A19) is panel D between the two rows of the grid (cut at 2.2 in, in the white band between the rows); E-G below
SPLIT = 2.2
P = Page(f'{VEC}/Figure2_vector_composite.pdf', H_BF + 0.05 + 2.52 + 0.05)
P.place(f'{DRV}/Figure2_compact_panels_B-F.pdf', 0, 0, 6.5, crop=(0, 0, 6.5 * PT, SPLIT * PT))
P.place(f2, A_POS[0], A_POS[1], A_POS[2], crop=(0, 0, 155, 80))
P.place('output/A19_Figure2G_TSS_vs_splicing_table.pdf', 0.0, SPLIT + 0.05, 6.5)
P.place(f'{DRV}/Figure2_compact_panels_B-F.pdf', 0, SPLIT + 0.05 + 2.52 + 0.05, 6.5, crop=(0, SPLIT * PT, 6.5 * PT, (H_BF - SPLIT) * PT))
P.save()

## Figure 3: A-B (driver) above C (ladder plot + TPM table)
P = Page(f'{VEC}/Figure3_vector_composite.pdf', 8.35)
h = P.place(f'{DRV}/Figure3_compact_panels_A-B.pdf', 0, 0, 6.5); h2 = P.place(f'{OUT}/Figure3_revised_panel_C_prevalence_ladder_plot.pdf', 0.0, h + 0.1, 6.5)
P.place(f'{OUT}/Figure3_revised_panel_C_prevalence_ladder_table.pdf', 0, h + 0.1 + h2 + 0.05, 6.5); P.save()

## Figure 4: A-C browser panels (submitted, scaled to 6.5 in), D-F compact
f4 = f'{SUB}/Fig4_CpuORF_v20260511pdf.pdf'
P = Page(f'{VEC}/Figure4_vector_composite.pdf', 7.6)   # V30 (2026-11-09): the D-F row is 2.1 in tall
y = 0
for crop in [(0, 0, 403, 105), (0, 108, 403, 96), (0, 207, 403, 121)]:
    h = P.place(f4, 0, y, 6.5, crop=crop); y += h + 0.05
P.place(f'{DRV}/Figure4_compact_panels_D-F.pdf', 0, y, 6.5)   # V30 (2026-11-09): one aligned row, legends inside the plots
P.save(); print('  Figure 4 D-F row at %.2f in' % y)
