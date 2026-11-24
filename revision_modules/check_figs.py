#!/usr/bin/env python3
"""Checks for the rebuilt figure PDFs: page size <= 6.5 x 9 in, fonts (Arial only), smallest word height (>= 8 pt),
and a list of words that look clipped at the page edge. usage: python3 check_figs.py file.pdf [...]"""
import re, subprocess, sys, html
PT = 72.0
def run(*a): return subprocess.run(a, capture_output=True, text=True).stdout
for f in sys.argv[1:]:
    info = run('pdfinfo', f); m = re.search(r'Page size:\s+([\d.]+) x ([\d.]+) pts', info)
    w, h = (float(m.group(1)) / PT, float(m.group(2)) / PT) if m else (None, None)
    pages = re.search(r'Pages:\s+(\d+)', info); pages = int(pages.group(1)) if pages else 0
    fonts = sorted(set(l.split()[0] for l in run('pdffonts', f).splitlines()[2:] if l.strip()))
    bbox = run('pdftotext', '-bbox', f, '-')
    words = [(float(a), float(b), float(c), float(d), html.unescape(t)) for a, b, c, d, t in re.findall(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">(.*?)</word>', bbox)]
    # word height ~ font size * 1.0-1.15 for Arial; flag heights below 8 pt that are not subscript digits or single glyphs of rotated titles
    small = [(round(y2 - y1, 1), t) for x1, y1, x2, y2, t in words if (y2 - y1) < 8.0 and (x2 - x1) > 1 and not (re.fullmatch(r'\d', t) or len(t) == 1)]
    small = sorted(set(small))[:8]
    edge = [t for x1, y1, x2, y2, t in words if x1 < 1 or (w and x2 > w * PT - 1) or y1 < 1 or (h and y2 > h * PT - 1)]
    size_ok = w is not None and w <= 6.51 and h <= 9.01
    fonts_ok = all(re.search(r'Arial|Helvetica', x) for x in fonts) if fonts else False
    print(f"{f.split('/')[-1]}: {w:.2f} x {h:.2f} in ({pages} p) {'OK' if size_ok else 'TOO LARGE'}; fonts {','.join(fonts) if fonts else 'none'} {'OK' if fonts_ok else 'CHECK'}; words {len(words)}; <8pt: {small if small else 'none'}; at edge: {edge[:6] if edge else 'none'}")
