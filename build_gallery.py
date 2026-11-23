"""Render historische Plots als beschriftete Galerietafeln ohne native Messraster."""

import argparse
import io
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
import pymupdf

COLLECTIONS = [
    ("praktikum", "Gasjet_PostProccesing/plots", "Praktikumsplots · Python"),
    ("abgabe", "Aufgabe6_Abgabe/abbildungen", "Aufgabe 6 · Abgabeabbildungen"),
    ("lang", "Aufgabe6_Abgabe/abbildungen_lang", "Aufgabe 6 · Langfassung"),
    ("minimal", "Aufgabe6_Abgabe/plots_minimal", "Aufgabe 6 · Vergleichsfassungen"),
    ("archiv", "Abbildungen", "Historische Plotfassungen · Archiv"),
]
PANEL_LIMIT = 6


def render(path):
    if path.suffix.lower() == ".pdf":
        with pymupdf.open(path) as document:
            if len(document) != 1:
                raise ValueError(f"Einseitige Plot-PDF erwartet: {path.name}")
            page = document[0]
            # Farbskalen sind schmal und bleiben lesbar. Breite eingebettete
            # Bilder werden vor dem Rendern auf grobe Farbzellen reduziert.
            for xref in {info[0] for info in page.get_images(full=True) if min(info[2:4]) >= 48}:
                image = document.extract_image(xref)
                with Image.open(io.BytesIO(image["image"])) as source:
                    reduced = source.convert("RGB")
                    reduced.thumbnail((30, 22), Image.Resampling.BOX)
                    reduced = reduced.resize(source.size, Image.Resampling.NEAREST)
                    stream = io.BytesIO()
                    reduced.save(stream, format="PNG")
                page.replace_image(xref, stream=stream.getvalue())
            scale = 700 / page.rect.width
            pix = page.get_pixmap(matrix=pymupdf.Matrix(scale, scale), alpha=False)
            return Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    with Image.open(path) as source:
        # Bei bereits flachen Archivbildern gibt es keinen getrennten Plotlayer.
        preview = source.convert("RGB")
        preview.thumbnail((85, 62), Image.Resampling.BOX)
        return preview.resize((preview.width * 6, preview.height * 6), Image.Resampling.BILINEAR)


def panel(items, title):
    size = 625
    height = 470
    sheet = Image.new("RGB", (2 * size, 52 + ((len(items) + 1) // 2) * height), "#eef4f7")
    draw = ImageDraw.Draw(sheet)
    draw.text((16, 12), title, fill="#173043", font=ImageFont.load_default(size=22))
    for i, path in enumerate(items):
        x, y = (i % 2) * size, 52 + (i // 2) * height
        draw.rounded_rectangle((x+6, y+6, x+size-6, y+height-7), radius=6, fill="white")
        figure = render(path)
        figure.thumbnail((size-30, height-58), Image.Resampling.LANCZOS)
        sheet.paste(figure, (x + (size - figure.width)//2, y+17))
        draw.text((x+16, y+height-37), path.stem, fill="#173043", font=ImageFont.load_default(size=15))
    return sheet


def build(root, target):
    target.mkdir(parents=True, exist_ok=True)
    total = 0
    for slug, location, title in COLLECTIONS:
        files = sorted(p for p in (root / location).iterdir()
                       if p.is_file() and p.suffix.casefold() in {".pdf", ".png"})
        if not files:
            raise ValueError(f"Keine Plotdateien in {location}.")
        for start in range(0, len(files), PANEL_LIMIT):
            index = start // PANEL_LIMIT + 1
            panel(files[start:start+PANEL_LIMIT], f"{title}  |  {index}").save(
                target / f"{slug}-{index:02d}.webp", format="WEBP", quality=77)
        total += len(files)
        print(f"{slug} {len(files)} Vorlagen")
    print(f"{total} Vorlagen insgesamt")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.resolve().is_relative_to(args.source_dir.resolve()):
        parser.error("Der Ausgabeordner muss außerhalb des Quellarchivs liegen.")
    build(args.source_dir, args.output_dir)
