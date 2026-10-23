#!/usr/bin/env python3

"""Assemble Figure 5 with the family-level panel in a two-row layout."""

from io import BytesIO
from pathlib import Path

from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import RectangleObject
from reportlab.pdfgen import canvas


BASE_DIR = Path(__file__).resolve().parent.parent
MAIN_SOURCE = BASE_DIR / "results/figures/Figure_5_unified_colours.pdf"
FAMILY_SOURCE = BASE_DIR / "results/figures/Figure_5_family_panel.pdf"
OUTPUT = BASE_DIR / "figures/Figure_5.pdf"

PAGE_WIDTH = 930
PAGE_HEIGHT = 900
MARGIN = 24
GAP_X = 18
GAP_Y = 16
TOP_HEIGHT = 540
BOTTOM_HEIGHT = PAGE_HEIGHT - 2 * MARGIN - GAP_Y - TOP_HEIGHT
CELL_WIDTH = (PAGE_WIDTH - 2 * MARGIN - GAP_X) / 2


def fresh_page(path: Path):
    return PdfReader(str(path)).pages[0]


def merge_crop(destination, source_path, crop, box):
    page = fresh_page(source_path)
    left, bottom, right, top = crop
    box_x, box_y, box_w, box_h = box
    page.mediabox = RectangleObject(crop)
    page.cropbox = RectangleObject(crop)

    scale = min(box_w / (right - left), box_h / (top - bottom))
    placed_w = (right - left) * scale
    placed_h = (top - bottom) * scale
    tx = box_x + (box_w - placed_w) / 2 - left * scale
    ty = box_y + (box_h - placed_h) / 2 - bottom * scale
    transform = Transformation().scale(scale).translate(tx, ty)
    destination.merge_transformed_page(page, transform, expand=False)


def make_label_overlay():
    stream = BytesIO()
    layer = canvas.Canvas(stream, pagesize=(PAGE_WIDTH, PAGE_HEIGHT))

    label_positions = {
        "a": (MARGIN + 2, PAGE_HEIGHT - MARGIN - 18),
        "b": (MARGIN + CELL_WIDTH + GAP_X + 2, PAGE_HEIGHT - MARGIN - 18),
        "c": (MARGIN + 2, MARGIN + BOTTOM_HEIGHT - 18),
        "d": (MARGIN + CELL_WIDTH + GAP_X + 2, MARGIN + BOTTOM_HEIGHT - 18),
    }

    # Cover panel letters retained in the source PDFs before adding a-d.
    layer.setFillColorRGB(1, 1, 1)
    layer.rect(68, 847, 35, 32, fill=1, stroke=0)
    layer.rect(580, 848, 40, 30, fill=1, stroke=0)
    layer.rect(483, 320, 34, 30, fill=1, stroke=0)
    layer.rect(483, 28, 34, 30, fill=1, stroke=0)

    layer.setFillColorRGB(0, 0, 0)
    layer.setFont("Helvetica-Bold", 18)
    for label, (x, y) in label_positions.items():
        layer.drawString(x, y, label)

    layer.save()
    stream.seek(0)
    return PdfReader(stream).pages[0]


def main():
    for path in (MAIN_SOURCE, FAMILY_SOURCE):
        if not path.exists():
            raise FileNotFoundError(path)

    writer = PdfWriter()
    page = writer.add_blank_page(width=PAGE_WIDTH, height=PAGE_HEIGHT)

    bottom_y = MARGIN
    top_y = MARGIN + BOTTOM_HEIGHT + GAP_Y
    left_x = MARGIN
    right_x = MARGIN + CELL_WIDTH + GAP_X

    # Original Figure 5: a occupies the left half; b and c occupy the right half.
    merge_crop(
        page, MAIN_SOURCE,
        crop=(0, 0, 415, 693),
        box=(left_x, top_y, CELL_WIDTH, TOP_HEIGHT),
    )
    # Family-level panel used as panel b of Figure 5.
    merge_crop(
        page, FAMILY_SOURCE,
        crop=(0, 0, 374, 515),
        box=(right_x, top_y, CELL_WIDTH, TOP_HEIGHT),
    )
    merge_crop(
        page, MAIN_SOURCE,
        crop=(455, 346, 908.348, 692.931),
        box=(left_x, bottom_y, CELL_WIDTH, BOTTOM_HEIGHT),
    )
    merge_crop(
        page, MAIN_SOURCE,
        crop=(455, 0, 908.348, 346),
        box=(right_x, bottom_y, CELL_WIDTH, BOTTOM_HEIGHT),
    )

    page.merge_page(make_label_overlay())
    writer.add_metadata({"/Title": "Figure_5", "/Creator": "Figure 5 assembly script"})
    with OUTPUT.open("wb") as output_handle:
        writer.write(output_handle)

    print(f"Saved: {OUTPUT}")


if __name__ == "__main__":
    main()
