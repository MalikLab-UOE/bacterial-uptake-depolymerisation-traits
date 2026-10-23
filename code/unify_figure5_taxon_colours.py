#!/usr/bin/env python3

"""Unify the taxon colours in the Illustrator-assembled Figure 5 PDF."""

import argparse
from pathlib import Path

from pypdf import PdfReader, PdfWriter
from pypdf.generic import StreamObject


REPLACEMENTS = {
    b"0.906 0.435 0.318": b"0.957 0.639 0.639",  # Bacteria, upper c
    b"1 0.431 0.412": b"0.957 0.639 0.639",      # Bacteria, lower c
    b"0.165 0.616 0.561": b"0.349 0.855 0.847",  # Fungi, upper c
    b"0 0.765 0.769": b"0.349 0.855 0.847",      # Fungi, lower c
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_pdf", type=Path)
    parser.add_argument("output_pdf", type=Path)
    args = parser.parse_args()

    reader = PdfReader(args.input_pdf)
    writer = PdfWriter(clone_from=reader)
    changed_streams = 0

    for obj in writer._objects:
        if not isinstance(obj, StreamObject):
            continue
        try:
            data = obj.get_data()
        except Exception:
            continue

        updated = data
        for old, new in REPLACEMENTS.items():
            updated = updated.replace(old, new)
        if updated != data:
            obj.set_data(updated)
            changed_streams += 1

    if changed_streams == 0:
        raise RuntimeError("No Figure 5 taxon-colour streams were found")

    args.output_pdf.parent.mkdir(parents=True, exist_ok=True)
    with args.output_pdf.open("wb") as output_handle:
        writer.write(output_handle)

    print(f"Updated {changed_streams} PDF streams")


if __name__ == "__main__":
    main()
