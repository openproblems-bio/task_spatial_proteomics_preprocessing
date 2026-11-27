#!/usr/bin/env python3
"""Render docs/metrics_catalogue.csv as a colour-coded Excel workbook.

The CSV is the source of truth; this only styles it. Run from the repo root:

    python docs/make_metrics_xlsx.py [OUTPUT.xlsx]
"""

from __future__ import annotations

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

CSV_PATH = Path("docs/metrics_catalogue.csv")

# Status drives the row's identity colour: what exists, what is new, what is
# still open, what was rejected.
STATUS_STYLE = {
    "IMPLEMENTED NOW":                            ("C8E6C9", "1B5E20", True),
    "pre-existing":                               ("DCE9F7", "13395E", False),
    "pre-existing (structurally dead)":           ("E2E2E2", "5A5A5A", False),
    "pre-existing (not wired into the workflow)": ("FFE7BA", "7A4E00", False),
    "CANDIDATE (not implemented)":                ("E6DAF2", "4A2C6B", False),
    "BLOCKED":                                    ("FFD5C2", "8A3A12", False),
    "NOT ADAPTING":                               ("F7D7DA", "8A1C24", False),
}

# Order matters twice over. The explicit non-directions come first, because
# their parentheticals can mention "high" or "low" in passing - "none
# (descriptive: high suits a lineage marker...)" is not a high-is-better metric.
# Then "weakly high" before "high", or it takes the firmer colour.
DIRECTION_STYLE = [
    ("none", "EDEDED", "6B6B6B"),
    ("n/a", "EDEDED", "6B6B6B"),
    ("constant", "EDEDED", "6B6B6B"),
    ("ambiguous", "FFF0D0", "7A5A00"),
    ("weakly high", "EAF6EB", "3E7B44"),
    ("weakly low", "EAF1FB", "3A5F80"),
    ("high", "D6EFD8", "1B5E20"),
    ("low", "D9E7FA", "13395E"),
]

# Values in the "separates" column that mean "no signal / not measured".
NO_SIGNAL = {"none", "n/a", "unknown", "not run", "cannot run", ""}

HEADER_FILL = PatternFill("solid", fgColor="21354B")
HEADER_FONT = Font(bold=True, color="FFFFFF", size=11)
THIN = Side(style="thin", color="C9CDD2")
BORDER = Border(bottom=THIN, right=THIN)

WIDTHS = {
    "metric_id": 46, "component": 24, "task_axis": 26, "what_it_measures": 58,
    "range": 14, "desired_direction": 30, "status": 34, "origin": 46,
    "inputs_needed": 18, "prior_knowledge_needed": 30,
    "separates_methods_on_exemplar_001": 24, "notes": 88,
}
WRAP = {"what_it_measures", "notes", "origin", "desired_direction",
        "prior_knowledge_needed", "status", "task_axis"}


def direction_style(value: str):
    """Fill and font for a direction, matching the most specific label first."""
    v = value.lower().strip()
    for key, fill, font in DIRECTION_STYLE:
        if key in v:
            return fill, font
    return None, None


def main() -> None:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("metrics_catalogue.xlsx")
    rows = list(csv.DictReader(CSV_PATH.open()))
    cols = list(rows[0].keys())

    wb = Workbook()
    ws = wb.active
    ws.title = "Metrics"

    ws.append([c.replace("_", " ") for c in cols])
    for cell in ws[1]:
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(vertical="center", wrap_text=True)
    ws.row_dimensions[1].height = 34

    for r in rows:
        ws.append([r[c] for c in cols])
        i = ws.max_row
        fill, font_colour, bold = STATUS_STYLE.get(r["status"], ("FFFFFF", "000000", False))

        for j, col in enumerate(cols, start=1):
            cell = ws.cell(row=i, column=j)
            cell.border = BORDER
            cell.alignment = Alignment(
                vertical="top", wrap_text=col in WRAP,
                horizontal="center" if col in {"range", "inputs_needed"} else "left",
            )
            # Row-wide tint from the status, so a glance groups the table.
            cell.fill = PatternFill("solid", fgColor=fill)
            cell.font = Font(color=font_colour, size=10,
                             bold=bold and col == "metric_id")

        # Columns that carry their own signal get a stronger, local colour.
        ws.cell(row=i, column=cols.index("metric_id") + 1).font = Font(
            color=font_colour, size=10, bold=True
        )
        ws.cell(row=i, column=cols.index("status") + 1).font = Font(
            color=font_colour, size=10, bold=True
        )

        d_cell = ws.cell(row=i, column=cols.index("desired_direction") + 1)
        d_fill, d_font = direction_style(r["desired_direction"])
        if d_fill:
            d_cell.fill = PatternFill("solid", fgColor=d_fill)
            muted = d_font in {"6B6B6B", "7A5A00"}
            d_cell.font = Font(color=d_font, size=10, bold=not muted, italic=muted)

        s_cell = ws.cell(row=i, column=cols.index("separates_methods_on_exemplar_001") + 1)
        sep = r["separates_methods_on_exemplar_001"].strip().lower()
        if sep in NO_SIGNAL:
            s_cell.fill = PatternFill("solid", fgColor="EDEDED")
            s_cell.font = Font(color="6B6B6B", size=10, italic=True)
        else:
            s_cell.fill = PatternFill("solid", fgColor="D6EFD8")
            s_cell.font = Font(color="1B5E20", size=10, bold=True)

    for j, col in enumerate(cols, start=1):
        ws.column_dimensions[get_column_letter(j)].width = WIDTHS.get(col, 20)
    ws.freeze_panes = "B2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(cols))}{ws.max_row}"

    _summary_sheet(wb, rows)
    wb.save(out_path)
    print(f"wrote {out_path}  ({len(rows)} metrics)")


def _summary_sheet(wb: Workbook, rows: list[dict]) -> None:
    ws = wb.create_sheet("Summary")

    def header(text, row, span=3):
        c = ws.cell(row=row, column=1, value=text)
        c.font = Font(bold=True, color="FFFFFF", size=11)
        c.fill = HEADER_FILL
        for j in range(2, span + 1):
            ws.cell(row=row, column=j).fill = HEADER_FILL

    ws.column_dimensions["A"].width = 46
    ws.column_dimensions["B"].width = 12
    ws.column_dimensions["C"].width = 74

    header("Status", 1)
    ws.cell(row=1, column=2, value="Count").font = Font(bold=True, color="FFFFFF")
    ws.cell(row=1, column=3, value="Meaning").font = Font(bold=True, color="FFFFFF")

    meaning = {
        "IMPLEMENTED NOW": "Added in this change; adapted from the sprot prototype or written from the catalogue",
        "pre-existing": "Already in the benchmark before this work",
        "pre-existing (structurally dead)": "Cannot vary given the methods implemented, so carries no signal",
        "pre-existing (not wired into the workflow)": "Component exists but the workflow never calls it",
        "CANDIDATE (not implemented)": "Worth adding; not done yet",
        "BLOCKED": "Cannot be implemented until missing input is available",
        "NOT ADAPTING": "Deliberately rejected; see notes for the reason",
    }
    r = 2
    for status, n in sorted(Counter(x["status"] for x in rows).items(), key=lambda kv: -kv[1]):
        fill, font, _ = STATUS_STYLE.get(status, ("FFFFFF", "000000", False))
        for j, v in enumerate([status, n, meaning.get(status, "")], start=1):
            c = ws.cell(row=r, column=j, value=v)
            c.fill = PatternFill("solid", fgColor=fill)
            c.font = Font(color=font, size=10, bold=j == 1)
            c.alignment = Alignment(vertical="top", wrap_text=j == 3)
            c.border = BORDER
        r += 1

    r += 1
    header("Directional metrics per pipeline axis", r)
    ws.cell(row=r, column=2, value="Before").font = Font(bold=True, color="FFFFFF")
    ws.cell(row=r, column=3, value="Now").font = Font(bold=True, color="FFFFFF")
    r += 1

    # A metric only helps rank a method if it has a direction and is live.
    def directional(x):
        d = x["desired_direction"].lower()
        return d.startswith(("high", "low")) and "dead" not in x["status"]

    before, now = defaultdict(int), defaultdict(int)
    for x in rows:
        if not directional(x):
            continue
        for axis in [a.strip() for a in x["task_axis"].split(",")]:
            if x["status"].startswith("pre-existing"):
                before[axis] += 1
                now[axis] += 1
            elif x["status"] == "IMPLEMENTED NOW":
                now[axis] += 1

    for axis in ["segmentation", "image correction", "aggregation", "normalization", "QC"]:
        b, n = before.get(axis, 0), now.get(axis, 0)
        ws.cell(row=r, column=1, value=axis).font = Font(size=10, bold=True)
        for j, v in ((2, b), (3, n)):
            c = ws.cell(row=r, column=j, value=v)
            c.alignment = Alignment(horizontal="center")
            # Red where an axis had nothing to rank on, green where it now does.
            c.fill = PatternFill("solid", fgColor="F7D7DA" if v == 0 else "C8E6C9")
            c.font = Font(color="8A1C24" if v == 0 else "1B5E20", size=10, bold=True)
        for j in range(1, 4):
            ws.cell(row=r, column=j).border = BORDER
        r += 1

    r += 1
    note = ws.cell(
        row=r, column=1,
        value=("A metric counts as directional only if it has a known good direction "
               "(high or low) and can actually vary. Descriptive metrics such as "
               "median_area, and metrics that are structurally constant such as the "
               "rank-correlation pair, are excluded."),
    )
    note.font = Font(italic=True, size=9, color="5A5A5A")
    note.alignment = Alignment(wrap_text=True, vertical="top")
    ws.merge_cells(start_row=r, start_column=1, end_row=r + 2, end_column=3)


if __name__ == "__main__":
    main()
