#!/usr/bin/env python3
"""Write the ``info.metrics`` block into every metric component's config.

OpenProblems requires each metric component to declare, per metric it reports:
``label, summary, description, references, links, min, max, maximize``
(``common/schemas/task_metric.yaml``). Maintaining that by hand across ~70
metric ids would guarantee drift from the implementations, so it is generated
from the registry in ``sprotmetrics``, which is where direction and bounds are
already recorded.

Two deliberate behaviours:

* Metrics with no direction are omitted. The schema has no way to express
  "descriptive", and ``maximize`` is required, so including them would hand the
  scoring stage a ranking criterion the metric cannot support. They are still
  computed and still appear in the score file - they are just not OpenProblems
  metrics.
* ``references`` and ``links`` are backfilled here. ``MetricInfo`` carries both
  fields, but ``to_openproblems_info()`` in sprotmetrics 0.1.0 does not emit
  them, and the schema requires both keys to be present. Where a metric records
  no reference of its own, the task's own links are used, which is honest: the
  metric originates with this benchmark rather than with a paper.

Run from the repo root:

    python scripts/create_component/generate_metric_info.py [--check]

``--check`` exits non-zero if any config is out of date, for CI.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import sprotmetrics as sm
import yaml

REPO = "https://github.com/openproblems-bio/task_spatial_proteomics_preprocessing"
TASK_LINKS = {"repository": REPO, "docs": f"{REPO}#readme"}

# Component directory -> the sprotmetrics family it reports.
FAMILIES = {
    "basic": "basic",
    "marker_distribution": "marker_distribution",
    "normalization_diagnostics": "normalization_diagnostics",
    "segmentation_quality": "segmentation_quality",
    "technical_covariates": "technical_covariates",
    "signal_background": "signal_background",
}


def entries_for(family: str) -> list[dict]:
    out = []
    for metric_id in sm.list_metrics(family=family, include_qc_variants=True).index:
        info = sm.get_metric(metric_id)
        if info.maximize is None:
            continue  # descriptive or ambiguous; see the module docstring
        entry = dict(info.to_openproblems())
        entry["references"] = (
            {"doi": list(info.references)} if info.references else {}
        )
        entry["links"] = dict(info.links) if info.links else dict(TASK_LINKS)
        # Schema order, so a regenerated file diffs cleanly against the last one.
        out.append({
            k: entry[k] for k in
            ["name", "label", "summary", "description", "references", "links",
             "min", "max", "maximize"] if k in entry
        })
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if any config would change")
    args = ap.parse_args()

    stale = []
    for component, family in FAMILIES.items():
        path = Path("src/metrics_proteomics") / component / "config.vsh.yaml"
        if not path.exists():
            print(f"!! missing {path}", file=sys.stderr)
            return 1

        config = yaml.safe_load(path.read_text())
        wanted = entries_for(family)
        current = (config.get("info") or {}).get("metrics")

        if current == wanted:
            print(f"   {component}: up to date ({len(wanted)} metrics)")
            continue
        stale.append(component)
        if args.check:
            print(f"!! {component}: out of date", file=sys.stderr)
            continue

        config.setdefault("info", {})["metrics"] = wanted
        path.write_text(yaml.safe_dump(config, sort_keys=False, width=100,
                                       default_flow_style=False, allow_unicode=True))
        print(f"   {component}: wrote {len(wanted)} metrics")

    if args.check and stale:
        print(f"\n{len(stale)} config(s) out of date; run this script without "
              f"--check", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
