#!/usr/bin/env python3
"""Build config/samples.csv from the files present in the data directory.

Sample identity is encoded in the filenames; there is no metadata file on disk.
Group is derived from the donor prefix (see AGENTS.md):
    ^Y -> Young   ^C -> Aged   ^P -> PD
so `YC4345` is Young, not Aged.

Re-running preserves any `use` column already edited by hand.

Usage: python src/make_sample_sheet.py [--data-dir GSE193688_RAW] [--pilot]
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

H5_RE = re.compile(r"^(GSM\d+)_(.+)_filtered_feature_bc_matrix\.h5$")
GROUP_PREFIX = (("Y", "Young"), ("C", "Aged"), ("P", "PD"))

PILOT_DONORS = {"C16007", "P02284"}

FIELDS = ["sample", "gsm", "donor", "group", "h5", "fragments", "peaks", "use"]


def group_of(donor: str) -> str:
    for prefix, group in GROUP_PREFIX:
        if donor.startswith(prefix):
            return group
    raise ValueError(f"cannot derive group from donor id {donor!r}")


def build(data_dir: Path, root: Path) -> list[dict[str, str]]:
    rows = []
    for h5 in sorted(data_dir.glob("*_filtered_feature_bc_matrix.h5")):
        m = H5_RE.match(h5.name)
        if not m:
            raise ValueError(f"unexpected filename: {h5.name}")
        gsm, donor = m.groups()
        sample = f"{gsm}_{donor}"
        frag = data_dir / f"{sample}_atac_fragments.tsv.gz"
        peaks = data_dir / f"{sample}_atac_peaks.bed.gz"
        for p in (frag, peaks):
            if not p.exists():
                raise FileNotFoundError(f"{p} missing for {sample}")
        rows.append(
            {
                "sample": sample,
                "gsm": gsm,
                "donor": donor,
                "group": group_of(donor),
                "h5": str(h5.relative_to(root)),
                "fragments": str(frag.relative_to(root)),
                "peaks": str(peaks.relative_to(root)),
                "use": "yes" if donor in PILOT_DONORS else "no",
            }
        )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="GSE193688_RAW")
    ap.add_argument("--out", default="config/samples.csv")
    ap.add_argument(
        "--all",
        action="store_true",
        help="mark every sample use=yes (ignores PILOT_DONORS)",
    )
    args = ap.parse_args()

    root = Path.cwd()
    rows = build(root / args.data_dir, root)
    out = root / args.out

    if out.exists():
        with out.open() as fh:
            previous = {r["sample"]: r.get("use", "no") for r in csv.DictReader(fh)}
        for row in rows:
            if args.all:
                row["use"] = "yes"
            elif row["sample"] in previous:
                row["use"] = previous[row["sample"]]

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    used = [r for r in rows if r["use"] == "yes"]
    counts = {}
    for r in used:
        counts[r["group"]] = counts.get(r["group"], 0) + 1
    print(f"wrote {out} ({len(rows)} samples, {len(used)} selected: {counts})")


if __name__ == "__main__":
    main()
