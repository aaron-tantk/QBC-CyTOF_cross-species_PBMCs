"""Shared path, input-alignment and public donor-label utilities (no RNG changes)."""
import argparse
import os
from pathlib import Path
import re


def parser(description, input_subdir="data", output_subdir="results"):
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--project-root", default=os.environ.get("QBC_PROJECT_ROOT"))
    p.add_argument("--input-dir", help=f"Overrides PROJECT_ROOT/{input_subdir}")
    p.add_argument("--output-dir", help=f"Overrides PROJECT_ROOT/{output_subdir}")
    p.set_defaults(input_subdir=input_subdir, output_subdir=output_subdir)
    return p


def paths(p, args):
    root = Path(args.project_root).expanduser().resolve() if args.project_root else None
    for key, sub in (("input_dir", args.input_subdir), ("output_dir", args.output_subdir)):
        value = getattr(args, key)
        if value is None and root is None:
            p.error(f"--{key.replace('_', '-')} or --project-root/QBC_PROJECT_ROOT is required")
        setattr(args, key, Path(value).expanduser().resolve() if value else root / sub)
    if not args.input_dir.is_dir():
        p.error(f"input directory does not exist: {args.input_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    return args


def neutral_donors(series, species, modality):
    """Relabel after source row selection, never before seeded group sampling.

    Sorted unique input strings define neutral indices. No identity lookup is
    returned, logged, or written. All linked tables must contain the full donor
    universe (or already use neutral labels) before selecting subsets.
    """
    if series.isna().any():
        raise ValueError("donor column contains missing values")
    values = series.astype(str)
    prefix = f"{species}_{modality}_"
    neutral = values.str.fullmatch(re.escape(prefix) + r"\d{2,}")
    if neutral.all():
        return values
    if neutral.any():
        raise ValueError("mixed neutral and original donor labels")
    mapping = {v: f"{prefix}{i:02d}" for i, v in enumerate(sorted(values.unique()), 1)}
    return values.map(mapping)


def aligned_annotations(cyto, ann, species):
    """Verify row-wise source annotation alignment before any shuffle/filter."""
    if len(cyto) != len(ann):
        raise ValueError(f"{species}: annotation row count mismatch")
    for table in (cyto, ann):
        table["donor"] = neutral_donors(table["donor"], species, "cytof")
    for col in ("donor", "stimulation"):
        if not (cyto[col].astype(str).to_numpy() == ann[col].astype(str).to_numpy()).all():
            raise ValueError(f"{species}: annotation {col} row-order mismatch")
    for col in ("CyTOF_annotation", "broad_lineage"):
        if col not in ann or ann[col].isna().any():
            raise ValueError(f"{species}: missing source {col}")


def species_arg(p):
    p.add_argument("--species", choices=("human", "rhesus", "both"), default="both")


def species_list(args):
    return ("human", "rhesus") if args.species == "both" else (args.species,)
