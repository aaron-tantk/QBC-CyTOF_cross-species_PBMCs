"""Adapt the released primary joint-CLR tables to the Figure 5 input schema."""
import argparse
from pathlib import Path
import pandas as pd


def prepare(processed, output):
    processed, output = Path(processed), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for species in ["human", "rhesus"]:
        md = pd.read_csv(processed/f"{species}_metadata.csv.gz")
        x = pd.read_csv(processed/f"{species}_joint_CLR.csv.gz", float_precision="round_trip")
        assert x.cell_id.equals(md.cell_id) and x.shape[1] == 175
        assert md.design.eq("primary").all()
        md = md.rename(columns={"source_reference_population": "population",
                                "source_reference_broad_lineage": "broad_lineage"})
        fields = ["cell_id", "reference_id", "species", "qbc_donor", "cytof_donor",
                  "condition", "population", "broad_lineage", "qbc_row_id", "cytof_row_id"]
        # No historical all-condition score or combined-stimulus label is carried over.
        md[fields].to_csv(output/f"{species}_metadata.csv.gz", index=False)
        x.to_csv(output/f"{species}_matrix.csv.gz", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.processed, args.output)
