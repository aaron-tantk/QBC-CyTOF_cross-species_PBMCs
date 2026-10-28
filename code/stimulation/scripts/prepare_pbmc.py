"""Make the source broad-lineage PBMC sensitivity input without matching/clustering."""
import pandas as pd
from portable import parser, paths, aligned_annotations, species_arg, species_list

EXCLUDE = {"Granulocyte/Myeloid", "Platelets", "Erythrocytes", "Unresolved"}


def main():
    p = parser(__doc__, "data", "data")
    species_arg(p)
    args = paths(p, p.parse_args())
    for sp in species_list(args):
        cyto = pd.read_csv(args.input_dir / f"cytof_{sp}_balanced.csv")
        ann = pd.read_csv(args.input_dir / f"cytof_{sp}_balanced_annotations.csv")
        aligned_annotations(cyto, ann, sp)
        cyto.loc[~ann["broad_lineage"].isin(EXCLUDE)].to_csv(
            args.output_dir / f"cytof_{sp}_balanced_pbmc_only.csv", index=False)


if __name__ == "__main__":
    main()
