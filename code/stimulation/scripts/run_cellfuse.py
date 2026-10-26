"""Run CellFuse (Kara Davis lab, Koladiya et al. 2025) exactly as published:
train supervised contrastive model on annotated CyTOF reference (18 shared anchors),
then predict broad-lineage labels for QBC query cells.
Hyperparameters from the CellFuse Quickstart vignette / DataPreparation guide.
Use --help. Requires the documented CellFuse Python source directory and a
separate compatible Python environment. No new training seed is imposed.
"""
import sys
import datetime
import os
from pathlib import Path
from portable import parser, paths

p = parser(__doc__, "cellfuse", "cellfuse")
p.add_argument("--species", choices=("human", "rhesus"), required=True)
p.add_argument("--tag", choices=("", "pbmc"), default="")
p.add_argument("--cellfuse-python", default=os.environ.get("CELLFUSE_PYTHON_DIR"),
               help="CellFuse checkout python/ directory containing python_funcs.py")
args = paths(p, p.parse_args())
if not args.cellfuse_python or not (Path(args.cellfuse_python) / "python_funcs.py").is_file():
    p.error("--cellfuse-python/CELLFUSE_PYTHON_DIR must contain python_funcs.py")
for sub in ("Saved_model", "pred"):
    (args.output_dir / sub).mkdir(exist_ok=True)
if not (args.input_dir / "data").is_dir():
    p.error("input directory must contain the prepared data/ CSV directory")

import matplotlib
matplotlib.use("Agg")
sys.path.insert(0, str(Path(args.cellfuse_python).resolve()))
import python_funcs as cf

species, tag = args.species, args.tag
ref_name = f"{species}_cytof_{tag}" if tag else f"{species}_cytof"

CF = str(args.output_dir)
INPUT = str(args.input_dir)
today = datetime.datetime.now().strftime("%Y-%m-%d")

print(f"=== CellFuse {ref_name} | model date {today} ===", flush=True)
cf.trainModel(
    dataset_name=ref_name,
    data_dir=f"{INPUT}/data",
    save_path=CF,
    device="cpu",
    lr=0.0009, margin=0.8, bs=256, epoch=50, k=5,
    min_delta=0.01, patience=5, val_step=5,
    output_dim=8, dropout_prob=0.7,
    activation_function="leaky_relu", alpha=0.01,
    cluster_column="cluster.orig",
)
print("=== training done, predicting ===", flush=True)
cf.predict_cells(
    dataset_name=ref_name,
    data_dir=f"{INPUT}/data",
    test_data_dir=f"{INPUT}/data",
    test_data=f"{species}_qbc",
    model_dir=f"{CF}/Saved_model",
    model_date=today,
    device="cpu",
    lr=0.001, margin=0.5, bs=256, epoch=50,
    knn_k=5, output_dim=8, dropout_prob=0.5,
    activation_function="leaky_relu",
    cluster_column="cluster.orig",
    output_dir=f"{CF}/pred",
)
print("=== done ===", flush=True)
