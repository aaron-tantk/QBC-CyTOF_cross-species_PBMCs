"""Validate saved data, index alignment and archived result summaries."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def sha256(path):
    result=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b""): result.update(block)
    return result.hexdigest()


def validate(data,repository,work):
    data,repository,work=map(Path,[data,repository,work])
    expected={"primary":{"human":20862,"rhesus":14414},
              "whole_blood_sensitivity":{"human":20550,"rhesus":13501}}
    counts={}
    for design,species_counts in expected.items():
        for species,n in species_counts.items():
            folder=data/"processed"/design
            md=pd.read_csv(folder/f"{species}_metadata.csv.gz")
            assert len(md)==n and md.cell_id.is_unique
            for stage in ["pre_joint_CLR","joint_CLR"]:
                x=pd.read_csv(folder/f"{species}_{stage}.csv.gz",float_precision="round_trip")
                assert x.shape==(n,175) and x.cell_id.equals(md.cell_id)
                assert sum(c.endswith("-QBC") for c in x)==134
                assert sum(c.endswith("-CyTOF") for c in x)==40
                assert np.isfinite(x.iloc[:,1:].to_numpy()).all()
            disagreement=int((md.source_reference_population!=md.audit_plurality_population).sum())
            counts[f"{species}_{design}"]={"n_rows":n,"source_plurality_disagreements":disagreement}
    excluded={"Granulocyte/Myeloid","Platelets","Erythrocytes","Unresolved"}
    for species in ["human","rhesus"]:
        upstream=data/"upstream_inputs"
        ann=pd.read_csv(upstream/f"cytof_{species}_balanced_annotations.csv.gz")
        qbc=pd.read_csv(upstream/f"qbc_{species}_30k.csv.gz")
        for design,prefix in [("primary","pbmc_"),("whole_blood_sensitivity","")]:
            matches=pd.read_csv(data/"matching"/f"{prefix}{species}_matches.csv")
            md=pd.read_csv(data/"processed"/design/f"{species}_metadata.csv.gz")
            qi=pd.Series(np.arange(len(qbc))).sample(frac=1,random_state=42).to_numpy()
            ci=np.flatnonzero(~ann.broad_lineage.isin(excluded)) if prefix else np.arange(len(ann))
            ci=pd.Series(ci).sample(frac=1,random_state=42).to_numpy()
            if not prefix: ci=ci[:30000]
            np.testing.assert_array_equal(qi[matches.row_idx_0.to_numpy(int)],md.qbc_row_id)
            np.testing.assert_array_equal(ci[matches.row_idx_1.to_numpy(int)],md.cytof_row_id)
    adapter=load("release_adapter",repository/"code/figure5/prepare_source.py")
    analysis=load("release_analysis",repository/"code/figure5/analysis_gmcsf.py")
    adapter.prepare(data/"processed/primary",work/"source")
    analysis.prepare(work/"source",work)
    # Data graphs are immutable inputs; only local work outputs are written.
    import shutil
    shutil.copytree(data/"figure5/graphs",work/"graphs",dirs_exist_ok=True)
    analysis.summarize(work)
    compared=0
    for path in (data/"results/figure5").glob("*.csv"):
        output=work/"summary"/path.name
        if not output.exists(): output=work/"inputs"/path.name
        pd.testing.assert_frame_equal(pd.read_csv(output),pd.read_csv(path),
                                      check_exact=False,rtol=1e-10,atol=1e-12)
        compared+=1
    network=load("release_networks",repository/"code/figure6/summarize_networks.py")
    network.summarize(data/"figure6/bootstrap_support",work/"networks")
    for path in (data/"results/figure6").glob("*.csv"):
        pd.testing.assert_frame_equal(pd.read_csv(work/"networks"/path.name),pd.read_csv(path))
    receipt={"status":"PASS","processed_data":counts,"matching_index_alignment":"PASS",
             "figure5_numerical_tables":compared,"network_support_tables":36,
             "network_summary_populations":12,"scientific_refitting_performed":False}
    (work/"validation.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))
    return receipt


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data",required=True,type=Path)
    parser.add_argument("--repository",required=True,type=Path)
    parser.add_argument("--work",required=True,type=Path)
    args=parser.parse_args()
    validate(args.data,args.repository,args.work)
