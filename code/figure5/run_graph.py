"""Bounded execution of one graph, with measured resource receipt."""
import argparse,json,os,resource,subprocess,time
from pathlib import Path

def main():
    ap=argparse.ArgumentParser()
    for key in ["rscript","rlib","code","input","out","graph"]:
        ap.add_argument("--"+key,required=True)
    a=ap.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    (out/"logs").mkdir(exist_ok=True)
    _,hard=resource.getrlimit(resource.RLIMIT_STACK)
    resource.setrlimit(resource.RLIMIT_STACK,(64*1024*1024,hard))
    env=os.environ.copy()
    env.update(TRAJECTORY_RLIB=a.rlib,R_LIBS_USER=a.rlib,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
    start=time.monotonic()
    with (out/"logs"/f"{a.graph}.log").open("w") as log:
        r=subprocess.run([a.rscript,"--vanilla",a.code,a.input,a.out,a.graph],
                         env=env,stdout=log,stderr=subprocess.STDOUT)
    usage=resource.getrusage(resource.RUSAGE_CHILDREN)
    receipt=dict(graph=a.graph,exit_code=r.returncode,elapsed_seconds=time.monotonic()-start,
                 peak_RSS_KiB=usage.ru_maxrss,CPU_seconds=usage.ru_utime+usage.ru_stime,
                 stack_MiB=64,threads=1)
    (out/"logs"/f"{a.graph}_resources.json").write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt,indent=2))
    if r.returncode:
        print((out/"logs"/f"{a.graph}.log").read_text()[-12000:])
    raise SystemExit(r.returncode)

if __name__=="__main__":main()
