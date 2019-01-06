import os, subprocess
from .spatial_slim_ts import *

# pop_width = 8.0, 
# numgens = 100, 
# sigma = 1.0, 
# K = 5.0,

def run_slim(script, seed = 23, 
             **kwargs):
    scriptbase = "_".join(script.split(".")[:-1])
    if not os.path.isdir(scriptbase):
        os.mkdir(scriptbase)
    base = os.path.join(scriptbase, 
            "_".join(["run"] + list(map(str, kwargs.values())) + [str(seed)]))
    treefile = base + ".trees"
    if os.path.isfile(treefile):
        print(treefile, "already exists.")
    else:
        logfile = base + ".log"
        slim_command = ["slim", "-s {}".format(seed)]
        slim_command += ["-d {}={}".format(k, v) for k, v in kwargs.items()]
        slim_command += ["-d \"OUTPATH='{}'\"".format(treefile), script]
        print(" ".join(slim_command))
        with open(logfile, "w") as log:
            subprocess.call(" ".join(slim_command), shell=True, stdout=log)
        if not os.path.isfile(treefile):
            raise ValueError("SLiM failed to produce output file {}".format(treefile))
    return(treefile)
