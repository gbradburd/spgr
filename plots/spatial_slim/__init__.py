import os, subprocess
from .spatial_slim_ts import *

def run_slim(script = "flat_map.slim",
             seed = 23, 
             pop_width = 8.0, 
             numgens = 100, 
             sigma = 1.0, **kwargs):
    scriptbase = "_".join(script.split(".")[:-1])
    if not os.path.isdir(scriptbase):
        os.mkdir(scriptbase)
    base = os.path.join(scriptbase, "run_{}_{}_{}".format(sigma, pop_width, numgens, seed))
    treefile = base + ".trees"
    if os.path.isfile(treefile):
        print(treefile, "already exists.")
    else:
        logfile = base + ".log"
        slim_command = ["slim", 
                        "-s {}".format(seed), 
                        "-d W={}".format(pop_width),
                        "-d NUMGENS={}".format(numgens),
                        "-d SIGMA={}".format(sigma),
                        "-d \"OUTPATH='{}'\"".format(treefile),
                        "flat_map.slim"]
        print(" ".join(slim_command))
        with open(logfile, "w") as log:
            subprocess.call(" ".join(slim_command), shell=True, stdout=log)
    return(treefile)
