import os
import numpy as np
from shutil import copyfile

Nb = 50
Nc = 20
dt = 0.0001
num_threads = 16
restart = 0
trace = range(10,20)
total_strain = 10.0
flowrate = 0.00096
nt_avg  = 8
bin_size = 2.0
fixSeed = 0
iterstart = 0
itermax = 2
c_norm = 0.1
for tr in trace:
  basename = str("%.3lf"%c_norm)+"_%d"%Nb+str("_%d"%tr)
  with open("qsub/disp_wi_aa_c_aa"+basename+".qsub","w") as f:     
        f.write("#!/bin/bash\n")
        f.write("#$ -cwd\n") # set the working directory to the current directory
        f.write("#$ -j y\n") # join output and error files
        f.write("#$ -S /bin/bash\n")
        f.write("#$ -p -1\n") # lower priority. Increase for higher priority
        f.write("#$ -q parallel2.q\n")
        f.write("#$ -pe default %d\n"%num_threads)
        # f.write("#$ -l h=\"compute-3-*\"\n") # specify the node, modify to compute-X-X.local 
        # This set ups an array job. Essentially we are running 100 identical tasks, and the only change
        # is the SGE_TASK_ID, which we are inputting as the trajectory number
        # The random seed will be different in each task, so this provides an ensemble of molecules
        f.write("#$ -o ./qsub/cluster_out\n") # set the output file directory
        f.write("\n")
        f.write("~/CApack/src/test_nov/DJ_meth/HI_jan_14b/run.out -b %d -a %d -c %lf -f %lf -e %lf -d %lf -n %d -r %d -s %d -g %lf -z %d -i %d -m  %d -p %d\n"%(Nb,Nc,c_norm,flowrate,total_strain,dt,num_threads,restart,nt_avg,bin_size,fixSeed,iterstart,itermax,tr))
## end
