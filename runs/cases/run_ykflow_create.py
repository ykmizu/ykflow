
import numpy as np
import os, sys

if __name__ == '__main__':
    inputs=sys.argv[1:]

    num = int(inputs[0])

    folder = "ykflow_2D_ROM_"+str(num)
    os.mkdir(folder)
    os.system("cp createParam.py naca.inp naca.job naca.geom naca_next_5.gri "+ str(folder))
    os.system("cp naca_FOM_dat/naca_* "+str(folder))
    os.chdir(os.getcwd()+"/"+folder)
    with open("YK_2D_"+str(num)+".batch", "w+") as inputWrite:
        inputWrite.write("#!/bin/bash\n")
        inputWrite.write("#SBATCH -o ykflow_2D_"+str(num)+"_ID_%j.out\n")
        inputWrite.write("#SBATCH --nodes=1\n");
        inputWrite.write("#SBATCH --exclude=unogpu[1-25],unoqs[1-8]\n");
        inputWrite.write("#SBATCH --job-name=yk_navier_stokes\n");
        inputWrite.write("#SBATCH --time=48:00:00\n");
        inputWrite.write("#SBATCH --account=FY140174\n");
        inputWrite.write("#SBATCH --partition=batch\n");
        inputWrite.write("#SBATCH --qos=normal\n");
        inputWrite.write("\n")

        inputWrite.write("../../../bin/yk_STLSPG -job naca -input naca\n");

    inputWrite.close()
    os.system("sbatch YK_2D_"+str(num)+".batch")
