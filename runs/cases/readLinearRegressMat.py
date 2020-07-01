import sys
import numpy as np
import time
import os
import petsc4py.PETSc as petsc
from sklearn import linear_model

def main(argv):

    mu = np.array([0.67, 8.37])
    inputs=argv
    M_l = float(inputs[0])
    M_r = float(inputs[1])
    dM = float(inputs[2])
    M_num = int(round((M_r-M_l)/dM))
    alfa_l = float(inputs[3])
    alfa_r = float(inputs[4])
    dalfa = float(inputs[5])
    alfa_num = int((alfa_r-alfa_l)/dalfa)
    mach = np.linspace(M_l, M_r, num= M_num+1)
    alpha = np.linspace(alfa_l, alfa_r, num = alfa_num+1)
    print(mach)
    print(alpha)
    mu_train = np.array(np.transpose(np.meshgrid(mach,alpha)).reshape\
                        (-1,2))
    print(mu_train)
    model = linear_model.LinearRegression()
    viewer = petsc.Viewer().createBinary('initialConditions.dat', mode ='r')
    A = petsc.Mat().load(viewer)
    A_dense = A.convert("dense")
    data = A_dense.getDenseArray()
    viewer.destroy()
    print(data)
    data = np.transpose(data)

    model.fit(mu_train, data)
    y0 = model.predict(np.array([mu])).flatten()
    print(y0)
    init = petsc.Vec().createWithArray(y0)

    viewerVec = petsc.Viewer().createBinary('predictedInit.dat', mode ='w')
    init.view(viewerVec)
    viewerVec.destroy()
if __name__ == '__main__':
    print("HELLLLLO\n")
    main(sys.argv[1:])
