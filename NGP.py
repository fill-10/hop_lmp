import numpy as np

# input must be numpy compatible
def NGP(ux, uy, uz, init_ux, init_uy, init_uz):
    # n dim vector:
    drsquare = (ux - init_ux )**2 +(uy - init_uy)**2 + (uz - init_uz)**2
    return  3/5* np.average(drsquare**2 ) / (np.average(drsquare))**2 -1


if  __name__ == '__main__':
    ux = np.array([1,1,1,1,1])
    uy = np.array([1,1,1,1,1])
    uz = np.array([1,1,1,1,1])
    iux = np.array([1,1,1,1,1])
    iuy = np.array([1,1,1,1,1])
    iuz = np.array([1,1,1,1,1])
    print(  NGP(ux, uy, uz, iux, iuy, iuz)  )
