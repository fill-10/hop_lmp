import numpy as np 

def MSD(ux, uy, uz, init_ux, init_uy, init_uz):
    # n dim vector:
    drsquare = (ux - init_ux )**2 +(uy - init_uy)**2 + (uz - init_uz)**2
    return np.average(drsquare)

