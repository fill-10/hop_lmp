import numpy as np
from numpy.linalg import norm

def pbc_dist(coor_1, coor_2, boxedges,dimension=3):
    # boxedges is a 3-d vector of box size, always positive. type: numpy.ndarray
    
    ##--- calculate pbc distance ---
    dv = abs(coor_1 - coor_2)
    
    for i in range(0, dimension):
        ##--- wrap  ---
        while dv[i] >= boxedges[i]:
            dv[i] -= boxedges[i]
        ##---
        dvifpbc = abs(boxedges - dv)
        dv[i] = min(dv[i], dvifpbc[i])

    return norm(dv)

def ifconn(coor_1, coor_2, boxedges, r_cut, dimension=3):
    dv = np.absolute( np.subtract(coor_1, coor_2) )
    for i in range(0, dimension):
        ##--- wrap distance by pbc ---
        while dv[i] >= boxedges[i]:
            dv[i] -= boxedges[i]
        ##-----------------
        dv[i] = min( dv[i], np.absolute( dv[i] - boxedges[i] ) )
        if dv[i] > r_cut:
            return False
    if np.dot(dv,dv) > r_cut**2:
        return False
    else:
        return True
   
# to test:
if __name__ == '__main__':
    boxlo = np.array([0.0,0.0,0.0])
    boxhi = np.array([2., 2., 1.])
    boxed = np.subtract(boxhi,boxlo)
    coor1 = np.array([5.9,0.01,0.02])
    coor2 = np.array([0.5,1.9,.2])
    from time import perf_counter
    start = perf_counter()
    for i in range(0,200000):
        so= ifconn( coor1,coor2, boxed, 0.82)
    finish = perf_counter()
    print( finish -start )
    print(so)
