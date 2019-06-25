import numpy as np
from pbc_dist import pbc_dist
def VANHOVE_S(ux, uy, uz, init_ux, init_uy, init_uz, maxdist =25.0, accuracy = 0.1):
    # input should be numpy.array or pandas.Series
    # N dim vector (N atoms)
    distances = np.sqrt( (ux-init_ux)**2 + (uy-init_uy)**2 + (uz-init_uz)**2 )    # element wise square: **2 or np.square(), element wise sqrt: np.sqrt()
    mybins = np.arange(0, maxdist, accuracy)
    return np.histogram(distances, bins=mybins, density=True)[0]*accuracy, (mybins[:-1]+accuracy/2) # normalized hist, center of each bin

def VANHOVE_D(x, y, z, init_x, init_y, init_z, boxedges, maxdist =25.0, accuracy = 0.1):
    Natom = min(len(x), len(y), len(z), len(init_x), len(init_y), len(init_z) )
    distances = []
    for i in range(0, Natom):
        for j in range(0, Natom):
            if j == i :
                continue
            else:
                distances.append( pbc_dist( np.array([ x[i],y[i],z[j] ]), np.array([init_x[j],init_y[j],init_z[j]]), boxedges)  )
                #distances.append(np.sqrt( (ux[i]-init_ux[j])**2 + (uy[i]-init_uy[j])**2 + (uz[i]-init_uz[j])**2 ) )
    mybins = np.arange(0, maxdist, accuracy)
    return np.histogram(distances, bins=mybins, density=True)[0]*accuracy, (mybins[:-1]+accuracy/2) 

if __name__=='__main__':
    y = np.array([-1,-0.5,4])
    y0 = np.array([-1,0.5,5])
    # 1.5, 6, 0.5,5.5, 5, 3.5
    
    x = x0 =np.array([0,0,0])
    
    z =  np.array([-0,0,0])
    z0 = z
    print(np.sqrt((x-x0)**2 +(y-y0)**2 +(z-z0)**2))
    print(VANHOVE_S(x,y,z,x0,y0,z0,6.9,0.5))
    print(VANHOVE_D(x,y,z,x0,y0,z0,[6,6,6],6.9,0.5))
    print(VANHOVE_D(x,y,z,x0,y0,z0,6.9,0.5))


