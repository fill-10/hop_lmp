import numpy as np
def VANHOVE_S(ux, uy, uz, init_ux, init_uy, init_uz, maxdist =25.0, accuracy = 0.1):
    # input should be numpy.array or pandas.Series
    # N dim vector (N atoms)
    distances = np.sqrt( (ux-init_ux)**2 + (uy-init_uy)**2 + (uz-init_uz)**2 )    # element wise square: **2 or np.square(), element wise sqrt: np.sqrt()
    mybins = np.arange(0, maxdist, accuracy)
    return np.histogram(distances, bins=mybins, density=True)[0]*accuracy, (mybins[:-1]+accuracy/2) # normalized hist, center of each bin

def VANHOVE_D(x, y, z, init_x, init_y, init_z, boxedges, maxdist =25.0, accuracy = 0.1):
    #vectorize
    Natom = len(x)
    L_dx = np.array([])
    L_dy = np.array([])
    L_dz = np.array([])
    distances = np.array([])
    
    for i in range(1, Natom):
        shift_x = np.append(x[i:], x[:i])
        shift_y = np.append(y[i:], y[:i])
        shift_z = np.append(z[i:], z[:i])
        dx = np.abs(shift_x - init_x)
        dy = np.abs(shift_y - init_y)
        dz = np.abs(shift_z - init_z)
        L_dx = np.append(L_dx, dx)
        L_dy = np.append(L_dy, dy)
        L_dz = np.append(L_dz, dz)

    dx_pbc = np.min( [L_dx, boxedges[0] - L_dx ], axis=0 )
    dy_pbc = np.min( [L_dy, boxedges[1] - L_dy ], axis=0 )
    dz_pbc = np.min( [L_dz, boxedges[2] - L_dz ], axis=0 )
   
    distances = np.sqrt(dx_pbc**2 + dy_pbc**2 + dz_pbc**2)
    
    mybins = np.arange(0, maxdist, accuracy)

    return np.histogram(distances, bins=mybins, density=True)[0]*accuracy, (mybins[:-1]+accuracy/2) 

if __name__=='__main__':
    y = np.array([-1,-0.5,4])
    y0 = np.array([-1,0.5,5])
    # 1.5, 6, 0.5,5.5, 5, 3.5
    
    x = x0 =np.array([0,0,0])
    
    z =  np.array([1,0,0])
    z0 = z+1
    #print(np.sqrt((x-x0)**2 +(y-y0)**2 +(z-z0)**2))
    #print(VANHOVE_S(x,y,z,x0,y0,z0,6.9,0.5))
    print(VANHOVE_D(x,y,z,x0,y0,z0,[6,6,6],6.9,0.5))


