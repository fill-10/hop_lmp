import numpy as np
import re
import pandas as pd

##--- center of mass ---
##--- for unwrapped data only ---
##--- this function returns unwrapped data ---
##--- input is a pd.DataFrame with ux uy uz inside ---

def COM(atom_group, mapping_col = 'type'):
    mass_table = {'H':1.008, 'Li':6.94,'LI':6.94, 'B':10.81, 'C':12.011, 'N':14.007, 'O':15.999, 'F':18.998, 'Na':22.990, 'NA':22.990, 'Mg':4.305, 'MG':4.305,'Al':26.982, 'AL':26.982,'Si':28.085, 'SI':28.085, 'P':30.974, 'S':32.06, 'Cl':35.45, 'CL':35.45, 'K':39.098, 'Br':79.904, 'BR':79.904}
    
    # map the masses
    masses = atom_group[mapping_col].str.replace('\d+', '').map(mass_table)
    total_mass = masses.sum()
    # COM
    we_ux = np.dot(atom_group['ux'] ,masses)/total_mass
    we_uy = np.dot(atom_group['uy'] ,masses)/total_mass
    we_uz = np.dot(atom_group['uz'] ,masses)/total_mass
    
    comux =  we_ux.sum()
    comuy =  we_uy.sum()
    comuz =  we_uz.sum()
    return comux, comuy, comuz

if __name__ == '__main__':
    atoms = pd.DataFrame([  [1, '1H2', 1, 4., 0., 0. ], \
                            [2, '2H', 1, 1., 1., 2., ], \
                            [3, 'H', 1, 1., 2., 1., ]  ] , \
                            index = [19,20,21], \
                            columns = ['atom', 'type', 'mol', 'ux', 'uy','uz'] \
                        )
    atoms = pd.DataFrame([  [401, 'N', 1,  32.770,  39.280,   6.190], \
                            [402, 'C1', 1, 33.270,  38.880,   7.400], \
                            [403, 'C2', 1, 32.030,  38.310,   5.610],  \
                            [404, 'N1', 1, 31.850,  37.400,   6.640] , \
                            [405, 'C3', 1, 32.580,  37.740,   7.770] , \
                            [413, 'H',  1, 34.090,  39.300,   7.960]  , \
                            [414, 'H1', 1, 31.720,  38.160,   4.590] , \
                            [415, 'H2', 1, 32.470,  37.120,   8.660]   ],
                            index = [19,20,21,22,23,24,25,26], \
                            columns = ['atom', 'type', 'mol', 'ux', 'uy','uz'] \
                        )

    com_ux, com_uy, com_uz= COM(atoms, 'type')
    print(com_ux, com_uy, com_uz)

