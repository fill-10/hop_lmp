import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

def Gs_0(x, msd_t_star):
    y = (1.5/np.pi/msd_t_star)**1.5 * np.exp( -1.5 * x**2 /msd_t_star) \
        * 4*np.pi * x**2
    return y

# settings
fn_prefix = '../AmC2_'
msd_t_star = 14.52
rstep = 0.1

# read
df = pd.read_csv('../AmC2_4pi_r2_vanhove_s.dat', sep=' ', escapechar='#')
df.columns = ['r', '4pir2vhs']

#normalize

xdata = df['r']
ydata = df['4pir2vhs']/np.sum(df['4pir2vhs'])/rstep

#gauss
yGs0 =  Gs_0( xdata, msd_t_star)
# normalize
yGs0 = yGs0/np.sum(yGs0)/rstep

np.savetxt( fn_prefix+'4pi_r2_vhs_norm.dat', np.transpose( [df['r'], ydata, yGs0] ), fmt=['%f', '%f', '%f'] , header='r Gs Gs0'   ) 

import matplotlib.pyplot as plt

plt.plot(xdata, ydata, 'r-', label='Gs')
plt.plot(xdata, yGs0, 'b-', label='Gs_0')
plt.legend()
plt.xlabel('r/A')
plt.ylabel('4PIr^2*G_s(r,t*)')
plt.show()

