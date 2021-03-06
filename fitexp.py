from scipy.optimize import curve_fit
from scipy.special import gamma
import numpy as np
import pandas as pd

# stretched exponential
def func(x, a0, a1, a2):
    return a0*np.exp(-(x/a1)**a2 )

# tau_c
try:
    df = pd.read_csv('../Tf2C2_Ctlong.dat', sep= ' ' ,escapechar='#')
    df.columns = ['t', 'Ct']
    df=df.dropna()
    popt, pcov= curve_fit(func, df['t'], df['Ct'], [1,1e4,0.1])
    perr = np.sqrt(np.diag(pcov))
    print('a0, a1, a2' , popt, '\n',pcov)

    # tau
    tau_c = popt[0]* popt[1] *gamma(1+1/popt[2])
    tau_c_low = (popt[0] - perr[0]) * (popt[1] - perr[1]) * gamma( 1 + 1/( popt[2] - perr[2]) )
    tau_c_high= (popt[0] + perr[0]) * (popt[1] + perr[1]) * gamma( 1 + 1/( popt[2] + perr[2]) )
    print('tau_c = ', tau_c)
    print( 'tau_c low and high:', (tau_c_low, tau_c_high) )
except:
    pass

# tau_s
try:
    df = pd.read_csv('../Tf2C2_St_new.dat', sep= ' ' ,escapechar='#')
    df.columns = ['t', 'St']
    df=df.dropna()
    popt, pcov= curve_fit(func, df['t'], df['St'], [1,1e4,0.1])
    perr = np.sqrt(np.diag(pcov))
    print('a0, a1, a2' , popt, '\n',pcov)

    # tau
    tau_s = popt[0]* popt[1] *gamma(1+1/popt[2])
    tau_s_low = (popt[0] - perr[0]) * (popt[1] - perr[1]) * gamma( 1 + 1/( popt[2] - perr[2]) )
    tau_s_high= (popt[0] + perr[0]) * (popt[1] + perr[1]) * gamma( 1 + 1/( popt[2] + perr[2]) )

    print('tau_s = ', tau_s)
    print( 'tau_s low and high:', (tau_s_low, tau_s_high) )
except:
    pass


# tau_q
try:
    df = pd.read_csv('../Tf2C2_fsqt_forward_3000.dat', sep= ' ' ,escapechar='#')
    df.columns = ['t', 'St']
    df=df.dropna()
    popt, pcov= curve_fit(func, df['t'], df['St'], [1,1e4,0.1])
    perr = np.sqrt(np.diag(pcov))
    print('a0, a1, a2' , popt, '\n',pcov)

    # tau
    tau_s = popt[0]* popt[1] *gamma(1+1/popt[2])
    tau_s_low = (popt[0] - perr[0]) * (popt[1] - perr[1]) * gamma( 1 + 1/( popt[2] - perr[2]) )
    tau_s_high= (popt[0] + perr[0]) * (popt[1] + perr[1]) * gamma( 1 + 1/( popt[2] + perr[2]) )

    print('tau_q = ', tau_s)
    print( 'tau_q low and high:', (tau_s_low, tau_s_high) )
except:
    pass

# tau_r1
try:
    df = pd.read_csv('../Tf2C2_RCF_1.dat', sep= ' ' ,escapechar='#')
    df.columns = ['t', 'RCF_1']
    df=df.dropna()
    popt, pcov= curve_fit(func, df['t'], df['RCF_1'], [1,1e4,0.1])
    perr = np.sqrt(np.diag(pcov))
    print('a0, a1, a2' , popt, '\n',pcov)

    # tau
    tau_C1 = popt[0]* popt[1] *gamma(1+1/popt[2])
    tau_C1_low = (popt[0] - perr[0]) * (popt[1] - perr[1]) * gamma( 1 + 1/( popt[2] - perr[2]) )
    tau_C1_high= (popt[0] + perr[0]) * (popt[1] + perr[1]) * gamma( 1 + 1/( popt[2] + perr[2]) )

    print('tau_r1 = ', tau_C1)
    print( 'tau_r1 low and high:', (tau_C1_low, tau_C1_high) )
except:
    pass

# tau_r2
try:
    df = pd.read_csv('../Tf2C2_RCF_2.dat', sep= ' ' ,escapechar='#')
    df.columns = ['t', 'RCF_2']
    df=df.dropna()
    popt, pcov= curve_fit(func, df['t'], df['RCF_2'], [1,1e4,0.1])
    perr = np.sqrt(np.diag(pcov))
    print('a0, a1, a2' , popt, '\n',pcov)

    # tau
    tau_C2 = popt[0]* popt[1] *gamma(1+1/popt[2])
    tau_C2_low = (popt[0] - perr[0]) * (popt[1] - perr[1]) * gamma( 1 + 1/( popt[2] - perr[2]) )
    tau_C2_high= (popt[0] + perr[0]) * (popt[1] + perr[1]) * gamma( 1 + 1/( popt[2] + perr[2]) )

    print('tau_r2 = ', tau_C2)
    print( 'tau_r2 low and high:', (tau_C2_low, tau_C2_high) )
except:
    pass

