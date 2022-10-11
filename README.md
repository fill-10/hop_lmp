# hop_lmp

This package:
1. reads in lammps trj, gro and pdb (unwrapped gro and pbc) into pandas dataframe and generate ions by centers of mass.
2. calculates association of the counterions (AN) and categorizes hopping types of the counterions (AN).
3. calcualtes non gaussian parameter and mean square displacement for AN. 
4. calculates van hove function and distribution of string like motions for AN.
5. calculates continous time autocorrelation function, intermittent time autocorrelation function, tau_c and tau_s.
6. Fs(q,t) and tau_q
7. saves trj in lammps format

Tested on python 3.6.4, NumPy 1.14.2, SciPy 1.0.0, pandas 0.22.0 and matplotlib 2.0.2.

Not complete yet but already useful.

Howto:

pdb.main_*.py is the 'template'.
gro.main_*.py is the 'template'.

Please cite:

Luo, X.; Liu, H.; Paddison, S.J.; Molecular Dynamics Simulations of Polymerized Ionic Liquids: Mechanism of Ion Transport with Different Anions, ACS Appl. Polym. Mater. 2021, 3, 1, 141–152
