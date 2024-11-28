# hop_lmp

This package is a little tool to handle gromacs/lammps/pdb trajectory files, with pandas and numpy. This allows: 
1. Straight-forward and easy manimpulation of atoms, selecting by residue name, atom name, id, x, y, z, etc, within pandas data frame;
2. Faster statistical calculations for thousands of frames by vectorization (numpy array calc).

This package can be used to:
1. reads and export lammps trj, gro and pdb (unwrapped gro and pbc) into/from pandas data frame (only orthogonal box),
2. generates ions by centers of mass.
3. calculates association of the counterions (AN) and categorizes hopping types of the counterions (AN).
4. calcualtes non gaussian parameter and mean square displacement for AN. 
5. calculates van hove function and distribution of string like motions for AN.
6. calculates continous time autocorrelation function, intermittent time autocorrelation function, tau_c and tau_s.
7. Fs(q,t) and tau_q
8. saves trj in lammps format
9. way more things to do!

Not completed yet but already useful.

Howto:

pdb.main_*.py is the 'template'.
gro.main_*.py is the 'template'.

Please cite:

Luo, X.; Liu, H.; Paddison, S.J.; Molecular Dynamics Simulations of Polymerized Ionic Liquids: Mechanism of Ion Transport with Different Anions, ACS Appl. Polym. Mater. 2021, 3, 1, 141–152
