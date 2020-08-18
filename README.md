# hop_lmp

This package:
1. reads in lammps trj into pandas dataframe and generate ions by centers of mass.

3. reads pdb from gmx trjconv -pbc nojump or -pbc whole.
4. calculates associating/coordinating counter ions/polymer chains per ion.
5. categorizes hopping types of the counterions (AN).
6. calcualtes non gaussian parameter and mean square displacement for AN ions.
7. output stats and lammps trj. 
8. calculates van hove function.
9. calculates intermittent time autocorrelation function and continous time autocorrelation function.
10. calculates tau_c and tau_s from 9.
10. extract string like motion of AN ions.

Tested on python 3.6 and 3.7.

Not complete yet but already useful.

Howto:

main_*.py is the 'template' to import ions from lammps.

pdb.main_*.py is the 'template'  to import trj and generate ions from pdb from gmx trjconv.

correct_image.py is the one to 'correct' the image flags for the ion atoms to make pseodu topology.

Todo:

1. fix bugs
2. read pdb cannot fix image errors for polyatomic ions.
