# hop_lmp

This package:
1. reads in lammps trj into pandas dataframe.
2. reads in ions from lammps 'fix' output of 'COM'.
3. calculate associating/coordinating counter ions/polymer chains per ion.
4. identify hopping types of the moving ions.

Tested on python 3.6.

Not complete yet but already useful.

Howto:

main.py is the 'template' to do all the above tasks.

correct_image.py is the one to correct the image flags for the ions.

Todo:

1. fix bugs

2. van hove function
