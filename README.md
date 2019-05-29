# hop_lmp

This package:
1. reads in lammps trj into pandas dataframe and generate ions by centers of mass.
2. reads in ions from lammps 'fix' output of 'COM'.
3. calculate associating/coordinating counter ions/polymer chains per ion.
4. categorize hopping types of the moving ions.
5. calcualte non gaussian parameters for AN ions.
6. output stat and lammps trj.

Tested on python 3.6.

Not complete yet but already useful.

Howto:

main.py is the 'template' to import ions directly from lmpfix.
main2.py is the 'template'  to import trj and generate ions from lmptrj.

correct_image.py is the one to 'correct' the image flags for the ion atoms to make pseodu topology.

Todo:

1. fix bugs

2. van hove function
