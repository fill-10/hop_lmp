from class_oneframe import oneframe
from read_1_frame import *

if __name__ == '__main__':
    pdbfilename = '../nvtlong_every100ps_0-50ns.pdb'
    f = open(pdbfilename, 'r')
    onef = oneframe()
    time, Natom, box, col, atoms, pos = read_1_pdb(f)
    
    onef.load_snap(time, box, atoms, col)
    onef.L_AN = onef.ion_gen(range(1,1+400), 1, 1, 1, [], 1, "sel[:]", 'type')
    onef.L_CT = onef.ion_gen(range(401,401+10), 2, 40, 1, [27, 1132], 8, "sel[ ((sel.index%29<=5 ) & (sel.index%29>=1)) | ( (sel.index%29>=12 ) & (sel.index%29<=14)) ]", 'type' )
    #print(onef.L_AN)
    print(onef.L_CT)

