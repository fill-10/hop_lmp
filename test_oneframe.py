from class_oneframe import oneframe
from read_1_frame import *

if __name__ == '__main__':
    pdbfilename = '../vanhove580_2200.pdb'
    f = open(pdbfilename, 'r')

    time, Natom, box, col, atoms, pos = read_1_pdb(f)
    time, Natom, box, col, atoms, pos = read_1_pdb(f)

    onef = oneframe()
    time, Natom, box, col, atoms, pos = read_1_pdb(f)
    onef.load_snap(time, box, atoms, col)
    onef.L_AN = onef.ion_gen(range(1,1+400), 1, 1, 1, [], 1, "sel[:]", 'type')
    #onef.L_CT = onef.ion_gen(range(401,401+10), 2, 40, 1, [27, 1132], 8, "sel[ ((sel.index%29<=5 ) & (sel.index%29>=1)) | ( (sel.index%29>=12 ) & (sel.index%29<=14)) ]", 'type' )
    del onef.L_atom
    onef.L_atom=[]

    secf = oneframe()
    time, Natom, box, col, atoms, pos = read_1_pdb(f)
    secf.load_snap(time, box, atoms, col)
    secf.L_AN = secf.ion_gen(range(1,1+400), 1, 1, 1, [], 1, "sel[:]", 'type')
    #secf.L_CT = secf.ion_gen(range(401,401+10), 2, 40, 1, [27, 1132], 8, "sel[ ((sel.index%29<=5 ) & (sel.index%29>=1)) | ( (sel.index%29>=12 ) & (sel.index%29<=14)) ]", 'type' )
    del secf.L_atom
    secf.L_atom = []

    f.close()
    
    onef.unwrap(onef.L_AN)
    secf.unwrap(secf.L_AN)
    fastpercent = secf.findfast(secf.L_AN, onef.L_AN, 4.6)
    onef.L_AN['fast'] = secf.L_AN['fast']
    stringlength_hist = secf.findstring(secf.L_AN, onef.L_AN, 2.5)
    print('first frame: \n', onef.L_AN.loc[onef.L_AN['fast']>0,['id','x','y','z','fast'] ])
    print('second frame: \n', secf.L_AN.loc[secf.L_AN['fast']>0,['id', 'x','y','z','fast','string'] ])
    print(fastpercent)
    print(stringlength_hist)
