def read_1_lmp(f, start=0, dim = 3):
    if start>0:
        f.seek(start)
    ###--- read controlling flags
    readtime = 0
    readnumber = 0
    readbox = 0
    readatom = 0

    item_atoms = []
    item_time = 0.0
    item_box = []
    Natom = 0
    atom_count  = 0
    col_dict = []
    # read head info
    rawline = '\n'
    while rawline !='':
        rawline = f.readline()
        cline = rawline.split()
        if cline[0] == 'ITEM:':  # control flags, find which to read
            if cline[1] == 'TIMESTEP':
                readtime = 1
                readnumber = 0
                readbox = 0
                readatom = 0
                continue
            if cline[1] == 'NUMBER':
                readnumber = 1
                readtime = 0
                readbox = 0
                readatom = 0
                continue
            if cline[1] == 'BOX':
                readbox = dim # as counter
                readtime = 0
                readnubmer = 0
                readatom = 0
                continue
            if cline[1] == 'ATOMS':
                readtime = 0
                readnumber = 0
                readbox = 0
                readatom = 1
                col_dict = cline[2:]
                if Natom>0:
                    atom_count = Natom
                    continue
                else:
                    raise Exception(f.tell(), 'Natom not given before reading atoms')
        # read info
        if readtime == 1 :
            item_time = float(cline[0])
            readtime = 0
        if readnumber == 1 :
            Natom = int(cline[0])
            atom_count = Natom
            readnumber = 0
        if readbox>0:
            # convert data type
            for k in range(0, len(cline)):
                cline[k] = float(cline[k])
            item_box += [cline]
            readbox -= 1
        if readatom and atom_count >0:
            # convert data types
            for k in range(0, len(col_dict)):
                if col_dict[k] in ['id', 'mol','type', 'ix', 'iy', 'iz']:
                    cline[k] = int(cline[k])
                elif col_dict[k] in ['q', 'x', 'y', 'z', 'ux', 'uy', 'uz']:
                    cline[k] = float(cline[k])
                elif col_dict[k] == 'element':
                    pass
                else:
                    cline[k] = float(cline[k])
            item_atoms += [cline]
            atom_count -= 1

        if readatom and atom_count<=0:
            return item_time, Natom, item_box, col_dict, item_atoms, f.tell()

if __name__ == '__main__' :
    filename = '../400K_corrected_lmp_B/VImC4_every100ps_0-50ns.lammpstrj'
    f = open(filename, 'r')
    time, Natom, box, col, atoms, pos = read_1_lmp(f)
    time, Natom, box, col, atoms, pos = read_1_lmp(f)
    print(time)
    print(box)
    print(col)
    print(atoms[1])
