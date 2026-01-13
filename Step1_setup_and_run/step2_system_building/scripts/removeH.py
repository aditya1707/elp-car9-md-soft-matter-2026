import sys

fn = sys.argv[1]

with open(fn, 'r') as file, open(fn.split('.')[0]+'_noH.pdb','w') as out:
    for line in file.readlines():
        l = line.split()
        if l[0] == 'ATOM':
            if len(l)>1:
                if l[3] == 'ACE' or l[3] == 'NME':
                    out.write(line)

                elif l[2][0] != 'H':
                    out.write(line)

        elif l[0] == 'TER' or l[0] == 'END':
            out.write(line)
