# python long_axis_vector.py <input_pdb>
import sys

import pandas as pd

pdb = sys.argv[1]

# print(pdb)

def lav(pdb):
    """Return long-axis vector from alpha carbons of first/last residue."""

    x1 = 0
    y1 = 0
    z1 = 0
    x2 = 0
    y2 = 0
    z2 = 0
    # Residue numbering in the input PDB starts at 2 for this workflow.
    start_res = str(2)
    # End residue is the last residue index minus one (final row is excluded).
    end_res = str(len(pdb.groupby(by=4))-1)

    for i in pdb.groupby(by=4):
        if i[0] == start_res:
            print(i[0], i[1])
            ind = i[1][i[1][2] == 'CA'].index[0]
            # Values are string tokens; eval preserves existing numeric parsing.
            x1 = eval(i[1][5][ind])
            y1 = eval(i[1][6][ind])
            z1 = eval(i[1][7][ind])
        if i[0] == end_res:
            print(i[0], i[1])
            ind = i[1][i[1][2] == 'CA'].index[0]
            # Values are string tokens; eval preserves existing numeric parsing.
            x2 = eval(i[1][5][ind])
            y2 = eval(i[1][6][ind])
            z2 = eval(i[1][7][ind])

    # print(x1, y1, z1)
    # print(x2, y2, z2)
    x = x2-x1
    y = y2-y1
    z = z2-z1
    # print(x, y, z)
    return x, y, z

def parse_pdb_lines(pdb_path):
    """Load and split PDB lines into token lists (excluding trailing TER/END)."""
    with open(pdb_path, 'r') as f:
        lines = f.readlines()

    tokens = []
    for line in lines:
        tokens.append(line.strip().split())

    # Remove the final 'TER' and 'END' records.
    return tokens[:-2]

def main():
    # prep pdb
    l1 = parse_pdb_lines(pdb)

    # pdb to pandas df
    pdb_df = pd.DataFrame(l1)

    x, y, z = lav(pdb_df)

    with open(pdb[:-4] + '_LA.txt', 'w') as f:
        f.write('PDB File: ' + pdb + '\n')
        f.write('Long axis vector (x y z) considering alpha carbons of 1st and last residue:\n')
        f.write("{{{:.4f} {:.4f} {:.4f}}}\n".format(x, y, z))


if __name__ == "__main__":
    main()
