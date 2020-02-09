"""
Code to parse 4u9v pdb file to detect falty protonated CYS.
Also check for SSBOND befor first ATOM line of the pdb. If detected the added H is deleted and the pdb is renumbered
"""

import os

if __name__ == '__main__':

    data_path = "./NATs_fix_cys"
    is_ssbond = False
    line_array = []
    list_files = [x for x in os.listdir(data_path) if x.endswith(".pdb")]

    for file_pdb in list_files:
        outfile = file_pdb[:-4] + "_cys.pdb"
        file_path = os.path.join("./NATs_fix_cys/" + file_pdb)
        with open(file_path, 'r') as f1, open(outfile, 'w') as f2:
            atom_cpt = 1
            for line in f1.readlines():
                if line.startswith("SSBOND"):
                    is_ssbond = True
                elif line.startswith("ATOM"):
                    line_array = line.split(" ")
                    if line_array[12].strip() == "01" or line_array[11].strip() == "01" \
                            or line_array[10].strip() == "01":
                        line_array[line_array.index("01") - 2] = "204"
                        line_array.remove("01")
                        if line_array[5].strip() == "HG":
                            atom_cpt = int(
                                line_array[next(i for i, v in enumerate(line_array) if str.isdigit(v))].strip())
                            continue
                    if line_array[12].strip() == "06" or line_array[11].strip() == "06" \
                            or line_array[10].strip() == "06":
                        line_array[line_array.index("06") - 2] = "209"
                        line_array.remove("06")
                        if line_array[5].strip() == "HG":
                            atom_cpt = int(
                                line_array[next(i for i, v in enumerate(line_array) if str.isdigit(v))].strip()) - 1
                            continue
                    line_array[next(i for i, v in enumerate(line_array) if str.isdigit(v))] = str(atom_cpt)
                    f2.write(" ".join(line_array))
                    atom_cpt += 1
                else:
                    line_array = line.split(" ")
                    f2.write(" ".join(line_array))
