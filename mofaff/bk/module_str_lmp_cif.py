from ase.io import read, write, vasp
import numpy as np
import fileinput

# Convert out_sort.dump into cif using ASE
def convert_dump_cif(input_file_dump, output_file_cif):
    frames = read(input_file_dump, index=-1, format='lammps-dump-text')
    write(output_file_cif, frames, format='cif')

# Convert MOF-H2O.cif into MOF.cif
def remove_water(input_filename,output_filename,n_mof_o):
    with open(input_filename, 'r') as file:
        lines = file.readlines()

    h_indices = [i for i, line in enumerate(lines) if ' H ' in line]
    o_indices = [i for i, line in enumerate(lines) if ' O ' in line]

    total_o = len(o_indices)
    num_o = total_o - n_mof_o
    num_h = 2*num_o

    h_w = h_indices[-num_h:] if num_h <= len(h_indices) else h_indices
    o_w = o_indices[-num_o:] if num_o <= len(o_indices) else o_indices
    
    indices_w = set(h_w + o_w)
    lines_mof = [line for i, line in enumerate(lines) if i not in indices_w]

    with open(output_filename, 'w') as file:
        file.writelines(lines_mof)

