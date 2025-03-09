from ase.io import read, write, vasp
import numpy as np
import fileinput
from ase import Atoms
from ase.geometry import find_mic

# Separate H2O from mof-H2O.cif
def separate_H2O(input_filename, output_filename, n_mof_o):
    atoms=read(input_filename)
    cell_vectors = atoms.cell
    cell_lengths = atoms.cell.lengths()
    cell_angles = atoms.cell.angles()

    # Get the coordinates of H2O from input_filename
    with open(input_filename, 'r') as file:
        lines = file.readlines()

    h_indices = [i for i, line in enumerate(lines) if ' H ' in line]
    o_indices = [i for i, line in enumerate(lines) if ' O ' in line]

    total_o = len(o_indices)
    num_o = total_o - n_mof_o
    num_h = 2*num_o

    h_w = h_indices[-num_h:] if num_h <= len(h_indices) else h_indices
    o_w = o_indices[-num_o:] if num_o <= len(o_indices) else o_indices

    #indices_w = set(h_w + o_w)
    #lines_w = [line for i, line in enumerate(lines) if i in indices_w]

    lines_o_w = [line for i, line in enumerate(lines) if i in o_w]
    lines_h_w = [line for i, line in enumerate(lines) if i in h_w]

    end_index = lines.index("  _atom_site_occupancy\n") + 1
    with open(output_filename, 'w') as file:
        file.writelines(lines[:end_index])
        file.writelines(lines_h_w)
        file.writelines(lines_o_w)

# Rearrange the O and H according to their distances (distinguish H2O molecules)
def rearrange_H2O(input_filename, output_filename):
    atoms=read(input_filename)
    O_indices = [atom.index for atom in atoms if atom.symbol == 'O']
    H_indices = [atom.index for atom in atoms if atom.symbol == 'H']
    with open(input_filename, 'r') as file:
        lines = file.readlines()
    loop_list=[];coor_list=[]
    for i, line in enumerate(lines):
        if "loop_\n" in line:
            loop_list.append(i)  # loop_list[-1] is the line number of "loop_"
    for i, line in enumerate(lines):
        if i > loop_list[-1] and "_atom_site" not in line:
            coor_list.append(i)

# make sure the lines under loop_ is:
#loop_
#  _atom_site_type_symbol
#  _atom_site_label
#  _atom_site_symmetry_multiplicity
#  _atom_site_fract_x
#  _atom_site_fract_y
#  _atom_site_fract_z
#  _atom_site_occupancy

    with open(output_filename, 'w') as file:
        for i, line in enumerate(lines):
            if i < coor_list[0]:
                file.write(line)

        for O_idx in O_indices:
            d_OH_dict = {}
            print(f"O   O{O_idx}   1.0", *atoms.get_scaled_positions()[O_idx], "1.0000", file=file)
            for H_idx in H_indices:
                d_OH = atoms.get_distance(O_idx, H_idx, mic=True, vector=False)
                d_OH_dict.update({H_idx: d_OH})
            two_H = sorted(d_OH_dict, key=d_OH_dict.get)[:2]
            print(f"H   O{two_H[0]}   1.0", *atoms.get_scaled_positions()[two_H[0]], "1.0000", file=file)
            print(f"H   O{two_H[1]}   1.0", *atoms.get_scaled_positions()[two_H[1]], "1.0000", file=file)
            #print(O_idx, two_H)
#            for key in two_H:
#                del d_OH_dict[key]
            H_indices = [idx for idx in H_indices if idx not in two_H]

# define a function to set bond angle
def set_bond_angle_pbc(atom_A, atom_B, atom_C, desired_angle, cell_parameters):
    def mic(vec, cell):
        return vec - np.round(vec @ np.linalg.inv(cell)) @ cell

    def calculate_bond_length(atom1, atom2, cell):
        return np.linalg.norm(mic(atom2 - atom1, cell))

    def calculate_angle(atom1, atom2, atom3, cell):
        vec1 = mic(atom1 - atom2, cell)
        vec2 = mic(atom3 - atom2, cell)
        cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        return np.degrees(np.arccos(cos_angle))

    def rotate_vector(vec, axis, angle):
        angle_rad = np.radians(angle)
        axis = axis / np.linalg.norm(axis)
        cos_angle = np.cos(angle_rad)
        sin_angle = np.sin(angle_rad)
        cross_prod_matrix = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
        rotation_matrix = cos_angle * np.eye(3) + sin_angle * cross_prod_matrix + (1 - cos_angle) * np.outer(axis, axis)
        return rotation_matrix @ vec

    cell = np.array(cell_parameters)

    old_bond_length_BC = calculate_bond_length(atom_B, atom_C, cell)
    old_angle_ABC = calculate_angle(atom_A, atom_B, atom_C, cell)

    vec_BC = mic(atom_C - atom_B, cell)
    vec_AB = mic(atom_A - atom_B, cell)

    normal = np.cross(vec_AB, vec_BC)
    new_vec_BC = rotate_vector(vec_BC, normal, desired_angle - old_angle_ABC)
    new_atom_C = atom_B + new_vec_BC
    new_atom_C = mic(new_atom_C, cell)

    new_bond_length_BC = calculate_bond_length(atom_B, new_atom_C, cell)
    new_angle_ABC = calculate_angle(atom_A, atom_B, new_atom_C, cell)

    return new_atom_C

# Adjust bond lengths and angles of H2O to TIP4P
def adjust_H2O_TIP4P(input_filename):
    atoms=read(input_filename)
    cell_parameters = atoms.get_cell()
    O_indices = [atom.index for atom in atoms if atom.symbol == 'O']
    for O_idx in O_indices:
        H1_idx = O_idx + 1; H2_idx = O_idx + 2
        atoms.set_distance(O_idx, H1_idx, 0.9572, fix=0, mic=True)
        new_H2_position = set_bond_angle_pbc(atoms[H1_idx].position, atoms[O_idx].position, atoms[H2_idx].position, 104.52, cell_parameters)
        atoms[H2_idx].position = new_H2_position
        atoms.set_distance(O_idx, H2_idx, 0.9572, fix=0, mic=True)
    write(input_filename, atoms, format='cif')

# Implement example:
#separate_H2O('mof_H2O_md.cif','H2O_md.cif', 32)
#rearrange_H2O('H2O_md.cif','H2O_adjust.cif')
#adjust_H2O_TIP4P('H2O_adjust.cif')
