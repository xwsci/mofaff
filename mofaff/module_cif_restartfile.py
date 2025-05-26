from ase.io import read, write, vasp
import numpy as np
import fileinput
from ase.geometry import find_mic

# TIP4P, calculate L_w
def TIP4P_Lw(Ow, H1, H2, L_Ow, input_cif_filename):
    structure = read(input_cif_filename)
    cell = structure.get_cell()
    pbc = structure.get_pbc()

    Ow = np.array(Ow)
    H1 = np.array(H1)
    H2 = np.array(H2)

    H1_mic, _ = find_mic(H1 - Ow, cell, pbc)
    H2_mic, _ = find_mic(H2 - Ow, cell, pbc)
    midpoint = (H1_mic + H2_mic) / 2 + Ow
    Ow_to_mid = midpoint - Ow
    unit_vector = Ow_to_mid / np.linalg.norm(Ow_to_mid)
    Lw = Ow + L_Ow * unit_vector
    return Lw.tolist()

# TIP4P read charge of Ow, Lw, and H_w
def TIP4P_charge(file_path = "./pseudo_atoms.def"):
    with open(file_path, 'r') as file:
        for line in file:
            words = line.split()
            if "charge" in words:
                column_number = words.index("charge") + 1
            if "Ow" in words:
                charge_Ow = words[column_number - 1]
            if "Lw" in words:
                charge_Lw = words[column_number - 1]
            if "Hw" in words:
                charge_Hw = words[column_number - 1]
    return charge_Ow, charge_Lw, charge_Hw

# Convert MOF-H2O.cif into H2O.dat -- gRASPA restartfile
def convert_H2O_raspa_restartfile(input_filename,output_filename): #,path_to_restartfile_head):
    atoms = read(input_filename)
    cell_vectors = atoms.cell
    cell_lengths = atoms.cell.lengths()
    cell_angles = atoms.cell.angles()

    # Find indices of oxygen and hydrogen atoms
    O_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
    H_indices = set(i for i, atom in enumerate(atoms) if atom.symbol == 'H')
    used_hydrogens = set() # To store used hydrogen indices
    
    # Calculate coordinates of each H2O and Write to output_filename
    with open(output_filename, 'w') as file:
        print("Cell info:\n========================================================================",file=file)
        print('number-of-unit-cells: 1 1 1',file=file)
        print('unit-cell-vector-a:',*[format(x, '.6f') for x in cell_vectors[0]],file=file)
        print('unit-cell-vector-b:',*[format(x, '.6f') for x in cell_vectors[1]],file=file)
        print('unit-cell-vector-c:',*[format(x, '.6f') for x in cell_vectors[2]],file=file)
        print('',file=file)
        print('cell-vector-a:',*[format(x, '.6f') for x in cell_vectors[0]],file=file)
        print('cell-vector-b:',*[format(x, '.6f') for x in cell_vectors[1]],file=file)
        print('cell-vector-c:',*[format(x, '.6f') for x in cell_vectors[2]],file=file)
        print('',file=file)
        print('cell-lengths:',*[format(x, '.6f') for x in cell_lengths],file=file)
        print('cell-angles:',*[format(x, '.6f') for x in cell_angles],file=file)
        print('\n',file=file)
        print('Maximum changes for MC-moves:\n========================================================================',file=file)
        print('Maximum-volume-change: 0.006250\nMaximum-Gibbs-volume-change: 0.025000\nMaximum-box-shape-change: 0.100000 0.100000 0.100000, 0.100000 0.100000 0.100000, 0.100000 0.100000 0.100000\n\n',file=file)
        print('Acceptance targets for MC-moves:\n========================================================================',file=file)
        print('Target-volume-change: 0.500000\nTarget-box-shape-change: 0.500000\nTarget-Gibbs-volume-change: 0.500000',file=file)
        print('\n',file=file)
        print(f'Components: 1 (Adsorbates {len(O_indices)}, Cations 0)',file=file)
        print('========================================================================',file=file)
        print('Components 0 (TIP4P)\n',file=file)
        print(f'Maximum-translation-change component 0: {cell_lengths[0]/10:.6f} {cell_lengths[1]/10:.6f} {cell_lengths[2]/10:.6f}',file=file)
        print('Maximum-translation-in-plane-change component 0: 0.000000,0.000000,0.000000',file=file)
        print('Maximum-rotation-change component 0: 0.523583 0.523583 0.523583',file=file)
        print('',file=file)
        print('Reactions: 0\n',file=file)
        print(f"Component: 0   Adsorbate {len(O_indices)} molecules of TIP4P\n------------------------------------------------------------------------",file=file)
        # H2O atomic position: Note: The atomic order in TIP4P.def MUST Be O_w L_w H_w H_w
        for i, o_idx in enumerate(O_indices):
            distances = []
            for h_idx in H_indices - used_hydrogens:
                dist = atoms.get_distance(o_idx, h_idx, mic=True)
                distances.append((dist, h_idx))
            nearest_two = sorted(distances, key=lambda x: x[0])[:2]
            #print(o_idx, nearest_two[0][1],*atoms[nearest_two[0][1]].position)
            used_hydrogens.update(h_idx for _, h_idx in nearest_two)
            print(f"Adsorbate-atom-position: {i} 0",*atoms[o_idx].position,file=file) # O
            L_w = TIP4P_Lw(atoms[o_idx].position.tolist(),atoms[nearest_two[0][1]].position.tolist(),atoms[nearest_two[1][1]].position.tolist(), L_Ow=0.15,input_cif_filename = input_filename)
            print(f"Adsorbate-atom-position: {i} 1 {L_w[0]} {L_w[1]} {L_w[2]}",file=file)
            H_count=2
            for dist, h_idx in nearest_two:
                print(f"Adsorbate-atom-position: {i} {H_count}", *atoms[h_idx].position,file=file) #, Distance: {dist}") # H
                H_count = H_count + 1

        # H2O parameters
        for i, o_idx in enumerate(O_indices): # velocity
            print(f"Adsorbate-atom-velocity: {i} 0 0  0  0",file=file)
            print(f"Adsorbate-atom-velocity: {i} 1 0  0  0",file=file)
            print(f"Adsorbate-atom-velocity: {i} 2 0  0  0",file=file)
            print(f"Adsorbate-atom-velocity: {i} 3 0  0  0",file=file)
        for i, o_idx in enumerate(O_indices): # force
            print(f"Adsorbate-atom-force: {i} 0 0  0  0",file=file)
            print(f"Adsorbate-atom-force: {i} 1 0  0  0",file=file)
            print(f"Adsorbate-atom-force: {i} 2 0  0  0",file=file)
            print(f"Adsorbate-atom-force: {i} 3 0  0  0",file=file)
        charge_Ow, charge_Lw, charge_Hw = TIP4P_charge()
        for i, o_idx in enumerate(O_indices): # charge
            print(f"Adsorbate-atom-charge: {i} 0 {charge_Ow}",file=file)
            print(f"Adsorbate-atom-charge: {i} 1 {charge_Lw}",file=file)
            print(f"Adsorbate-atom-charge: {i} 2 {charge_Hw}",file=file)
            print(f"Adsorbate-atom-charge: {i} 3 {charge_Hw}",file=file)
        for i, o_idx in enumerate(O_indices): # scaling
            print(f"Adsorbate-atom-scaling: {i} 0 1",file=file)
            print(f"Adsorbate-atom-scaling: {i} 1 1",file=file)
            print(f"Adsorbate-atom-scaling: {i} 2 1",file=file)
            print(f"Adsorbate-atom-scaling: {i} 3 1",file=file)
        for i, o_idx in enumerate(O_indices): # fixed
            print(f"Adsorbate-atom-fixed: {i} 0 0  0  0",file=file)
            print(f"Adsorbate-atom-fixed: {i} 1 0  0  0",file=file)
            print(f"Adsorbate-atom-fixed: {i} 2 0  0  0",file=file)
            print(f"Adsorbate-atom-fixed: {i} 3 0  0  0",file=file)


#L = TIP4P_Lw(Ow, H1, H2, L_Ow=0.15,input_cif_filename = 'H2O_adjust.cif')
#convert_H2O_raspa_restartfile('H2O_adjust.cif','restart')
