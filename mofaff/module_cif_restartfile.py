from ase.io import read, write, vasp
import numpy as np
import fileinput

# TIP4P, calculate L_w
def TIP4P_Lw(Ow, H1, H2, L_Ow=0.15):  # Consider adjusting L_Ow slightly
    Ow = np.array(Ow)
    H1 = np.array(H1)
    H2 = np.array(H2)
    # Calculate the midpoint of the H-H vector
    midpoint = (H1 + H2) / 2
    # Calculate the vector from oxygen to the midpoint
    Ow_to_mid = midpoint - Ow
    # Normalize the Ow-to-mid vector
    unit_vector = Ow_to_mid / np.linalg.norm(Ow_to_mid)
    # Calculate the Lw position
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
def convert_H2O_raspa_restartfile(input_filename,output_filename1,output_filename2,n_mof_o): #,path_to_restartfile_head):
    # Get the head of the existing restartfile
#    with open(path_to_restartfile_head, 'r') as file:
#        lines = file.readlines()
#    index = next((i for i, line in enumerate(lines) if "molecules of TIP4P" in line), None)
#    if index is not None:
#        with open(output_filename2, 'w') as file:
#            file.writelines(lines[:index])
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
    with open(output_filename1, 'w') as file:
        file.writelines(lines[:end_index])
        file.writelines(lines_h_w)
        file.writelines(lines_o_w)

    atoms = read(output_filename1)

    # Find indices of oxygen and hydrogen atoms
    O_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
    H_indices = set(i for i, atom in enumerate(atoms) if atom.symbol == 'H')
    used_hydrogens = set() # To store used hydrogen indices
    
    # Calculate coordinates of each H2O and Write to output_filename2
    with open(output_filename2, 'w') as file:
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
        print(f'Components: 1 (Adsorbates {O_indices[-1]-O_indices[0]+1}, Cations 0)',file=file)
        print('========================================================================',file=file)
        print('Components 0 (TIP4P)\n',file=file)
        print(f'Maximum-translation-change component 0: {cell_lengths[0]/10:.6f} {cell_lengths[1]/10:.6f} {cell_lengths[2]/10:.6f}',file=file)
        print('Maximum-translation-in-plane-change component 0: 0.000000,0.000000,0.000000',file=file)
        print('Maximum-rotation-change component 0: 0.523583 0.523583 0.523583',file=file)
        print('',file=file)
        print('Reactions: 0\n',file=file)
        print(f"Component: 0   Adsorbate {O_indices[-1]-O_indices[0]+1} molecules of TIP4P\n------------------------------------------------------------------------",file=file)
        # H2O atomic position: Note: The atomic order in TIP4P.def MUST Be O_w L_w H_w H_w
        for o_idx in O_indices:
            distances = []
            for h_idx in H_indices - used_hydrogens:
                dist = atoms.get_distance(o_idx, h_idx, mic=True)
                distances.append((dist, h_idx))
            nearest_two = sorted(distances, key=lambda x: x[0])[:2]
            #print(o_idx, nearest_two[0][1],*atoms[nearest_two[0][1]].position)
            used_hydrogens.update(h_idx for _, h_idx in nearest_two)
            print(f"Adsorbate-atom-position: {o_idx-O_indices[0]} 0",*atoms[o_idx].position,file=file) # O
            L_w = TIP4P_Lw(atoms[o_idx].position.tolist(),atoms[nearest_two[0][1]].position.tolist(),atoms[nearest_two[1][1]].position.tolist())
            print(f"Adsorbate-atom-position: {o_idx-O_indices[0]} 1 {L_w[0]} {L_w[1]} {L_w[2]}",file=file)
            H_count=2
            for dist, h_idx in nearest_two:
                print(f"Adsorbate-atom-position: {o_idx-O_indices[0]} {H_count}", *atoms[h_idx].position,file=file) #, Distance: {dist}") # H
                H_count = H_count + 1

        # H2O parameters
        for o_idx in O_indices: # velocity
            print(f"Adsorbate-atom-velocity: {o_idx-O_indices[0]} 0 0  0  0",file=file)
            print(f"Adsorbate-atom-velocity: {o_idx-O_indices[0]} 1 0  0  0",file=file)
            print(f"Adsorbate-atom-velocity: {o_idx-O_indices[0]} 2 0  0  0",file=file)
            print(f"Adsorbate-atom-velocity: {o_idx-O_indices[0]} 3 0  0  0",file=file)
        for o_idx in O_indices: # force
            print(f"Adsorbate-atom-force: {o_idx-O_indices[0]} 0 0  0  0",file=file)
            print(f"Adsorbate-atom-force: {o_idx-O_indices[0]} 1 0  0  0",file=file)
            print(f"Adsorbate-atom-force: {o_idx-O_indices[0]} 2 0  0  0",file=file)
            print(f"Adsorbate-atom-force: {o_idx-O_indices[0]} 3 0  0  0",file=file)
        charge_Ow, charge_Lw, charge_Hw = TIP4P_charge()
        for o_idx in O_indices: # charge
            print(f"Adsorbate-atom-charge: {o_idx-O_indices[0]} 0 {charge_Ow}",file=file)
            print(f"Adsorbate-atom-charge: {o_idx-O_indices[0]} 1 {charge_Lw}",file=file)
            print(f"Adsorbate-atom-charge: {o_idx-O_indices[0]} 2 {charge_Hw}",file=file)
            print(f"Adsorbate-atom-charge: {o_idx-O_indices[0]} 3 {charge_Hw}",file=file)
        for o_idx in O_indices: # scaling
            print(f"Adsorbate-atom-scaling: {o_idx-O_indices[0]} 0 1",file=file)
            print(f"Adsorbate-atom-scaling: {o_idx-O_indices[0]} 1 1",file=file)
            print(f"Adsorbate-atom-scaling: {o_idx-O_indices[0]} 2 1",file=file)
            print(f"Adsorbate-atom-scaling: {o_idx-O_indices[0]} 3 1",file=file)
        for o_idx in O_indices: # fixed
            print(f"Adsorbate-atom-fixed: {o_idx-O_indices[0]} 0 0  0  0",file=file)
            print(f"Adsorbate-atom-fixed: {o_idx-O_indices[0]} 1 0  0  0",file=file)
            print(f"Adsorbate-atom-fixed: {o_idx-O_indices[0]} 2 0  0  0",file=file)
            print(f"Adsorbate-atom-fixed: {o_idx-O_indices[0]} 3 0  0  0",file=file)

