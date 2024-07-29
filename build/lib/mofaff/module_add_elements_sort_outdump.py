import re
import os
import subprocess

def sort_atoms_in_snapshot(input_path, output_path, mapping):
    atom_line_regex = re.compile(r'^\d+\s+\d+')

    def sort_atoms(snapshot_atoms):
        atom_tuples = [
            (int(line.split()[1]), int(line.split()[0]), ' '.join(line.split()[2:])) for line in snapshot_atoms
        ]
        sorted_atoms = sorted(atom_tuples, key=lambda x: (x[0], x[1]))
        return [
            f"{atom_id} {atom_type} {mapping.get(atom_type, 'Unknown')} {rest_of_line}" for atom_type, atom_id, rest_of_line in sorted_atoms
        ]

    with open(input_path, 'r') as input_file:
        snapshot_atoms = []
        processing_atoms = False
        file_content = []

        for line in input_file:
            if "ITEM: ATOMS" in line:
                if snapshot_atoms:
                    file_content.extend(sort_atoms(snapshot_atoms) + [''])
                    snapshot_atoms = []
                processing_atoms = True
                columns = line.strip().split()[2:]
                if 'type' in columns and 'element' not in columns:
                    type_index = columns.index('type')
                    columns.insert(type_index + 1, 'element')
                file_content.append("ITEM: ATOMS " + ' '.join(columns))
                continue

            if processing_atoms and line.startswith("ITEM:"):
                processing_atoms = False
                if snapshot_atoms:
                    file_content.extend(sort_atoms(snapshot_atoms) + [''])
                    snapshot_atoms = []

            if processing_atoms and atom_line_regex.match(line):
                snapshot_atoms.append(line.strip())
            else:
                file_content.append(line.strip())

        if snapshot_atoms:
            file_content.extend(sort_atoms(snapshot_atoms) + [''])

    with open(output_path, 'w') as output_file:
        for item in file_content:
            output_file.write(item + '\n')

# Read the LAMMPS data from file
def clean_file(output_file_path):
    with open(output_file_path, 'r+') as f:
        lammps_out = f.read()
    with open(output_file_path, 'r+') as file:
        lines = file.readlines()
        file.seek(0)
        file.writelines(line for line in lines if line.strip())
        file.truncate()
    
