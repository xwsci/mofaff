import shutil
import re
import numpy as np

def get_type(element, mapping):
    for key, value in mapping.items():
        if value == element:
            return key
    return None

def raspa_to_lmp(raspa_data,output,N_mof_O, mapping, hybrid_FF):
    shutil.copy(raspa_data, output)
    with open(output, 'r') as file:
        lines = file.readlines()

    with open(output, 'w') as file:
        between_masses_atoms = False
        atoms_section = False
        for line in lines:
            if 'Masses' in line:
                between_masses_atoms = True
                continue
            if 'Atoms' in line:
                between_masses_atoms = False
                atoms_section = True
                file.write(line)
                continue
            if between_masses_atoms:
                continue
            if not re.search(r'(bond|angle|dihedral|improper|C_co2|O_co2|N_n2|N_com|Lw|Zr|atom types)', line):
                file.write(line)

    with open(output, 'r') as file:
        lines = file.readlines()

    found_atoms_1 = False
    found_atoms_2 = False
    count_1=0
    count_2=1
    n_mol=1

    with open(output, 'w') as file:
        for line in lines:
            if 'Atoms' in line:
                found_atoms_1 = True
            if found_atoms_1 and len(line.split()) > 3:
                count_1 = count_1 + 1

        for line in lines:
            if 'atoms' in line:
                file.write(f'{count_1} atoms\n{len(mapping)} atom types\n')
            else:
                for key, element in mapping.items():
                    line = re.sub(f'# {element}(?!\w)', f'#{element} {key}', line)
                line = line.replace('# Hw', f'#Hw {get_type("H", mapping)}')
                line = line.replace('# Ow', f'#Ow {get_type("O", mapping)}')
                if 'Atoms' in line:
                    found_atoms_2 = True
                if found_atoms_2 and len(line.split()) > 3:
                    parts = line.split()
                    if parts[-2] != '#Ow' and parts[-2] != '#Hw':
                        line = f"{count_2} {n_mol} {parts[-1]} {parts[4]} {parts[5]} {parts[6]} {parts[-2]}\n"                    
                    elif parts[-2] == '#Ow':
                         n_mol += 1
                         line = f"{count_2} {n_mol} {parts[-1]} {parts[4]} {parts[5]} {parts[6]} {parts[-2]}\n"                         
                    elif parts[-2] == '#Hw':
                         line = f"{count_2} {n_mol} {parts[-1]} {parts[4]} {parts[5]} {parts[6]} {parts[-2]}\n"
                    count_2 = count_2 + 1
                file.write(line)

    if hybrid_FF == True:    
        with open(output, 'a') as file: # add Bonds
            print('',file=file)
            print('Bonds\n',file=file)
            start_atom_number = count_1 + 1 - (n_mol - 1)*3
            count_3 = 1
            for x in np.arange(start_atom_number, count_1, 3):
                print(f"{count_3} 1 {x} {x+1} # h2o_{round(float(count_3/2)+0.1)}",file=file)
                count_3 += 1
                print(f"{count_3} 2 {x} {x+2} # h2o_{round(float(count_3/2)+0.1)}",file=file)
                count_3 += 1

            # add angles
            print('',file=file)
            print('Angles\n',file=file)
            count_4 = 1
            for x in np.arange(start_atom_number, count_1, 3):
                print(f"{count_4} 1 {x+1} {x} {x+2} # h2o_{count_4}",file=file)
                count_4 += 1

        with open(output, 'r') as file:
            lines = file.readlines()

        for i, line in enumerate(lines):
            if 'atoms' in line:
                lines.insert(i + 1, f'{count_3-1} bonds\n{count_4-1} angles\n')
            if 'atom types' in line:
                lines.insert(i + 1, '2 bond types\n1 angle types\n')

        with open(output, 'w') as file:
            file.writelines(lines)

#type_to_element = {
#    1: 'C',
#    2: 'F',
#    3: 'H',
#    4: 'N',
#    5: 'Nb',
#    6: 'Ni',
#    7: 'O',
#}
#
#N_mof_O = 32
#remove_H2O = True #False
#
#raspa_to_lmp('mof_raspa.data','conf.lmp', N_mof_O, type_to_element, remove_H2O)

