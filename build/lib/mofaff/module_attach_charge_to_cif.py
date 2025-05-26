
def attach_charge(filename,charge_path):
    # delete all lines in between the last loop_ and label:
    with open(filename, "r") as file:
        lines = file.readlines()
    start = max(i for i, line in enumerate(lines) if line.strip() == "loop_")
    end = lines.index("  _atom_site_label\n", start)
    lines = lines[:start+1] + lines[end:]
    with open(filename, "w") as file:
        file.writelines(lines)

    # delete all _atom lines after label 
    with open(filename, "r") as file:
        lines = file.readlines()
    index = lines.index("  _atom_site_label\n")
    lines = lines[:index + 1] + [line for line in lines[index + 1:] if "_atom" not in line]
    with open(filename, "w") as file:
        file.writelines(lines)

    # insert the lines we need and tune the columns accordingly
    with open(filename, 'r') as file:
        lines = file.readlines()
    with open(filename, 'w') as file:
        for line in lines:
            line = line.replace('_atom_site_label', '_atom_site_label\n  _atom_site_type_symbol\n  _atom_site_fract_x\n  _atom_site_fract_y\n  _atom_site_fract_z\n  _atom_site_charge')
            file.writelines(line)

    found_charge = False
    with open(filename, 'r') as file:
        lines = file.readlines()
    with open(filename, 'w') as file:
        for line in lines:
            if '_atom_site_charge' in line:
                found_charge = True
            if found_charge and len(line.split()) > 3:
                parts = line.split()
                line = f"{parts[1]} {parts[0]} {parts[3]} {parts[4]} {parts[5]}\n"
            file.write(line)

    # add charge to cif
    with open(filename, 'r') as f1, open(charge_path, 'r') as f2:
        file1_lines = f1.readlines()
        file2_lines = f2.readlines()

    output_lines = []
    add_charge = False
    charge_index = 0

    for line in file1_lines:
        line = line.strip()
        if line:  # Only process non-empty lines
            if add_charge and charge_index < len(file2_lines):
                line += ' ' + file2_lines[charge_index].strip()
                charge_index += 1
            output_lines.append(line)
            if "_atom_site_charge" in line:
                add_charge = True

    with open(filename.split('.cif')[0]+'_charged.cif', 'w') as f_out:
        for line in output_lines:
            f_out.write(line + '\n')



