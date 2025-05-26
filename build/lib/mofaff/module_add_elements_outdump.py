def out_dump_add_element(input_file, output_file, mapping):
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('ITEM: ATOMS'):
                fout.write(line.strip() + ' element\n')
            elif line.startswith('ITEM:'):
                fout.write(line)
            else:
                data = line.split()
                if len(data) >= 5:
                    atom_type = int(data[2])
                    element = mapping.get(atom_type, 'Unknown')
                    fout.write(" ".join(data) + f" {element}\n")
                else:
                    fout.write(line)


