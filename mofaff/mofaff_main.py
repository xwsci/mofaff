import os
import subprocess
import numpy as np
#import pacmof
from . import module_add_elements_outdump
from . import module_attach_charge_to_cif
from . import module_str_lmp_cif
from . import module_str_raspa_lmp
from . import module_cif_restartfile
from . import module_adjust_h2o

def main(input_path='input'):
    # Read parameters from input
    settings = {}
    type_to_element = {}
    reading_dict = False
    dict_name = ''
    with open(input_path, 'r') as file:  # Make sure to provide the correct path to 'input'
        for line in file:
            # Strip comments
            line = line.split('#', 1)[0].strip()
            if line.endswith('{'):
                # Start reading a dictionary
                reading_dict = True
                dict_name = line.split('=')[0].strip()
                continue
            elif reading_dict:
                if line == '}':
                    # Stop reading the dictionary
                    reading_dict = False
                    continue
                # Read dictionary entries
                key_value = line.split(':')
                key = int(key_value[0].strip())  # Assuming keys are integers as per the input example
                value = key_value[1].strip().strip(',').strip("'")  # Strip extra commas and single quotes
                type_to_element[key] = value
            elif "=" in line:
                key, value = line.split('=', 1)
                value = value.strip()
                # Try to convert to int, float, boolean, or use as string
                try:
                    settings[key.strip()] = int(value)
                except ValueError:
                    try:
                        settings[key.strip()] = float(value)
                    except ValueError:
                        if value.lower() == 'true':
                            settings[key.strip()] = True
                        elif value.lower() == 'false':
                            settings[key.strip()] = False
                        else:
                            settings[key.strip()] = value

    pressure = settings['pressure']
    N_cycle = settings['N_cycle']
    MLP_path = settings['MLP_path']
    N_mof_O = settings['N_mof_O']
    start_cycle = settings['start_cycle']
    charge_path = settings['charge_path']
    remove_H2O = settings['remove_H2O']
    adjust_H2O = settings['adjust_H2O']
    hybrid_FF = settings['hybrid_FF']

    # Read commands_LAMMPS and commands_gRASPA from input
    commands_LAMMPS = ""
    commands_gRASPA = ""
    with open(input_path, 'r') as file:
        start_collecting_LAMMPS = False
        start_collecting_gRASPA = False
        for line in file:
            line = line.strip()
            if line.startswith('commands_LAMMPS'):
                start_collecting_LAMMPS = True
                start_collecting_gRASPA = False
                continue
            elif line.startswith('commands_gRASPA'):
                start_collecting_gRASPA = True
                start_collecting_LAMMPS = False
                continue
            elif line.strip() == "":
                start_collecting_LAMMPS = False
                start_collecting_gRASPA = False
                continue
            if start_collecting_LAMMPS:
                if not line.startswith('"""'):
                    commands_LAMMPS += line + "\n"
            elif start_collecting_gRASPA:
                if not line.startswith('"""'):
                    commands_gRASPA += line + "\n"

    ## Prepare for the workflow
    direct = os.getcwd()
    subprocess.run(['cp', 'conf.lmp', 'conf.lmp.bk']) # Backup initial conf.lmp
# restart setting
    if start_cycle > 1:
        os.chdir(f"{direct}/cycle_{start_cycle-1}/2.gcmc/Movies/System_0/")
        last_snapshot = max((f for f in os.listdir('.') if f.startswith("result_") and f.endswith(".data")), key=lambda x: int(x[7:-5]), default="No files found")
        subprocess.run(['cp', last_snapshot, os.path.join(direct, 'mof_raspa.data')])
        os.chdir(direct)
        subprocess.run(['sed', '-i', r's/\(TIP4P \|mof_md_charged\.cif \)//g', 'mof_raspa.data'])
        module_str_raspa_lmp.raspa_to_lmp('mof_raspa.data','conf.lmp', N_mof_O, type_to_element, hybrid_FF)

    # Modify in.lammps file for pressure (Change if this is a NPT calculation)
#    with open('in.lammps', 'r') as file:
#        content = file.read()
#    content = content.replace("p_start", f"{pressure/100000:.6f}")
#    content = content.replace("p_end", f"{pressure/100000:.6f}")
#    with open('in.lammps', 'w') as file:
#        file.write(content)

    # Modify gcmc/simulation.input
    with open(f'{direct}/gcmc/simulation.input', 'r') as file:
        lines = file.readlines()
    with open(f'{direct}/gcmc/simulation.input', 'w') as file:
        for line in lines:
            if remove_H2O == False and line.strip().startswith('RestartFile'):
                line = line.replace(line.strip(), "RestartFile yes")
            if line.startswith('Pressure'):
                file.write(f'Pressure {pressure}\n')
            else:
                file.write(line)

    # Start the cycle of simulations
    for n in np.arange(start_cycle,N_cycle+1):
        print(f"Cycle: {n}")
        # Setup cycle_{n} directories
        cycle_dir = os.path.join(direct, f'cycle_{n}')
        os.makedirs(cycle_dir, exist_ok=True)
        os.chdir(cycle_dir)
    
        # LAMMPS simulations
        os.makedirs(os.path.join(cycle_dir, '1.lmp'), exist_ok=True)
        os.chdir(os.path.join(direct, f'cycle_{n}', '1.lmp'))
        subprocess.run(['cp', os.path.join(direct, 'conf.lmp'), os.path.join(direct, 'in.lammps'), os.path.join(direct, f'cycle_{n}', '1.lmp')])
        subprocess.run(['ln', '-s', MLP_path])
        #subprocess.run(commands_LAMMPS, shell=True, executable='/bin/bash', cwd=os.path.join(direct, f'cycle_{n}', '1.lmp'))
        subprocess.run(commands_LAMMPS, shell=True, executable='/bin/bash')
    
        module_add_elements_outdump.out_dump_add_element('out.dump', 'out_elements.dump', type_to_element)
        module_str_lmp_cif.convert_dump_cif('out_elements.dump','mof_H2O_md.cif')
    
        # gRASPA GCMC calculations
        os.chdir(os.path.join(direct, f'cycle_{n}'))
        os.makedirs('2.gcmc', exist_ok=True)
        os.chdir(os.path.join(direct, f'cycle_{n}', '2.gcmc'))
        subprocess.run(['cp', '../../gcmc/force_field.def', '.'])
        subprocess.run(['cp', '../../gcmc/force_field_mixing_rules.def', '.'])
        subprocess.run(['cp', '../../gcmc/pseudo_atoms.def', '.'])
        subprocess.run(['cp', '../../gcmc/simulation.input', '.'])
        subprocess.run(['cp', '../../gcmc/TIP4P.def', '.'])
        subprocess.run(['cp', '../1.lmp/mof_H2O_md.cif', '.'])
        module_str_lmp_cif.remove_water('mof_H2O_md.cif','mof_md.cif',N_mof_O) # output mof_md.cif with H2O removed
        # Add charge into mof_md.cif
        # data = pacmof.get_charges_single_serial('mof_md.cif', create_cif=True) # run PACMOF for each cycle, output mof_md_charged.cif
        module_attach_charge_to_cif.attach_charge('mof_md.cif', charge_path) # output mof_md_charged.cif

        # generate restartfile for H2O if remove_H2O == False
        if remove_H2O == False:
            module_adjust_h2o.separate_H2O('mof_H2O_md.cif','H2O_md.cif', N_mof_O)
            if adjust_H2O == "TIP4P":
                module_adjust_h2o.rearrange_H2O('H2O_md.cif','H2O_adjust.cif')
                module_adjust_h2o.adjust_H2O_TIP4P('H2O_adjust.cif')
                module_cif_restartfile.convert_H2O_raspa_restartfile('H2O_adjust.cif','restartfile_ini_H2O')
            else:
                module_cif_restartfile.convert_H2O_raspa_restartfile('H2O_md.cif','restartfile_ini_H2O')
            subprocess.run('mkdir -p RestartInitial/System_0', shell=True, executable='/bin/bash')
            subprocess.run(['cp', './restartfile_ini_H2O', './RestartInitial/System_0/restartfile'])

        subprocess.run(commands_gRASPA, shell=True, executable='/bin/bash')
        # Output pressure and H2O loading
        subprocess.run('cat Output/System_0_*.data >> result', shell=True, executable='/bin/bash')
        H2O_loading = float(subprocess.getoutput("grep -A 20 'LOADING: mol/kg' result | grep 'Overall' | tail -1 | awk '{print $3}' | sed 's/,//g'"))
        with open(os.path.join(direct, 'workflow.out'), 'a') as f:
            f.write(f"{n} {H2O_loading}\n")
    
        # Copy the last snapshot of GCMC to direct/conf.lmp
        os.chdir(f"{direct}/cycle_{n}/2.gcmc/Movies/System_0/")
        last_snapshot = max((f for f in os.listdir('.') if f.startswith("result_") and f.endswith(".data")), key=lambda x: int(x[7:-5]), default="No files found")
        subprocess.run(['cp', last_snapshot, os.path.join(direct, 'mof_raspa.data')])
        os.chdir(direct)
        subprocess.run(['sed', '-i', r's/\(TIP4P \|mof_md_charged\.cif \)//g', 'mof_raspa.data'])
        module_str_raspa_lmp.raspa_to_lmp('mof_raspa.data','conf.lmp', N_mof_O, type_to_element, hybrid_FF)
    
        # Check for convergence
    #    if n > 10:
    #        with open(os.path.join(direct, 'workflow.out'), 'r') as f:
    #            lines = f.readlines()
    #            if abs(float(lines[-1].split()[1])-float(lines[-2].split()[1])) < 0.01:
    #                break
    
    # After finishing all cycles
    print("Workflow complete!")
    

if __name__ == '__main__':
    import sys
    input_path = sys.argv[1] if len(sys.argv) > 1 else 'input'
    main(input_path)


