# MOFAFF

A package to calculate MOF Adsorption with Framework Flexibility (MOFAFF)

## Requirement
1. ML potential trained using DeepMD-kit
2. LAMMPS with MLP installed
3. gRASPA for GCMC simulations

## Installation

pip install .

or 

pip install . --upgrade

#Future:
#pip install mofaff


## Usage

import mofaff
mofaff.mofaff_main('path_to_input')


