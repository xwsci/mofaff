# MOFAFF

A package to calculate MOF adsorption with framework flexibility (MOFAFF). This package iteratively runs machine learning potential (MLP) based MD (LAMMPS) and MC (gRASPA) simulations to calculate gas uptake within porous materials, considering framework flexibility. (Refer to workflow.png for more details)

## Requirement
1. numpy>=1.24.4; ase>=3.22.1; python>=3.8
2. ML potential trained using DeepMD-kit
3. LAMMPS with MLP installed
4. gRASPA for GCMC simulations

## Installation

pip install .

or 

pip install . --upgrade

#Future:
#pip install mofaff


## Usage

import mofaff

mofaff.mofaff_main('path_to_input')


