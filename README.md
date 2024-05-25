# MOFAFF

A package to calculate MOF adsorption with framework flexibility (MOFAFF). This package iteratively runs machine learning potential (MLP) based MD (LAMMPS) and MC (gRASPA) simulations to calculate gas uptake within porous materials, considering framework flexibility.

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


