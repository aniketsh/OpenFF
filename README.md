# OpenFF

openff_topology.py generates topology of a small molecule using Open Forcefield
Usage python openff_topology.py -i input.sdf -hmr True -mol_name MOL

-hmr : Hydrogen Mass Partition for faster time step (4fs)
-mol_name : Molecule name in topology file
-i : input molecule file with explicit hydrogens 

environment.yaml file contains all packages required to run this code
To generate the environment : conda env create -f environment.yml
