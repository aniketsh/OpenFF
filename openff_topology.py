# All imports
import warnings, os, subprocess
warnings.filterwarnings('ignore')

# Openbabel is required to generate temp.pdb and temp.sdf file from user input file as
# openff requires this two as inputs
from openbabel import openbabel

# RDKIT is required only for the formal charge calculation and partial charge correction
# As sometimes the Ambertools antechamber produces non-integer formal charge
from rdkit import Chem

# OpenFF imports
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.topology import Molecule
from openmm.app import PDBFile

# PARMED : required to convert topology file formats
import parmed as pmd

# Load forcefield
force_field = ForceField("openff_unconstrained-2.0.0.offxml")

# Process user inputs
import argparse


####################################################################################################

def parameterize(input_molecule,hmr=False,input_molecule_name="MOL"):

    print("------------------------------------------------------")
    print("processing : %s"%input_molecule)
    print("setting HMR to : %s"%hmr)
    print("Molecule name will be set to : %s"% input_molecule_name)
    print("------------------------------------------------------")

    # Check the molecule extension prepare the workdir 
    current_path=os.getcwd()
    dir_name=os.path.split(input_molecule)[-1]
    dir_name=dir_name.replace(".mol2","")
    dir_name=dir_name.replace(".sdf","")
    if os.path.isdir(os.path.join(current_path,dir_name)):
        subprocess.getoutput("rm -fr %s"%(os.path.join(current_path,dir_name)))
    os.mkdir(os.path.join(current_path,dir_name))
    
    # With openbabel convert the molecule into sdf and pdb

    mol = openbabel.OBMol()

    obConversion1 = openbabel.OBConversion()
    obConversion2 = openbabel.OBConversion()

    if ".sdf" in input_molecule:
        obConversion1.SetInAndOutFormats("sdf", "pdb")
        obConversion2.SetInAndOutFormats("sdf", "sdf")
    elif ".mol2" in input_molecule:
        obConversion1.SetInAndOutFormats("mol2", "pdb")
        obConversion2.SetInAndOutFormats("sdf", "sdf")
    else:
        print("Provide valid input molecule .mol2 or .sdf with a single molecule entry")

    temp_pdb=os.path.join(current_path,dir_name,"temp.pdb")
    temp_sdf=os.path.join(current_path,dir_name,"temp.sdf")
        
    obConversion1.ReadFile(mol, input_molecule)
    obConversion1.WriteFile(mol, temp_pdb)

    obConversion2.ReadFile(mol, input_molecule)
    obConversion2.WriteFile(mol, temp_sdf)

    # Generate the topology
    try:
        ligand_off_molecule = Molecule(temp_sdf)
        ligand_pdbfile = PDBFile(temp_pdb)
        ligand_system = force_field.create_openmm_system(ligand_off_molecule.to_topology())
        ligand_structure = pmd.openmm.load_topology(ligand_pdbfile.topology,ligand_system,xyz=ligand_pdbfile.positions)
    except:
        print("Can  not parameterize, check input")

    # Make Hydrogens Heavy for 4fs timestep
    if hmr:
        pmd.tools.HMassRepartition(ligand_structure,3).execute()

    # Change the residue name, dafault MOL
    for atom in ligand_structure.atoms:
        atom.residue.name=input_molecule_name
    
    # Get formal charge
    if ".sdf" in input_molecule:
        rdkit_mol=Chem.MolFromMolFile(input_molecule)
    elif ".mol2" in input_molecule:
        rdkit_mol=Chem.MolFromMol2File(input_molecule)
    else:
        print("Provide a valid .sdf or .mol2 molecule with single molecule entry")
    formal_charge=Chem.GetFormalCharge(rdkit_mol)
        
    # Round up the formal charge and check sum of charges from the generated openff topology
    sum_partial_charge=0
    for atom in ligand_structure:
        atom.charge=round(atom.charge,3)
        sum_partial_charge+=atom.charge
        
    if not sum_partial_charge == formal_charge:
        print(sum_partial_charge,formal_charge)
        diff=formal_charge-sum_partial_charge
        ligand_structure[0].charge=ligand_structure[0].charge-diff
    else:
        print("No charge correction needed")
        
    # Write the output topologies
    try:
        os.mkdir(os.path.join(current_path,dir_name,"for_amber"))
        ligand_structure.save(os.path.join(current_path,dir_name,"for_amber","MOL.prmtop"))
        ligand_structure.save(os.path.join(current_path,dir_name,"for_amber","MOL.inpcrd"))
        ligand_structure.save(os.path.join(current_path,dir_name,"for_amber","MOL.pdb"))
    except:
        pass
    
    try:
        os.mkdir(os.path.join(current_path,dir_name,"for_gromacs"))
        ligand_structure.save(os.path.join(current_path,dir_name,"for_gromacs","MOL.top"))
        ligand_structure.save(os.path.join(current_path,dir_name,"for_gromacs","MOL.pdb"))
        ligand_structure.save(os.path.join(current_path,dir_name,"for_gromacs","MOL.gro"))
    except:
        pass


#prepare_workdir(input_molecule,hmr,input_molecule_name)
parser = argparse.ArgumentParser(
    description='OpenFF topology generator for small molecules',
)


parser.add_argument('-i',type=str,default='arg_default',nargs='?',help='input_molecule_file',required=True)
parser.add_argument('-mol_name',type=str,default='MOL',help='Molecule Name')
parser.add_argument('-hmr',type=bool,default=False,help='Makes hydrogens heavy')

args = parser.parse_args()

if args.i:
    if ".sdf" in args.i or ".mol2" in args.i:
        input_molecule=args.i
        input_molecule_name=args.mol_name
        print(args.hmr)
        parameterize(input_molecule,hmr=args.hmr,input_molecule_name=input_molecule_name)       

