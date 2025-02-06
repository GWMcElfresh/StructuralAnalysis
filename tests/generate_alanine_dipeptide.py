from openmm.app import *
from openmm import *
from openmm.unit import *
import os

# Build alanine dipeptide molecule manually with terminal groups
pdb_str = """\
ATOM      1  N   NALA    1      -0.677  -0.395   0.000  1.00  0.00           N  
ATOM      2  CA  NALA    1       0.677  -0.395   0.000  1.00  0.00           C  
ATOM      3  C   NALA    1       1.210   0.935   0.000  1.00  0.00           C  
ATOM      4  O   NALA    1       0.438   1.897   0.000  1.00  0.00           O  
ATOM      5  CB  NALA    1       1.210  -1.312   1.200  1.00  0.00           C  
ATOM      6  HN  NALA    1      -1.210  -1.312  -1.200  1.00  0.00           H  
ATOM      7  HA  NALA    1       1.210  -1.312  -1.200  1.00  0.00           H  
ATOM      8  HB1 NALA    1       2.354  -1.312   1.200  1.00  0.00           H  
ATOM      9  HB2 NALA    1       0.677  -2.354   1.200  1.00  0.00           H  
ATOM     10  HB3 NALA    1       1.210  -1.312   2.354  1.00  0.00           H  
ATOM     11  N   CALA    2       2.677   1.210   0.000  1.00  0.00           N  
ATOM     12  CA  CALA    2       3.210   2.540   0.000  1.00  0.00           C  
ATOM     13  C   CALA    2       4.677   2.815   0.000  1.00  0.00           C  
ATOM     14  O   CALA    2       5.210   3.815   0.000  1.00  0.00           O  
TER
END
"""

# Save the PDB file
with open("alanine.pdb", "w") as f:
    f.write(pdb_str)

# Load the PDB into OpenMM
pdb = PDBFile("alanine.pdb")

# Load force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Use Modeller to add missing hydrogens and terminal groups
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)

# Create system
system = forcefield.createSystem(modeller.topology, 
                                 nonbondedMethod=PME, 
                                 nonbondedCutoff=1*nanometer, 
                                 constraints=HBonds)

# Create simulation
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(modeller.topology, system, integrator)

# Assign initial positions
simulation.context.setPositions(modeller.positions)

# Save Amber topology and coordinate files
prmtop = AmberPrmtopFile.createFromSystem(system, modeller.topology)
inpcrd = AmberInpcrdFile.writeFile(modeller.positions, velocities=None, boxVectors=None)

# Save files
prmtop.writeFile("./test.prmtop")
inpcrd.writeFile("./test.inpcrd")

print("Amber topology (test.prmtop) and coordinate file (test.inpcrd) successfully created!")
