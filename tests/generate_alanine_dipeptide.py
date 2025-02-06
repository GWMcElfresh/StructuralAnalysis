from openmm.app import *
from openmm import *
from openmm.unit import *

# Build alanine dipeptide with explicit capping groups
pdb_str = """\
ATOM      1  C   ACE     1      -0.677  -0.395   0.000  1.00  0.00           C  
ATOM      2  O   ACE     1       0.677  -0.395   0.000  1.00  0.00           O  
ATOM      3  CH3 ACE     1       1.210   0.935   0.000  1.00  0.00           C  
ATOM      4  N   ALA     2       2.210   1.210   0.000  1.00  0.00           N  
ATOM      5  CA  ALA     2       3.210   2.540   0.000  1.00  0.00           C  
ATOM      6  C   ALA     2       4.677   2.815   0.000  1.00  0.00           C  
ATOM      7  O   ALA     2       5.210   3.815   0.000  1.00  0.00           O  
ATOM      8  CB  ALA     2       3.210   1.210   1.200  1.00  0.00           C  
ATOM      9  N   ALA     3       5.677   2.210   0.000  1.00  0.00           N  
ATOM     10  CA  ALA     3       6.210   3.540   0.000  1.00  0.00           C  
ATOM     11  C   ALA     3       7.677   3.815   0.000  1.00  0.00           C  
ATOM     12  O   ALA     3       8.210   4.815   0.000  1.00  0.00           O  
ATOM     13  CB  ALA     3       6.210   2.210   1.200  1.00  0.00           C  
ATOM     14  C   NME     4       8.677   2.815   0.000  1.00  0.00           C  
ATOM     15  H   NME     4       9.210   3.815   0.000  1.00  0.00           H  
TER
END
"""

# Save PDB file
with open("alanine_dipeptide.pdb", "w") as f:
    f.write(pdb_str)

# Load the PDB into OpenMM
pdb = PDBFile("alanine_dipeptide.pdb")
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Use Modeller to add missing hydrogens
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)  # No pH adjustment

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
