from openmm.app import PDBFile, Modeller, PDBFile, AmberPrmtopFile, AmberInpcrdFile, Simulation
from openmm.app.forcefield import ForceField
from openmm import *
from openmm.unit import *

pdb = PDBFile('ala-dipeptide.pdb')  # You can download this PDB or create it using a tool like Avogadro

modeller = Modeller(pdb.topology, pdb.positions)

# Add hydrogens and solvent if needed
modeller.addHydrogens(forcefield=ForceField('amber14-all.xml', 'amber14/tip3pfb.xml'))
modeller.addSolvent(forcefield=ForceField('amber14-all.xml', 'amber14/tip3pfb.xml'), padding=1.0*nanometer)

# Create system
system = forcefield.createSystem(modeller.topology, 
                                 nonbondedMethod=PME, 
                                 nonbondedCutoff=1*nanometer, 
                                 constraints=HBonds)

# Setup simulation
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# Write output files
PDBFile.writeFile(modeller.topology, modeller.positions, open("alanine_dipeptide.pdb", "w"))

prmtop = AmberPrmtopFile.createFromSystem(system, modeller.topology)
inpcrd = AmberInpcrdFile.writeFile(modeller.positions, velocities=None, boxVectors=None)

prmtop.writeFile("test.prmtop")
inpcrd.writeFile("test.inpcrd")

print("Amber topology (test.prmtop) and coordinate file (test.inpcrd) successfully created!")
