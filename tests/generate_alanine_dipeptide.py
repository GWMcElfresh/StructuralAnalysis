from openmm.app import *
from openmm import *
from openmm.unit import *

# Create alanine dipeptide (ACE-ALA-NME)
residues = ['ACE', 'ALA', 'NME']
topology = Topology()
chain = topology.addChain()

for res in residues:
    topology.addResidue(res, chain)

# Load force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Generate positions with Modeller
modeller = Modeller(topology, [])
modeller.addHydrogens(forcefield)  # Adds missing hydrogens

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
