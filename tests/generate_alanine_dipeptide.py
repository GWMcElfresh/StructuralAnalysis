from openmm.app import ForceField, PDBFile, Modeller, AmberPrmtopFile, AmberInpcrdFile, Simulation, PME, HBonds
from openmm import *
from openmm.unit import *

pdb = PDBFile('ala-dipeptide.pdb')  # You can download this PDB or create it using a tool like Avogadro

modeller = Modeller(pdb.topology, pdb.positions)

# Load the force field
ff = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# Add hydrogens and solvent if needed
modeller.addHydrogens(forcefield=ff)
modeller.addSolvent(forcefield=ff, padding=1.0 * nanometer)

# Create system
system = ff.createSystem(modeller.topology, 
                         nonbondedMethod=PME, 
                         nonbondedCutoff=1 * nanometer, 
                         constraints=HBonds)

# Setup simulation
integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# Write output files
with open("alanine_dipeptide.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
  
print("PDB file successfully created!")
