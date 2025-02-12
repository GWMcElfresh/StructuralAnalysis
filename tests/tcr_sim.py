import os
from openmm.app import *
from openmm import *
from openmm.unit import *
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt


def run_md_simulation():
    """Runs a short OpenMM MD simulation and checks for successful execution."""
    pdb_file = 'AF-Q30677-F1-model_v4.pdb'
    output_dcd = 'output.dcd'
    output_log = 'output.log'

    # Ensure PDB file exists
    assert os.path.exists(pdb_file), f"Missing PDB file: {pdb_file}"

    # Load PDB
    pdb = PDBFile(pdb_file)

    # Load force field
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield=ForceField('amber14-all.xml', 'amber14/tip3pfb.xml'))
    modeller.addSolvent(forcefield=ForceField('amber14-all.xml', 'amber14/tip3pfb.xml'), padding=1.0*nanometer)
    # Create system
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1 * nanometer,
        constraints=HBonds
    )

    # Set up integrator
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

    # Set up simulation
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    # Minimize energy
    simulation.minimizeEnergy()

    # Assign velocities
    simulation.context.setVelocitiesToTemperature(300*kelvin)

    # Set up reporters (write trajectory and log output)
    simulation.reporters.append(DCDReporter(output_dcd, 1000))  
    simulation.reporters.append(StateDataReporter(output_log, 1000, step=True, 
                                                  potentialEnergy=True, temperature=True))

    # Run a short MD simulation (e.g., 10,000 steps for testing)
    simulation.step(10000)

    # Ensure output files exist
    assert os.path.exists(output_dcd), "Trajectory file (DCD) was not generated!"
    assert os.path.exists(output_log), "Log file was not generated!"

    # Check log file has content
    with open(output_log, 'r') as log:
        log_content = log.read()
        assert "Step" in log_content, "Simulation log does not contain expected output!"

    del simulation
    PDBFile.writeFile(modeller.topology, modeller.positions, open("topology.pdb", "w"))
    print("MD simulation ran successfully!")



    # Load the trajectory using the generated PDB file as topology
    u = mda.Universe("topology.pdb", "output.dcd")  

    # Select backbone Cα atoms
    protein = u.select_atoms("protein and name CA")
    n_residues = len(protein.residues)

    # Define contact threshold (e.g., 5.0 Å)
    contact_threshold = 5.0  

    # Initialize contact matrix
    contact_map = np.zeros((n_residues, n_residues))

    # Compute contacts over trajectory
    for ts in u.trajectory:
        distances = np.linalg.norm(protein.positions[:, None, :] - protein.positions[None, :, :], axis=-1)
        contact_map += (distances < contact_threshold).astype(int)  # Binary contact map

    # Normalize by number of frames
    contact_map /= len(u.trajectory)

    # Plot contact map
    plt.figure(figsize=(8, 6))
    plt.imshow(contact_map, cmap="hot", origin="lower")
    plt.colorbar(label="Fraction of frames in contact")
    plt.xlabel("Residue Index")
    plt.ylabel("Residue Index")
    plt.title("Residue Contact Map")
    plt.show()


run_md_simulation()
