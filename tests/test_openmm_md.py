import pytest
import os
from openmm.app import *
from openmm import *
from openmm.unit import *

@pytest.fixture
def run_md_simulation():
    """Runs a short OpenMM MD simulation and checks for successful execution."""
    pdb_file = 'alanine_dipeptide.pdb'
    output_dcd = 'output.dcd'
    output_log = 'output.log'

    # Ensure PDB file exists
    assert os.path.exists(pdb_file), f"Missing PDB file: {pdb_file}"

    # Load PDB
    pdb = PDBFile(pdb_file)

    # Load force field
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # Create system
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1 * nanometer,
        constraints=HBonds
    )

    # Set up integrator
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

    # Set up simulation
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)

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

    print("MD simulation ran successfully!")

# Run the test
@pytest.mark.parametrize("test_name", ["OpenMM MD Simulation"])
def test_md_run(run_md_simulation):
    """Test if the MD simulation completes and generates expected files."""
    pass  # The fixture runs automatically and asserts success
