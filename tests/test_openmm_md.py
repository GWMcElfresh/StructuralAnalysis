import openmm as mm
from openmm.app import *
from openmm import unit
import pytest

# Test function to validate OpenMM setup
def test_openmm_md():
    #minimal system
    prmtop = AmberPrmtopFile('test.prmtop')
    inpcrd = AmberInpcrdFile('test.inpcrd')

    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*unit.nanometer, constraints=HBonds)
    
    integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)

    platform = mm.Platform.getPlatformByName('CPU')
    simulation = Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)

    #test: minimize energy and run
    simulation.minimizeEnergy()
    simulation.step(100)

    # Passes if no exceptions occur
    assert True
