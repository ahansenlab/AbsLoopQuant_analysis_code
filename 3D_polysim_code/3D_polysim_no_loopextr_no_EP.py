# The following code was adapted from the Open2C polychrom package's example code (https://github.com/open2c/polychrom/blob/master/examples/example/example.py).
# In this condition:
# [No] Loop Extrusion
# [No] Enhancer-promoter attraction

import os, sys
import polychrom
from polychrom import simulation, starting_conformations, forces, forcekits
import openmm
import os
from polychrom.hdf5_format import HDF5Reporter
import numpy as np
from datetime import datetime

gpu = '0'

N = 70000
d = 0.3
l = int((N/d) ** (1/3))

num_chains = 1
chain = [(0,70000,0)]

collision_rate = 1
errtol = 0.01
truncated_potentials = 3

interaction_matrix = np.array([[0, 0],[0, 0]])
monomer_types = np.zeros(N, dtype='int')

polymer = starting_conformations.grow_cubic(N, l)

date = datetime.today().strftime('%Y%m%d')
folder_name = f'/mnt/coldstorage/jjusuf/simabsloopquant_data/region_3_simulations/sim_{date}_no_LE_no_EP_clean'

reporter = HDF5Reporter(folder=folder_name, max_data_length=100, overwrite=True)
sim = simulation.Simulation(
    platform="CUDA",
    integrator="variableLangevin",
    error_tol=errtol,
    GPU=gpu,
    collision_rate=collision_rate,
    N=N,
    save_decimals=2,
    PBCbox=False,
    reporters=[reporter],
)

sim.set_data(polymer, center=True)  # loads a polymer, puts a center of mass at zero

sim.add_force(forces.spherical_confinement(sim, density=d, k=1))

sim.add_force(
    forcekits.polymer_chains(
        sim,
        chains=chain, # makes circular chains
        bond_force_func=forces.harmonic_bonds, # adds harmonic bonds
        bond_force_kwargs={
            'bondLength':1.0,
            'bondWiggleDistance':0.1, # Bond distance will fluctuate +- 0.1 on average
            # This is fine because our simulation is "soft" on average (all forces are soft-core)
            # And this is potentially desirable to make polymer chain softer to avoid jerking 
            # when a new bond is initialized,
            'override_checks':True,
         },
        angle_force_func=forces.angle_force,
        angle_force_kwargs={
            'k':0.05, # we are making a very flexible polymer, basically not necessary here,
            'override_checks':True,
        },
        nonbonded_force_func=forces.heteropolymer_SSW,
        nonbonded_force_kwargs={
            # "trunc": 3.0,  # this will let chains cross sometimes
            #'trunc':10.0, # this will resolve chain crossings and will not let chain cross anymore
            "interactionMatrix": interaction_matrix,
            "monomerTypes": monomer_types,
            "extraHardParticlesIdxs": [],
            "attractionEnergy": 0,
            "attractionRadius": 2,
            "repulsionRadius": 1.05,  # this is from old code
            "repulsionEnergy": truncated_potentials,
        },
        except_bonds=True,
        #extra_bonds=extra_bonds,
        override_checks=True,
    )
)

for _ in range(1000000):  # Do 1000000 blocks
    sim.do_block(2360)  # Of 236*10 timesteps each. Data is saved automatically.

reporter.dump_data()  # always need to run in the end to dump the block cache to the disk