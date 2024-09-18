# The following code was adapted from the simulation code used in JH Yang, HB Brandao, and AS Hansen, DNA double-strand break end synapsis by DNA loop extrusion. Nature Communications 14, 1913 (2023).
# In this condition:
# [Yes] Loop Extrusion
# [Yes] Enhancer-promoter attraction (3 kBT attractive energy)

from __future__ import absolute_import, division, print_function
import numpy as np
import sys
import os
import time
import tempfile
import logging
import warnings

import pickle
import os
import time
import numpy as np
import polychrom

from polychrom import polymerutils
from polychrom import forces
from polychrom import forcekits
from polychrom.simulation import Simulation
from polychrom.starting_conformations import grow_cubic
from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file

import os 
import shutil

import pyximport; pyximport.install()
pyximport.install(setup_args={'include_dirs': np.get_include()})
from DSB_smcTranslocator_v2 import smcTranslocatorDirectional 

import warnings
import h5py 
import glob
import re

from itertools import product
from scipy.ndimage.filters import gaussian_filter
from scipy.sparse import coo_matrix
from scipy.ndimage.filters import gaussian_filter1d
from scipy.stats import expon

from datetime import datetime

from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file
import h5py as hp
from matplotlib import pyplot as plt

from collections.abc import Iterable

import openmm
import simtk.unit


logging.basicConfig(level=logging.INFO)

class DSB_Simulation(Simulation):
    def do_block(
        self,
        steps=None,
        check_functions=[],
        get_velocities=False,
        save=True,
        save_extras={},
        positions_to_sample=None,
    ):
        """performs one block of simulations, doing steps timesteps,
        or steps_per_block if not specified.

        Parameters
        ----------

        steps : int or None
            Number of timesteps to perform.
        increment : bool, optional
            If true, will not increment self.block and self.steps counters
        """

        if not self.forces_applied:
            if self.verbose:
                logging.info("applying forces")
                sys.stdout.flush()
            self._apply_forces()
            self.forces_applied = True

        a = time.time()
        self.integrator.step(steps)  # integrate!

        self.state = self.context.getState(
            getPositions=True, getVelocities=get_velocities, getEnergy=True
        )

        b = time.time()
        coords = self.state.getPositions(asNumpy=True)
        newcoords = coords / simtk.unit.nanometer
        newcoords = np.array(newcoords, dtype=np.float32)
        if self.kwargs["save_decimals"] is not False:
            newcoords = np.round(newcoords, self.kwargs["save_decimals"])

        self.time = self.state.getTime() / simtk.unit.picosecond

        # calculate energies in KT/particle
        eK = self.state.getKineticEnergy() / self.N / self.kT
        eP = self.state.getPotentialEnergy() / self.N / self.kT
        curtime = self.state.getTime() / simtk.unit.picosecond

        msg = "block %4s " % int(self.block)
        msg += "pos[1]=[%.1lf %.1lf %.1lf] " % tuple(newcoords[0])

        check_fail = False
        for check_function in check_functions:
            if not check_function(newcoords):
                check_fail = True

        if np.isnan(newcoords).any():
            raise IntegrationFailError("Coordinates are NANs")
        if eK > self.eK_critical:
            print(eK)
            print(self.eK_critical)
            raise EKExceedsError("Ek={1} exceeds {0}".format(self.eK_critical, eK))
        if (np.isnan(eK)) or (np.isnan(eP)):
            raise IntegrationFailError("Energy is NAN)")
        if check_fail:
            raise IntegrationFailError("Custom checks failed")

        dif = np.sqrt(np.mean(np.sum((newcoords - self.get_data()) ** 2, axis=1)))
        msg += "dr=%.2lf " % (dif,)
        self.data = coords
        msg += "t=%2.1lfps " % (self.state.getTime() / simtk.unit.picosecond)
        msg += "kin=%.2lf pot=%.2lf " % (eK, eP)
        msg += "Rg=%.3lf " % self.RG()
        msg += "SPS=%.0lf " % (steps / (float(b - a)))

        if (
            self.integrator_type.lower() == "variablelangevin"
            or self.integrator_type.lower() == "variableverlet"
        ):
            dt = self.integrator.getStepSize()
            msg += "dt=%.1lffs " % (dt / simtk.unit.femtosecond)
            mass = self.system.getParticleMass(1)
            dx = simtk.unit.sqrt(2.0 * eK * self.kT / mass) * dt
            msg += "dx=%.2lfpm " % (dx / simtk.unit.nanometer * 1000.0)

        logging.info(msg)

        if positions_to_sample is None:
            result = {
                "pos": newcoords,
                "potentialEnergy": eP,
                "kineticEnergy": eK,
                "time": curtime,
                "block": self.block,
            }
        else: 
            
            result = {
                "pos": [newcoords[pos] for pos in positions_to_sample],
                "potentialEnergy": eP,
                "kineticEnergy": eK,
                "time": curtime,
                "block": self.block,
            }
    
        if get_velocities:
            result["vel"] = self.state.getVelocities() / (
                simtk.unit.nanometer / simtk.unit.picosecond
            )
        result.update(save_extras)
        if save:
            for reporter in self.reporters:
                reporter.report("data", result)

        self.block += 1
        self.step += steps

        return result    
    
def run_simulation(idx, N_monomers, translocator_initialization_steps, \
                  smcStepsPerBlock, \
                  smcBondDist, smcBondWiggleDist,\
                  steps_per_block, volume_density, block, \
                  chain, extra_bonds, save_folder, save_every_x_blocks, \
                  total_saved_blocks, restartBondUpdaterEveryBlocks, GPU_choice = 0, 
                   overwrite=False,density=0.3,positions_to_sample=None,monomer_types=None,interaction_matrix=None, colrate=0.3,errtol=0.01,trunc=3,
                  initial_conformation=None,block_to_save_all = None,save_length=10000):
    
    print(f"Starting simulation. Doing {total_saved_blocks} steps, saving every {save_every_x_blocks} block(s).")    
    print(f"A total of {total_saved_blocks} steps will be performed")
    
    # assertions for easy managing code below 
    assert restartBondUpdaterEveryBlocks % save_every_x_blocks == 0 
    assert (total_saved_blocks * save_every_x_blocks) % restartBondUpdaterEveryBlocks == 0 
    savesPerBondUpdater = restartBondUpdaterEveryBlocks // save_every_x_blocks
    BondUpdaterInitsTotal  = (total_saved_blocks ) * save_every_x_blocks // restartBondUpdaterEveryBlocks
    print("BondUpdater will be initialized {0} times".format(BondUpdaterInitsTotal))
    
    
    # create SMC translocation object   
    SMCTran = initModel(idx)                        
    SMCTran.steps(translocator_initialization_steps)  # steps to "equilibrate" SMC dynamics

    # now feed bond generators to BondUpdater 
    BondUpdater = simulationBondUpdater(SMCTran)   
       
    # clean up the simulation directory
    folder = save_folder
    
    # set initial configuration of the polymer    
    if initial_conformation is None:
        box = (N_monomers / volume_density) ** 0.33  
        data = grow_cubic(N_monomers, int(box) - 2, method='linear')  # creates a compact conformation, but one can take pre-existing configurations
    else:
        data = initial_conformation
        
    # create the reporter class - this saves polymer configurations and other information you specify
    reporter = HDF5Reporter(folder=folder, max_data_length=save_length,overwrite=True,blocks_only=True)
    
    # Iterate over various BondUpdaterInitializations
    for BondUpdaterCount in range(BondUpdaterInitsTotal):
        
        # simulation parameters are defined here     
        sim = DSB_Simulation(
                platform="cuda",
                integrator="variableLangevin", 
                error_tol=errtol, 
                GPU = "{}".format(GPU_choice), 
                collision_rate=colrate, 
                N = len(data),
                reporters=[reporter],
                PBCbox=False,
                precision="mixed")

        sim.set_data(data)  # loads a polymer, puts a center of mass at zero
        sim.add_force(forces.spherical_confinement(sim,density=density,k=1))
        # -----------Adding forces ---------------
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
                    "repulsionEnergy": trunc,
                },
                except_bonds=True,
                #extra_bonds=extra_bonds,
                override_checks=True,
            )
        )

        sim.step = block
        kbond = sim.kbondScalingFactor / (smcBondWiggleDist ** 2)
        bondDist = smcBondDist * sim.length_scale
        activeParams = {"length":bondDist,"k":kbond}
        inactiveParams = {"length":bondDist, "k":0}
        BondUpdater.setParams(activeParams, inactiveParams)

        # this step actually puts all bonds in and sets first bonds to be what they should be
        BondUpdater.setup( BondUpdaterCount,
                     bondForce=sim.force_dict['harmonic_bonds'],
                    smcStepsPerBlock=smcStepsPerBlock,
                    blocks=restartBondUpdaterEveryBlocks) 
        print("Restarting BondUpdater")
        
        # Minimize energy since we're reconfiguring the random walk into crumples and SMC-dependent loops
        sim.local_energy_minimization() 

        # Iterate over simulation time steps within each BondUpdater 
        for b in range(restartBondUpdaterEveryBlocks):
            # BondUpdater updates bonds at each time step.
            curBonds, pastBonds = BondUpdater.step(sim.context)  
            if (b % save_every_x_blocks == 0):  
                # save SMC positions and monomer positions 
                total_blocks_so_far = BondUpdaterCount*restartBondUpdaterEveryBlocks + b
                if total_blocks_so_far in block_to_save_all:
                    sim.do_block(steps=steps_per_block, save_extras={"SMCs":curBonds,"SMC_step":sim.step},
                             positions_to_sample=[i for i in range(N_monomers)])
                else:
                    sim.do_block(steps=steps_per_block, save_extras={"SMCs":curBonds,"SMC_step":sim.step},
                             positions_to_sample=positions_to_sample)
            else:
                sim.integrator.step(steps_per_block)  # do steps without getting the positions from the GPU (faster)

        data = sim.get_data()  # get the polymer positions (i.e. data) 
        block = sim.step       # get the time step
        del sim               # delete the simulation
        time.sleep(0.2)      # wait 200ms to let garbage collector do its magic
    # dump data to  output file!
    reporter.dump_data() 
    done_file = open(os.path.join(folder,'sim_done.txt'),"w+")
    done_file.close()
    
def initModel(idx):

    #unchanging parameters
    BELT_ON=0
    BELT_OFF=1
    switchRate = 0 
    SWITCH_PROB= switchRate # switching rate
    PUSH=0
    PAIRED=0
    SLIDE=1
    SLIDE_PAUSEPROB=0.99 
    loop_prefactor=1.5 
    FULL_LOOP_ENTROPY=1 
    FRACTION_ONESIDED=0

    # Extruder dynamics frequency of sampling parameters 
    numExtruderSteps = 500 # steps taken for each simulation sample
    numInitializationSteps = 10000 # how long we take to equilibrate the simulation

    # Polymer and extruder dynamics parameters #
    processivity = lof[idx][3] # processivity
    separations = lof[idx][4] # separation
    longlived_fraction = lof[idx][8] # fraction of all LEFs that are longlived
    PAUSEPROB= pause_prob # motor pause probability
    
    normal_sep = separations/(1-longlived_fraction)
    smcNum = int(chrom_size//normal_sep) # number of SMCs loaded
    if longlived_fraction==0:
        longlived_smcNum = 0
    else:
        longlived_sep = separations/longlived_fraction
        longlived_smcNum = int(chrom_size//longlived_sep)
        
    SWITCH =  np.ones(chrom_size,dtype=np.double)*SWITCH_PROB

    birthArray = np.ones(chrom_size)/chrom_size
    deathArray = np.zeros(chrom_size, dtype=np.double) + 1. / (0.5*processivity/(1-PAUSEPROB))  # times 0.5 to account for two-sided extrusion 
    longlived_deathArray = np.zeros(chrom_size, dtype=np.double) + 1. / (0.5*processivity*lof[idx][9]/(1-PAUSEPROB)) 
    deathArray[boundary_coordinates] = 1./ (0.5*processivity/(1-PAUSEPROB))/ lof[idx][7]
    longlived_deathArray[boundary_coordinates] = 1./ (0.5*processivity*lof[idx][9]/(1-PAUSEPROB)) / lof[idx][7]

    boundary_flag = (boundaryStrengthsL>0) + (boundaryStrengthsR>0)


    pauseArray = PAUSEPROB*np.ones(chrom_size, dtype=np.double) 
    
    slidePauseArray = np.zeros(chrom_size, dtype=np.double) + SLIDE_PAUSEPROB
    oneSidedArray = np.zeros(smcNum, dtype=np.int64)
    longlived_oneSidedArray = np.zeros(longlived_smcNum, dtype=np.int64)
    belt_on_array = np.zeros(smcNum, dtype=np.double) + BELT_ON
    belt_off_array = np.zeros(smcNum, dtype=np.double) + BELT_OFF
    spf=slidePauseArray*(1.-(1.-SLIDE_PAUSEPROB)*np.exp(-1.*loop_prefactor))
    spb=slidePauseArray*(1.-(1.-SLIDE_PAUSEPROB)*np.exp(loop_prefactor))
    
    ################### TAD BOUNDARY CONDITION###################
    # stallLeft is from the LEF perspective, boundaryStrength is from the CTCF perspective
    stallLeftArray = boundaryStrengthsR
    stallRightArray =  boundaryStrengthsL
    ##################################################################

    transloc = smcTranslocatorDirectional(birthArray, deathArray, longlived_deathArray, stallLeftArray, stallRightArray, pauseArray,
                                         smcNum, longlived_smcNum, oneSidedArray,longlived_oneSidedArray, FRACTION_ONESIDED, slide=SLIDE,
                                         slidepauseForward=spf, slidepauseBackward=spb, switch=SWITCH, pushing=PUSH,
                                        belt_on=belt_on_array, belt_off=belt_off_array,SLIDE_PAUSEPROB=SLIDE_PAUSEPROB) 
    return transloc


class simulationBondUpdater(object):
    """
    This class precomputes simulation bonds for faster dynamic allocation. 
    """

    def __init__(self, smcTransObject):#, plectonemeObject):
        """
        :param smcTransObject: smc translocator object to work with
        """
        self.smcObject = smcTransObject
        self.allBonds = []

    def setParams(self, activeParamDict, inactiveParamDict):
        """
        A method to set parameters for bonds.
        It is a separate method because you may want to have a Simulation object already existing

        :param activeParamDict: a dict (argument:value) of addBond arguments for active bonds
        :param inactiveParamDict:  a dict (argument:value) of addBond arguments for inactive bonds

        """
        self.activeParamDict = activeParamDict
        self.inactiveParamDict = inactiveParamDict


    def setup(self, BondUpdaterCount,  bondForce, smcStepsPerBlock, blocks = 100):
        """
        A method that milks smcTranslocator object
        and creates a set of unique bonds, etc.

        :param bondForce: a bondforce object (new after simulation restart!)
        :param blocks: number of blocks to precalculate
        :param smcStepsPerBlock: number of smcTranslocator steps per block
        :return:
        """

        if len(self.allBonds) != 0:
            raise ValueError("Not all bonds were used; {0} sets left".format(len(self.allBonds)))

        self.bondForce = bondForce

#         # update smc translocator object
#         if BondUpdaterCount == BondUpdaterCountBeforeDSB:
            
#             stall_prob_left,stall_prob_right = self.smcObject.get_stall_prob()
#             loading_prob = self.smcObject.get_emission()
#             # make double strand break site impermeable
#             stall_prob_left[DSB_coordinates[1::2]] = 1 # even elements of the DSB coordinates
#             stall_prob_right[DSB_coordinates[::2]] = 1 # odd elements of the DSB coordinates
#             # boosting by DSB
#             processivity = lof[idx][3] 
#             PAUSEPROB=pause_prob # motor pause probability
#             deathArray = np.zeros(chrom_size, dtype=np.double) + 1. / (0.5*processivity/(1-PAUSEPROB))  # times 0.5 to account for two-sided extrusion 
#             longlived_deathArray = np.zeros(chrom_size, dtype=np.double) + 1. / (0.5*processivity*lof[idx][9]/(1-PAUSEPROB)) 
#             deathArray[boundary_coordinates] = 1./ (0.5*processivity/(1-PAUSEPROB))/ lof[idx][7]
#             longlived_deathArray[boundary_coordinates] = 1./ (0.5*processivity*lof[idx][9]/(1-PAUSEPROB)) / lof[idx][7]
#             deathArray[DSB_coordinates] =  1./ (0.5*processivity/(1-PAUSEPROB))/ lof[idx][10] 
#             longlived_deathArray[DSB_coordinates] =1./ (0.5*processivity*lof[idx][9]/(1-PAUSEPROB)) / lof[idx][10]

#             # DSB flag to ensure correct valency at DSB site: DSB can only stabilize one LEF
#             DSB_flag = np.zeros(chrom_size, int)
#             DSB_flag[DSB_coordinates] = 1 
#             self.smcObject.updateStallprob(stall_prob_left,stall_prob_right)
#             self.smcObject.updateDSBFlag(DSB_flag)
#             self.smcObject.updateFalloffProb(deathArray,longlived_deathArray)
#             # implement targeted loading mechanism
#             loading_prob[DSB_coordinates]= 1/chrom_size*lof[idx][11]
#             self.smcObject.updateEmissionProb(loading_prob)
#             pauseArray = PAUSEPROB*np.ones(chrom_size, dtype=np.double) 
#             self.smcObject.updatePauseProb(pauseArray) # increase sampling rate to not miss any synapsis events
            

        #precalculating all bonds
        allBonds = []
        
        left, right = self.smcObject.getSMCs()
        # add SMC bonds
        bonds = [(int(i), int(j)) for i,j in zip(left, right)]

        left, right = self.smcObject.getlonglivedSMCs()
        # add longlivedSMC bonds
        bonds += [(int(i), int(j)) for i,j in zip(left, right)]
        self.curBonds = bonds.copy()
        
        allBonds.append(bonds)
        for dummy in range(blocks):
            self.smcObject.steps(smcStepsPerBlock)
            left, right = self.smcObject.getSMCs()
            # add SMC bonds
            bonds = [(int(i), int(j)) for i,j in zip(left, right)]
            
            left, right = self.smcObject.getlonglivedSMCs()
            # add longlivedSMC bonds
            bonds += [(int(i), int(j)) for i,j in zip(left, right)]

            allBonds.append(bonds)
        
        self.allBonds = allBonds
        self.uniqueBonds = list(set(sum(allBonds, [])))
        
        allBonds.pop(0)
        

        #adding forces and getting bond indices
        self.bondInds = []
        
    
        for bond in self.uniqueBonds:
            paramset = self.activeParamDict if (bond in self.curBonds) else self.inactiveParamDict
            ind = bondForce.addBond(bond[0], bond[1], **paramset) # changed from addBond
            self.bondInds.append(ind)
        self.bondToInd = {i:j for i,j in zip(self.uniqueBonds, self.bondInds)}
        return self.curBonds,[]


    def step(self, context, verbose=False):
        """
        Update the bonds to the next step.
        It sets bonds for you automatically!
        :param context:  context
        :return: (current bonds, previous step bonds); just for reference
        """
        if len(self.allBonds) == 0:
            raise ValueError("No bonds left to run; you should restart simulation and run setup  again")

        pastBonds = self.curBonds
        self.curBonds = self.allBonds.pop(0)  # getting current bonds
        bondsRemove = [i for i in pastBonds if i not in self.curBonds]
        bondsAdd = [i for i in self.curBonds if i not in pastBonds]
        bondsStay = [i for i in pastBonds if i in self.curBonds]
        if verbose:
            print("{0} bonds stay, {1} new bonds, {2} bonds removed".format(len(bondsStay),
                                                                            len(bondsAdd), len(bondsRemove)))
        bondsToChange = bondsAdd + bondsRemove
        bondsIsAdd = [True] * len(bondsAdd) + [False] * len(bondsRemove)
        for bond, isAdd in zip(bondsToChange, bondsIsAdd):
            ind = self.bondToInd[bond]
            paramset = self.activeParamDict if isAdd else self.inactiveParamDict
            self.bondForce.setBondParameters(ind, bond[0], bond[1], **paramset)  # actually updating bonds
        self.bondForce.updateParametersInContext(context)  # now run this to update things in the context
        return self.curBonds, pastBonds
    
## Simulation parameters and model setup

# define the region (custom code written by James)
region_size = 2000
CTCF_sites_L = np.array([574, 694, 866, 1241, 1390, 1580, 1752, 1800])
CTCF_sites_R = np.array([200, 330, 724, 1425, 1433, 1604])
sticky_elements = np.array([250, 372, 540, 745, 775, 833, 961, 1202, 1330, 1640, 1722])
num_CTCF_sites_L = len(CTCF_sites_L)
num_CTCF_sites_R = len(CTCF_sites_R)
num_sticky_elements = len(sticky_elements)
stall_probs_CTCF_L = np.array([0.6, 0.8 , 0.95, 0.1 , 0.6, 0.6, 0.8, 0.1  ])
stall_probs_CTCF_R = np.array([0.9, 0.3 , 0.95, 0.4 , 0.3, 0.4 ])

# now prepare the variables needed to actually run the simulation
chrom_size = 70000
num_regions = chrom_size // region_size

region_starts = np.arange(num_regions) * region_size
L = np.repeat(region_starts,num_CTCF_sites_L) + np.tile(CTCF_sites_L, num_regions)
R = np.repeat(region_starts,num_CTCF_sites_R) + np.tile(CTCF_sites_R, num_regions)
boundaryPauseProbL = np.tile(stall_probs_CTCF_L, num_regions)
boundaryPauseProbR = np.tile(stall_probs_CTCF_R, num_regions)

boundaryStrengthsL = np.zeros(chrom_size)
boundaryStrengthsR = np.zeros(chrom_size)
for n,s  in zip(L,boundaryPauseProbL):
    boundaryStrengthsL[n] = s 
for n,s  in zip(R,boundaryPauseProbR):
    boundaryStrengthsR[n] = s
boundaryPauseProb = 'variable'

sticky_elements_all = np.repeat(region_starts,num_sticky_elements) + np.tile(sticky_elements, num_regions)
monomer_types = np.zeros(chrom_size, dtype='int')
monomer_types[sticky_elements_all] = 1
EP_interaction_energy = 3
interaction_matrix = np.array([[0, 0],[0, EP_interaction_energy]])

# continue defining simulation parameters and model setup

initialization_steps = 10000
steps_per_sample = 100

## Dummy leftover variables 
num_cores = 1
samps_per_core = 1
####

processivities = [300]
separations = [240]

ctcf_boost_factors = [4]
longlived_fraction = [0]
longlived_boost_factor = [20]
dsb_boost_factor = [0]
targeted_loading_factor = [1]
llp = [2]
h = [0.4]

boundary_coordinates = []

# create chain with breaks
# chain segment lengths 
num_chains = 1
chain = [(x*chrom_size,(x+1)*chrom_size,0) for x in range(num_chains)]

# # get initial conformation
# conformation_folder = f'Data/3D_PolymerSimulationEquilibrationrRun/blocks_490000-499999.h5::499999'

# u = load_URI(conformation_folder)  
# data = u['pos']

positions_to_sample = np.arange(chrom_size)  # save every monomer position

# saving for polymer simulation
GPU_choice = 1

steps_per_block = 236 # number of polymer simulation steps per block (~1/40 of 9441)
total_saved_blocks =  2000000 # total number of blocks saved
save_every_x_blocks = 10 # number of blocks before a 3D configuration is saved
total_blocks_done = total_saved_blocks * save_every_x_blocks
restartBondUpdaterEveryBlocks = 100 

# block number to store all monomer positions
block_to_save_all = []

# parameters for SMC translocator
# simulation parameters for smc bonds 
smcBondWiggleDist = 0.2
smcBondDist = 0.5
smcStepsPerBlock = 1#int(1) # this scales the number of SMC substeps per simulation time step

stiff = 2
volume_density = 0.3
block = 0  # starting block 
N_monomers = chain[-1][1]

collision_rate = 1
truncated_potentials = 3

pause_prob = 1-steps_per_block/9441/4 # probability motor pauses, resulting in 125bp/s extruding speed (one motor)
# parameters for SMC translocator
translocator_initialization_steps = 10000 # for SMC translocator

date = datetime.today().strftime('%Y%m%d')
save_folder = f'/mnt/coldstorage/jjusuf/simabsloopquant_data/region_3_simulations/sim_{date}_dens{volume_density}_N{N_monomers}_C{num_chains}_EPattractiveEnergy_{EP_interaction_energy}_collisionrate_{collision_rate}_trunc{truncated_potentials}_stepsPerBlock{steps_per_block}_proc{processivities[0]}_sep{separations[0]}_boost{ctcf_boost_factors[0]}_strong_TAD_high_density'
name_format_string = save_folder+'/sim_tests_proc{}_sep{}_lp{}_H{}_boundaryPauseProb{}_ctcf{}_longlivedfraction{}_longlivedfactor{}_dsb{}_superloading{}_nreps{}'

lof_nums = list(product(llp,boundaryPauseProb,processivities,separations,[''],h,ctcf_boost_factors,longlived_fraction,longlived_boost_factor,dsb_boost_factor,targeted_loading_factor))
lof = [(x[0],x[1],name_format_string.format(x[2],x[3],x[0],x[5],x[1],x[6],x[7],x[8],x[9],x[10],int(samps_per_core*num_cores) ),
        x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10]) for x in lof_nums]

if not os.path.exists(save_folder):
    os.makedirs(save_folder)

extra_bonds = None
idx = 0 # dummy for now
run_simulation(idx, N_monomers, translocator_initialization_steps, \
                  smcStepsPerBlock, \
                  smcBondDist, smcBondWiggleDist,\
                  steps_per_block, volume_density, block, \
                  chain, extra_bonds, save_folder, save_every_x_blocks, \
                  total_saved_blocks, restartBondUpdaterEveryBlocks, GPU_choice = GPU_choice, 
                  overwrite=True,density=0.3,positions_to_sample=positions_to_sample,monomer_types=monomer_types, \
                  interaction_matrix=interaction_matrix,colrate=collision_rate,errtol=0.01, \
                  trunc=truncated_potentials,block_to_save_all=block_to_save_all,save_length=1000)