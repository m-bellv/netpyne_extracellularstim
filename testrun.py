#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 12:23:47 2022

@author: mbellvila
"""

#Testing extracellular stim with morphology files - only using 2 neuron types

from netpyne import sim
from netpyne import specs
import math
import numpy as np
import neuron
import os
from re import search
import sys
import tkinter as tk
from tkinter import filedialog

# object of class NetParams to store the network parameters
netParams = specs.NetParams()
# object of class SimConfig to store the simulation configuration
simConfig = specs.SimConfig()


# Simulation parameters
simConfig.duration = 3000  # Duration of the simulation, in ms
simConfig.dt = 0.05  # Internal integration timestep to use
# Seeds for randomizers (connectivity, input stimulation and cell locations)
simConfig.seeds = {'conn': 1, 'stim': 1, 'loc': 1}
simConfig.createNEURONObj = 1  # create HOC objects when instantiating network
# create Python structure (simulator-independent) when instantiating network
simConfig.createPyStruct = 1
simConfig.verbose = 0  # Whether to write diagnostic information on events


# Recording
simConfig.recordCells = []  # list of cells to record from
# 'V':{'sec':'soma','loc':0.5,'var':'v'}}
simConfig.recordTraces = {'V': {'sec': 'soma', 'loc': 0.5, 'var': 'v'}}
#simConfig.recordLFP = [[0.5*netParams.sizeX, 0.215*netParams.sizeY, 0.5*netParams.sizeZ]]

simConfig.recordStim = False  # record spikes of cell stims
# Step size in ms to save data (eg. V traces, LFP, etc)
simConfig.recordStep = 0.1
simConfig.hParams = {'celsius': 34, 'v_init': -80}


# Saving
simConfig.filename = 'M1_ynorm_izhi'  # Set file output name
#simConfig.saveDataInclude = (['simData'])
simConfig.saveFileStep = 3000  # step size in ms to save data to disk
simConfig.savePickle = False  # save to pickle file
simConfig.saveJson = False  # save to json file
simConfig.saveMat = False  # save to mat file
simConfig.saveTxt = False  # save to txt file
simConfig.saveDpk = False  # save to .dpk pickled file
simConfig.saveHDF5 = False  # save to HDF5 file

# Analysis and plotting
simConfig.analysis['plotRaster'] = {'orderInverse': True, 'figSize': (
    12, 10), 'lw': 0.3, 'markerSize': 6, 'marker': '.', 'dpi': 300, 'saveFig': 'raster.png'}  # Whether or not to plot a raster
simConfig.analysis['plotTraces'] = {'include': [('IT_L23', 10), ('IT_L23', 1), (
    'IT_L23', 18), ('IT_L23', 11), ]}  # plot recorded traces for this list of cells
simConfig.analysis['plot2Dnet'] = {'showConns': False}
simConfig.analysis['plotLFP'] = False


###############################################################################
# NETWORK PARAMETERS
###############################################################################


# General network parameters
netParams.scale = 1  # Scale factor for number of cells
# x-dimension (horizontal length) size in um
netParams.sizeX = 10*math.sqrt(netParams.scale)
# y-dimension (vertical height or cortical depth) size in um
netParams.sizeY = 10*math.sqrt(netParams.scale)
# z-dimension (horizontal depth) size in um
netParams.sizeZ = 10*math.sqrt(netParams.scale)


# General connectivity parameters
netParams.scaleConnWeight = 0.0003  # Connection weight scale factor
# Connection weight scale factor for NetStims
netParams.scaleConnWeightNetStims = 0.25
netParams.defaultDelay = 2.0  # default conn delay (ms)
netParams.propVelocity = 100.0  # propagation velocity (um/ms)
# length constant (lambda) for connection probability decay (um)
netParams.probLambda = 100.0
netParams.defaultThreshold = 0.0


print('Select folder containing morphology files')
#Morphology files obtained from: https://github.com/suny-downstate-medical-center/netpyne/tree/development/examples/M1detailed

root = tk.Tk()
root.withdraw()

morph_folder_path = filedialog.askdirectory()
root.update() 
print(morph_folder_path)

sys.path.append(morph_folder_path)
netParams.loadCellParamsRule(
    label='IT_L23', fileName=morph_folder_path +'/IT2_reduced_cellParams.pkl')
netParams.loadCellParamsRule(
    label='PV_L23', fileName= morph_folder_path + '/PV_simple_cellParams.pkl')



#Importing cell params and adding extracellular mechanism
netParams.cellParams['IT_L23']['conds'] = {'cellType': 'IT_L23'}
ITL23seg = []
ITL23dendlist = []
for secName in netParams.cellParams['IT_L23']['secs']:
    if search('dend', secName):
        ITL23seg.append(secName)
        ITL23dendlist.append(secName)
    if search('soma', secName):
        ITL23seg.append(secName)
    if search('axon', secName):
        ITL23seg.append(secName)
for segname in ITL23seg:
    netParams.cellParams['IT_L23']['secs'][segname]['mechs']['extracellular'] = {
    }
netParams.cellParams['IT_L23']['conds']['ynorm'] = [0, 1]


netParams.cellParams['PV_L23']['conds'] = {'cellType': 'PV_L23'}
PVL23seg = []
for secName in netParams.cellParams['PV_L23']['secs']:
    if search('dend', secName):
        PVL23seg.append(secName)
    if search('soma', secName):
        PVL23seg.append(secName)
    if search('axon', secName):
        ITL23seg.append(secName)
for segname in PVL23seg:
    netParams.cellParams['PV_L23']['secs'][segname]['mechs']['extracellular'] = {
    }
netParams.cellParams['PV_L23']['conds']['ynorm'] = [0, 1]


netParams.popParams['ITL23'] = {'cellModel': 'HH_reduced', 'cellType': 'IT_L23', 'numCells': 20}
netParams.popParams['PVL23'] = {'cellModel': 'HH_simple', 'cellType': 'PV_L23', 'numCells': 20}


# Synaptic mechanism parameters
netParams.synMechParams['AMPA'] = {
    'mod': 'MyExp2SynBB', 'tau1': 0.05, 'tau2': 5.3, 'e': 0}  # AMPA
netParams.synMechParams['NMDA'] = {
    'mod': 'MyExp2SynNMDABB', 'tau1NMDA': 15, 'tau2NMDA': 150, 'e': 0}  # NMDA
netParams.synMechParams['exc'] = {
    'mod': 'Exp2Syn', 'tau1': 0.05, 'tau2': 5.3, 'e': 0}
netParams.synMechParams['GABAB'] = {
    'mod': 'MyExp2SynBB', 'tau1': 3.5,    'tau2': 260.9,  'e': -93}
netParams.synMechParams['GABAA'] = {
    'mod': 'MyExp2SynBB', 'tau1': 0.07,   'tau2': 18.2,   'e': -80}
netParams.synMechParams['GABAASlow'] = {
    'mod': 'MyExp2SynBB', 'tau1': 2,      'tau2': 100,    'e': -80}
netParams.synMechParams['GABAASlowSlow'] = {
    'mod': 'MyExp2SynBB', 'tau1': 200,    'tau2': 400,    'e': -80}

ESynMech = ['AMPA', 'NMDA']
SOMESynMech = ['GABAASlow', 'GABAB']
SOMISynMech = ['GABAASlow']
PVSynMech = ['GABAA']


# Stimulation parameters
netParams.stimSourceParams['background_E'] = {
    'type': 'NetStim',  'rate': 20, 'noise': 1.0}  # background inputs to Exc
netParams.stimSourceParams['background_I'] = {
    'type': 'NetStim',  'rate': 10, 'noise': 1.0}  # background inputs to Inh

netParams.stimTargetParams['bgE->IT'] = {'source': 'background_E', 'conds': {'pop': ['ITL23']},
                                         'secList': ITL23seg, 'allSegs': True, 'synMech': ESynMech, 'weight': 0.1, 'delay': '2+normal(5,3)'}

netParams.stimTargetParams['bgI->PV'] = {'source': 'background_E', 'conds': {'pop': ['PVL23']},
                                         'secList': PVL23seg, 'allSegs': True, 'synMech': ESynMech, 'weight': 0.05, 'delay': '2+normal(5,3)'}


# List of connectivity rules/params
netParams.ItoIweight = 0.5

netParams.addConnParams(None, {'preConds': {'cellType': ['IT_L23']},
                               'postConds': {'cellType': 'IT_L23'},
                               'synMech': ESynMech,
                               'probability': 0.1,
                               'weight': 0.2,
                               'synMechWeightFactor': [0.5, 0.5],
                               'delay': 0})


netParams.addConnParams(None, {'preConds': {'cellType': 'IT_L23'},
                               'postConds': {'cellType': 'PV_L23'},
                               'synMech': ESynMech,
                               'probability': 0.1,
                               'synMechWeightFactor': [0.5, 0.5],
                               'weight': 0.2,
                               'delay': 0})


netParams.addConnParams(None, {'preConds': {'cellType': 'PV_L23'},
                               'postConds': {'ynorm': [0, 1]},
                               'synMech': PVSynMech,
                               'probability': '1.0 * exp(-dist_3D/probLambda)',
                               'weight': 'ItoIweight',
                               'delay': 0})


# Dictionary of annotations
netParams.annots = {}


sim.create(netParams=netParams)
sim.net.defineCellShapes()

# The parameters of the extracellular point current source
acs_params = {'position': [0, 0, 0.],  # um
              'amp': 10000000000000000.,  # uA,
              'stimstart': 1000,  # ms
              'stimend': 2000,  # ms
              'frequency': 10,  # Hz
              'sigma': 0.3  # decay constant
              }


def insert_v_ext(cell, v_ext, t_ext):

    cell.t_ext = neuron.h.Vector(t_ext)
    cell.v_ext = []
    for v in v_ext:
        cell.v_ext.append(neuron.h.Vector(v))

    # play v_ext into e_extracellular reference
    i = 0
    cell.v_ext[i].play(cell.secs['soma']['hObj'](
        0.5)._ref_e_extracellular, cell.t_ext)


def make_extracellular_stimuli(acs_params, cell):
    """ Function to calculate and apply external potential """
    x0, y0, z0 = acs_params['position']
    ext_field = np.vectorize(lambda x, y, z: 1 / (4 * np.pi *
                                                  acs_params['sigma'] * np.sqrt((x0 - x)**2 + (y0 - y)**2 + (z0 - z)**2)))

    stimstart = acs_params['stimstart']
    stimend = acs_params['stimend']
    stimdif = stimend-stimstart

    # MAKING THE EXTERNAL FIELD
    n_tsteps = int(stimdif / simConfig.dt + 1)
    n_start = int(stimstart/simConfig.dt)
    n_end = int(stimend/simConfig.dt + 1)
    t = np.arange(start=n_start, stop=n_end) * simConfig.dt
    pulse = acs_params['amp'] * 1000. * \
        np.sin(2 * np.pi * acs_params['frequency'] * t / 1000)

    #v_cell_ext = np.zeros((cell.secs['soma']['hObj'].nseg, n_tsteps))
    v_cell_ext = np.zeros((1, n_tsteps))
    #v_cell_ext[:, :] = ext_field(cell.getSomaPos()[0], cell.getSomaPos()[1], cell.getSomaPos()[2]).reshape(cell.secs['soma']['hObj'].nseg, 1) * pulse.reshape(1, n_tsteps)
    v_cell_ext[:, :] = ext_field(cell.getSomaPos()[0], cell.getSomaPos(
    )[1], cell.getSomaPos()[2]).reshape(1, 1) * pulse.reshape(1, n_tsteps)
    insert_v_ext(cell, v_cell_ext, t)

    return ext_field, pulse

#Add extracellular stim
for c in range(len(sim.net.cells)):
    ext_field, pulse = make_extracellular_stimuli(acs_params, sim.net.cells[c])


sim.simulate()
sim.analyze()