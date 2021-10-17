# -*- coding: utf-8 -*-
"""
Created on Sun Jul 11 19:22:10 2021

@author: ar-h1
"""

import numpy as np
import matplotlib.pyplot as plt 
from astropy import units as u
from scipy import interpolate
from scipy.signal import butter, sosfiltfilt
from scipy.signal import savgol_filter
import os
from spectres import spectres


def bin_spectrum_fast(wl_native, spectrum_native, wl_start, wl_end, R_bin):
    # Create binned wavelength grid at resolution R_bin
    delta_log_wl_bins = 1.0/R_bin
    N_wl_bins = (np.log(wl_end) - np.log(wl_start)) / delta_log_wl_bins
    N_wl_bins = np.around(N_wl_bins).astype(np.int64)
    log_wl_binned = np.linspace(np.log(wl_start), np.log(wl_end), N_wl_bins)
    wl_binned = np.exp(log_wl_binned)
    spectrum_binned = spectres(wl_binned, wl_native, spectrum_native)
    return wl_binned, spectrum_binned


SpeciesName = 'Fe'
NewRes = 2e5

SavePath = 'F:/KELT-9b_CARMENES_emission_data/ScriptsForPetitRADTRANS/ModelsFromLoop/Fe/R2e5'

# if not os.path.exists(SavePath):
#     os.makedirs(SavePath)

Fc_list = [0.25, 0.5, 0.75, 1.0]
FeH_list = [-1.0, 0.0, 1.0, 1.7, 2.0, 2.3]  ### metallicities (0.1, 1, 10, 50, 100, 200; × solar) 
CToO_list = [0.35, 0.55, 0.7, 0.75, 1.0, 1.5] ### Note that 0.55 is solar 

# Fc_list = [0.25]
# FeH_list = [-1.0]  ### metallicities (0.1, 1, 10, 50, 100, 200; × solar) 
# CToO_list = [0.35] ### Note that 0.55 is solar 

for FcIndex in range(len(Fc_list)):
    for FeHIndex in range(len(FeH_list)):
        for CToOIndex in range(len(CToO_list)):
            
            Fc = Fc_list[FcIndex]
            FeH = FeH_list[FeHIndex]
            CToO = CToO_list[CToOIndex]
            
            if FeH < 0:
                ModelDescription = '%.2f_%.1f_%.2f'%(Fc, FeH, CToO)
            
            if FeH >= 0:
                ModelDescription = '%.2f_+%.1f_%.2f'%(Fc, FeH, CToO)
                
                
            pRT_spec = np.load('F:/KELT-9b_CARMENES_emission_data/ScriptsForPetitRADTRANS/ModelsFromLoop/Fe/R1e6/K9b_%s_%s.npy'%(SpeciesName,ModelDescription))
            
            
            wl_binned, spectrum_binned = bin_spectrum_fast(pRT_spec[:,0],pRT_spec[:,1],np.min(pRT_spec[:,0]),np.max(pRT_spec[:,0]),NewRes)
            
            # plt.plot(pRT_spec[:,0],pRT_spec[:,1])
            # plt.plot(wl_binned, spectrum_binned)
            
            wl_binned = wl_binned[1:-1]
            spectrum_binned = spectrum_binned[1:-1]
            
            SaveArray = np.hstack((wl_binned[:,np.newaxis],spectrum_binned[:,np.newaxis]))

            np.save('%s/K9b_%s_%s.npy'%(SavePath, SpeciesName, ModelDescription),SaveArray)
            
            ##raise Exception
            
            