#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 18:14:37 2021

@author: andrew
"""

import numpy as np
import matplotlib.pyplot as plt 
from astropy import units as u
from scipy import interpolate
from scipy.signal import butter, sosfiltfilt
from scipy.signal import savgol_filter
import os
from scipy.signal import find_peaks

def FindEnvelopeWithBins(wave,depth,PointsPerBin=800,numbins=None,PolyOrder=50):
    
    if numbins == None:
        numbins = int(len(wave)/PointsPerBin)
    
    

    bins = np.linspace(min(wave), max(wave), numbins+1)
    digitized = np.digitize(wave, bins)
    binned_mins = [np.nanmin(depth[digitized==i]) for i in range(1, len(bins))]  
    
    BinMiddles = (bins[0:-1] + bins[1:])/2
    
    coeff = np.polyfit(BinMiddles,binned_mins,PolyOrder)
    LowOrderPolyFit = np.polyval(coeff, wave)
    
    return LowOrderPolyFit, BinMiddles, binned_mins


a = np.load('/media/andrew/easystore2/KELT-9b_CARMENES_emission_data/ScriptsForPetitRADTRANS/KELT-9b_emission_for_NEID_prop/Fe/NEID_res/K9b_Fe_0.50_+0.0_0.55_NEID_res.npy')



env,binmiddles,binnedmins = FindEnvelopeWithBins(wave = np.arange(len(a[:,0])),depth=a[:,1],PointsPerBin=100, numbins=None, PolyOrder=3)

flat = a[:,1] - env
flat = flat[1:-1]

x = flat

peaks, _ = find_peaks(x, prominence=0.5e-5)

plt.plot(x)

plt.plot(peaks, x[peaks], "x")

plt.plot(np.zeros_like(x), "--", color="gray")

plt.show()

# 

# 

# PointsPerBin=100
# numbins=None
# PolyOrder=1

# if numbins == None:
#     numbins = int(len(wave)/PointsPerBin)



# bins = np.linspace(min(wave), max(wave), numbins+1)
# digitized = np.digitize(wave, bins)
# binned_mins = [np.nanmin(depth[digitized==i]) for i in range(1, len(bins))]  

# BinMiddles = (bins[0:-1] + bins[1:])/2

# coeff = np.polyfit(BinMiddles,binned_mins,PolyOrder)
# LowOrderPolyFit = np.polyval(coeff, wave)



