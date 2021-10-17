# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 10:04:04 2021

@author: ar-h1
"""
import numpy as np
import matplotlib.pyplot as plt 
from spectres import spectres
from astropy import units as u
import os
#from scipy.signal import find_peaks
#from astropy.modeling import models



# sampling wavelength grid uniformly in velocity space
def define_wavegrid(wavelmin,wavelmax,res):
   #c=299792458.
   dx=np.log(1.+1./res)
   x=np.arange(np.log(wavelmin),np.log(wavelmax),dx)
   wavelength=np.exp(x)
   #waveno=1e4/wavelength
   return wavelength#,waveno


def bin_spectrum_fast(wl_native, spectrum_native, wl_start, wl_end, R_bin):
    # Create binned wavelength grid at resolution R_bin
    delta_log_wl_bins = 1.0/R_bin
    N_wl_bins = (np.log(wl_end) - np.log(wl_start)) / delta_log_wl_bins
    N_wl_bins = np.around(N_wl_bins).astype(np.int64)
    log_wl_binned = np.linspace(np.log(wl_start), np.log(wl_end), N_wl_bins)
    wl_binned = np.exp(log_wl_binned)
    spectrum_binned = spectres(wl_binned, wl_native, spectrum_native)
    return wl_binned, spectrum_binned


#ModelName = 'Fe_0.50_+0.0_0.55'
#ModelName = 'Fe_0.50_-1.0_0.55'
#ModelName = 'Al_0.50_+0.0_0.55'
#ModelName = 'Al_0.50_+2.3_0.55'
#ModelName = 'CaII_0.50_+0.0_0.55'
#ModelName = 'CO_all_iso_0.50_+0.0_0.55'
#ModelName = 'Cr_0.50_+0.0_0.55'
#ModelName = 'FeH_main_iso_0.50_+0.0_0.55'
#ModelName = 'K_0.50_+0.0_0.55'
#ModelName = 'Mg_0.50_+0.0_0.55'
#ModelName = 'Na_0.50_+0.0_0.55'
#ModelName = 'Ti_0.50_+0.0_0.55'
#ModelName = 'TiO_48_Exomol_McKemmish_0.50_+0.0_0.55'
#ModelName = 'TiO_48_Plez_0.50_+0.0_0.55'
#ModelName = 'TiO_all_iso_Plez_0.50_+0.0_0.55'
#ModelName = 'VO_0.50_+0.0_0.55'
#ModelName = 'FeII_0.50_+0.0_0.55'
#ModelName = 'Ca_0.50_+0.0_0.55'

#ModelDescription = '0.50_-1.0_0.55'
#ModelDescription = '0.50_+1.0_0.55'
#ModelDescription = '0.50_+1.7_0.55'
#ModelDescription = '0.50_+2.0_0.55'
#ModelDescription = '0.50_+2.3_0.55'


Species = 'Fe'

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


            #HighResModelPath = 'F:/KELT-9b_CARMENES_emission_data/ModelSpectra/ModelDescription_0.50_+0.0_0.55/RotationalBroadening'
            #HighResModelPath = '/media/andrew/easystore/KELT-9b_CARMENES_emission_data/ModelSpectra/ModelDescription_0.50_+0.0_0.55/RotationalBroadening'
            #HighResModelPath = '/media/andrew/easystore/KELT-9b_CARMENES_emission_data/ModelSpectra/ModelDescription_%s/RotationalBroadening'%(ModelDescription)
            HighResModelPath = 'F:/KELT-9b_CARMENES_emission_data/ScriptsForPetitRADTRANS/ModelsFromLoop/Fe/R2e5_RotBroad6.63'

            
            
            
            
            #HighResModelPath = 'F:\KELT-9b_CARMENES_emission_data\ModelSpectra\ModelDescription_0.50_+2.3_0.55\RotationalBroadening'
            #HighResModelPath = 'F:/KELT-9b_CARMENES_emission_data/ModelSpectra/ModelDescription_%s/RotationalBroadening'%(ModelDescription)
            
            
            
            LowResSavePath = 'F:/KELT-9b_CARMENES_emission_data/ScriptsForPetitRADTRANS/ModelsFromLoop/Fe/CarmenesRes_RotBroad6.63'
            # if not os.path.exists(LowResSavePath):
            #     os.makedirs(LowResSavePath)
            
            
            
            
            
            
            
            
            
            HighResModel = np.load('%s/K9b_%s_%s_Vrot6.63_pRT_flux_per_Hz.npy'%(HighResModelPath,Species,ModelDescription))
            SaveName = 'K9b_%s_%s_Vrot6.63_CarmenesRes_pRT_flux_per_Hz.npy'%(Species,ModelDescription)
            
            #HighResModel = np.load('/media/andrew/easystore/KELT-9b_CARMENES_emission_data/ModelSpectra/BB_10170K_R1e6_erg_Per_s_Per_cm2_Per_Hz.npy')
            
            
            InitialR = np.mean(HighResModel[1:,0]/np.diff(HighResModel[:,0]))
            
            CARMENES_VIS_Res = 94600
            CARMENES_NIR_Res = 80400
            
            SeperatingWavelength = 0.96
            
            VIS_wl_binned, VIS_spectrum_binned = bin_spectrum_fast(HighResModel[:,0],HighResModel[:,1],HighResModel[0,0],HighResModel[-1,0],CARMENES_VIS_Res)
            
            VIS_wl_binned = VIS_wl_binned[1:-1]
            VIS_spectrum_binned = VIS_spectrum_binned[1:-1]
            
            RbinnedVis = VIS_wl_binned[1:]/np.diff(VIS_wl_binned)
            
            NIR_wl_binned, NIR_spectrum_binned = bin_spectrum_fast(HighResModel[:,0],HighResModel[:,1],HighResModel[0,0],HighResModel[-1,0],CARMENES_NIR_Res)
            
            NIR_wl_binned = NIR_wl_binned[1:-1]
            NIR_spectrum_binned = NIR_spectrum_binned[1:-1]
            
            NeededVisIndices = np.where(VIS_wl_binned<=SeperatingWavelength)
            NeededNirIndices = np.where(NIR_wl_binned>SeperatingWavelength)
            
            NumVisPoints = len(NeededVisIndices[0])
            NumNirPoints = len(NeededNirIndices[0])
            
            LowResSpec = np.zeros((NumVisPoints+NumNirPoints,2))
            
            LowResSpec[0:NumVisPoints,0] = VIS_wl_binned[NeededVisIndices]
            LowResSpec[0:NumVisPoints,1] = VIS_spectrum_binned[NeededVisIndices]
            
            LowResSpec[NumVisPoints:,0] = NIR_wl_binned[NeededNirIndices]
            LowResSpec[NumVisPoints:,1] = NIR_spectrum_binned[NeededNirIndices]
            
            
            
            
            RLowResCombined = LowResSpec[1:,0]/np.diff(LowResSpec[:,0])
            
            
            np.save('%s/%s'%(LowResSavePath,SaveName),LowResSpec)
            #np.save('/media/andrew/easystore/KELT-9b_CARMENES_emission_data/ModelSpectra/BB_10170K_CarmenesRes_erg_Per_s_Per_cm2_Per_Hz.npy',LowResSpec)
            
            RLowResCombinedVisPart = LowResSpec[1:NumVisPoints,0]/np.diff(LowResSpec[0:NumVisPoints,0])
            
            RLowResCombinedNirPart = LowResSpec[NumVisPoints+1:,0]/np.diff(LowResSpec[NumVisPoints:,0])
            
            
            plt.plot(HighResModel[:,0],HighResModel[:,1],label='high res')
            plt.plot(LowResSpec[:,0],LowResSpec[:,1],label='low res')
            plt.legend()
            
            
            plt.figure()
            plt.title('Combined effective resolution')
            plt.plot(RLowResCombined)
            
            assert len(np.unique(LowResSpec[:,0])) == len(LowResSpec[:,0]) ## check no repeated wavelengths 
            
            assert False not in (LowResSpec[:,0] == np.sort(LowResSpec[:,0]))