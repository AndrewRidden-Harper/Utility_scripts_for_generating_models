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

def butterworth(x, order, freq, filt_type='highpass'):
	"""
	Applies a high-pass Butterworth filter, with a given order and 
	cut-off frequency, to the given model.
	"""
	butterfilt = butter(order, freq, btype=filt_type, output='sos')
	x_filtered = sosfiltfilt(butterfilt, x)
	return x_filtered

def remove_env(wave, spec, px):
	"""
	Subtracts the lower envelope from a model spectrum by finding 
	the minimum value in the given stepsize, then interpolating.
	"""
	low_wave, low_spec = [], []
	for i in range(len(spec)/px - 1):
		idx = np.nanargmin(spec[i*px:(i+1)*px])
		low_spec.append(spec[idx+i*px])
		low_wave.append(wave[idx+i*px])
	interp = interp1d(low_wave, low_spec, fill_value='extrapolate')
	envelope = interp(wave)
	corrected = spec - envelope
	return corrected

def FindEnvelopeWithBins(wave,depth,PointsPerBin=800,numbins=None,PolyOrder=50):
    
    if numbins == None:
        numbins = int(len(wave)/PointsPerBin)
    
    

    bins = np.linspace(min(wave), max(wave), numbins+1)
    digitized = np.digitize(wave, bins)
    binned_mins = [np.min(depth[digitized==i]) for i in range(1, len(bins))]  
    
    BinMiddles = (bins[0:-1] + bins[1:])/2
    
    coeff = np.polyfit(BinMiddles,binned_mins,PolyOrder)
    LowOrderPolyFit = np.polyval(coeff, wave)
    
    return LowOrderPolyFit, BinMiddles, binned_mins

c=299792.458 #speed of light (km/s)

SystemicVelocity = -20.567  ## kms for KELT-9b from NASA Exoplanet Archive. Gaudi et al. 2017 
PlanetInclination = 87.2*np.pi/180   ## for KELT-9b, Ahlers et al. 2020 
Rplanet = 1.891*u.Rjup 
Rstar = 2.36*u.Rsun

RplanetOverRstarSquared = ((Rplanet/Rstar)**2).decompose().value


#FirstPartOfLoadPath = '../CrossCorrelationDataAndProcessing'
#FirstPartOfLoadPath = 'F:'

#ModelForCrossCor = 'KELT9b_Fe_0.50_+0.0_0.55_Vrot6.63_CarmenesRes'



ModelSpecies = 'Fe'


#ModelDescription = '0.50_-1.0_0.55'
#ModelDescription = '0.50_+1.7_0.55'
#ModelDescription = '0.50_+2.3_0.55'

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

            #ModelForCrossCor = 'KELT9b_%s_%s_Vrot6.63_CarmenesRes'%(ModelSpecies,ModelDescription)
            ModelForCrossCor =  'K9b_%s_%s_Vrot6.63_CarmenesRes_pRT_flux_per_Hz.npy'%(ModelSpecies,ModelDescription)
            
            
            SavePath = '../ModelSpectraLLFromLoop'
            
            if not os.path.exists(SavePath):
                os.makedirs(SavePath)
            
            
            ##StellarWaveFlux = np.load('%s/KELT-9b_CARMENES_emission_data/ModelSpectra/BB_10170K_CarmenesRes_erg_Per_s_Per_cm2_Per_Hz.npy'%(FirstPartOfLoadPath))         
            StellarWaveFlux = np.load('../ModelSpectra/BB_10170K_CarmenesRes_erg_Per_s_Per_cm2_Per_Hz.npy')         
            
            
            
            
            
            #PlanetWaveFlux = np.load('%s/KELT-9b_CARMENES_emission_data/ModelSpectra/%s_pRT_flux_per_Hz.npy'%(FirstPartOfLoadPath,ModelForCrossCor))
            
            #PlanetWaveFlux = np.load('%s/KELT-9b_CARMENES_emission_data/ModelSpectra/ModelDescription_%s/RotationalBroadening/CARMENES_resolution/%s_pRT_flux_per_Hz.npy'%(FirstPartOfLoadPath,ModelDescription,ModelForCrossCor))
            #PlanetWaveFlux = np.load('../ModelSpectra/ModelDescription_%s/RotationalBroadening/CARMENES_resolution/%s_pRT_flux_per_Hz.npy'%(ModelDescription,ModelForCrossCor))
            
            
            
            PlanetWaveFlux = np.load('ModelsFromLoop/%s/CarmenesRes_RotBroad6.63/%s'%(ModelSpecies, ModelForCrossCor))

            
            
            
            #F:\KELT-9b_CARMENES_emission_data\ModelSpectra\ModelDescription_0.50_+0.0_0.55\RotationalBroadening\CARMENES_resolution
            
            PlanetWaveFlux = PlanetWaveFlux[np.argsort(PlanetWaveFlux[:,0])]
            StellarWaveFlux = StellarWaveFlux[np.argsort(StellarWaveFlux[:,0])]
            
            StellarInterpObject = interpolate.interp1d(StellarWaveFlux[:,0], StellarWaveFlux[:,1], bounds_error=False, fill_value='extrapolate')
            InterpedStellarFlux = StellarInterpObject(PlanetWaveFlux[:,0])
            
            
            
            
            ratio = (RplanetOverRstarSquared)*(PlanetWaveFlux[:,1]/InterpedStellarFlux)
            
            
            #mspec_bf = butterworth(ratio, 1, 0.01)
            #mspec_bf = butterworth(ratio, 1, 0.03)
            #mspec_bf = butterworth(ratio, 1, 0.05)
            
            
            
            # plt.figure()
            # plt.plot(mspec_bf*1e6)
            
            ###############################################
            
            ModelSpectrumFlattening_PointsPerBin = 50 
            ModelSpectrumFlattening_PolyOrder = 2
            NumParts = 300
            
            OffsetTranDepNoNan = ratio
            
            
            xIndexArray = np.arange(len(OffsetTranDepNoNan))
            
            #NumParts = 300 ## Original 
            #NumParts = 600
            
            PointsPerPart = int(len(OffsetTranDepNoNan)/NumParts)
            
            fl = np.copy(OffsetTranDepNoNan)
            totalenv = np.copy(OffsetTranDepNoNan)
            
            
            for PartIndex in range(NumParts):
                
                IndexLims = np.array([PartIndex,PartIndex+1])*PointsPerPart
                
                PointsInPart = IndexLims[1] - IndexLims[0]
                
                print('%d points in part %d'%(PointsInPart,PartIndex))
                
                if PartIndex == (NumParts - 1):
                    
                    if IndexLims[1] > len(OffsetTranDepNoNan):
                        print('Last part needs fewer points to cover spectrum')
                        IndexLims[1] = len(OffsetTranDepNoNan)
                        
                    if IndexLims[1] < len(OffsetTranDepNoNan) - 1:
                        print('Points are missing at the end of the last part')
                        IndexLims[1] = len(OffsetTranDepNoNan)
                        PointsInPart = IndexLims[1] - IndexLims[0]
                        print('%d points in part %d'%(PointsInPart,PartIndex))
                
                xpart = xIndexArray[IndexLims[0]:IndexLims[1]]
                rpart = OffsetTranDepNoNan[IndexLims[0]:IndexLims[1]]
                
            
            
                env,binmiddles,binnedmins = FindEnvelopeWithBins(xpart,rpart,PointsPerBin=ModelSpectrumFlattening_PointsPerBin, numbins=None, PolyOrder=ModelSpectrumFlattening_PolyOrder)
                
                # if (LoadFileName == 'GJ486b_O2_700K') & (ResStr == 'R70k') & (PartIndex == 331):
                #     env,binmiddles,binnedmins = FindEnvelopeWithBins(xpart,rpart,PointsPerBin=ModelSpectrumFlattening_PointsPerBin, numbins=None, PolyOrder=1)
                #     env = np.ones_like(rpart)*np.min(rpart)
            
                
                
                # plt.figure()
                # plt.plot(xpart,rpart)    
                # plt.plot(xpart,env)
                # plt.title('%s/Part %d of %d'%(LoadFileName,PartIndex,NumParts))
            
            
               
                # fl[IndexLims[0]:IndexLims[1]] = rpart - env #+ 1
                totalenv[IndexLims[0]:IndexLims[1]] = env
                # plt.plot(binmiddles,binnedmins,'o')
                
                
                #fl[IndexLims[0]:IndexLims[1]] -= np.sort(fl[IndexLims[0]:IndexLims[1]])[100]
                #fl[IndexLims[0]:IndexLims[1]] -= np.min(fl[IndexLims[0]:IndexLims[1]])
                
            SmoothTotalEnv1 = savgol_filter(totalenv, 101, 3)
            
            # SmoothTotalEnv2 = savgol_filter(totalenv, 1001, 3)
            # SmoothTotalEnv3 = savgol_filter(totalenv, 5001, 3)
            # SmoothTotalEnv4 = savgol_filter(totalenv, 10001, 3)
            
            
            
            # plt.figure()
            # plt.plot(totalenv,label='unsmoothed')
            #plt.plot(SmoothTotalEnv1,label='1')
            # plt.plot(SmoothTotalEnv2,label='2')
            # plt.plot(SmoothTotalEnv3,label='3')
            # plt.plot(SmoothTotalEnv4,label='4')
            # plt.legend()
            
            
            SmoothTotalEnv = SmoothTotalEnv1
            #SmoothTotalEnv = totalenv
            
            
            fl -= SmoothTotalEnv
            
            plt.figure()
            plt.plot(ratio)
            plt.plot(SmoothTotalEnv1,label='1')
            
            plt.figure()
            plt.plot(fl*1e6)
            
            
            
            FlatTranDepArray = np.zeros((len(fl),2))
            FlatTranDepArray[:,0] = PlanetWaveFlux[:,0]
            FlatTranDepArray[:,1] = fl
            
            
            #np.save('%s/Flat_%s_%s'%(SavePath,LoadFileName,ResStr),FlatTranDepArray)
            np.save('%s/ContSub_%s.npy'%(SavePath,ModelForCrossCor),FlatTranDepArray)



