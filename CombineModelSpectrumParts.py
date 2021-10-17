#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 18:22:48 2020

@author: andrew
"""

import numpy as np
import matplotlib.pyplot as plt 
import astropy.io.fits as pyfits
from astropy import units as u
from scipy import interpolate
import spectres
from astropy.modeling.models import BlackBody
from astropy import constants as const
import time 
import os 
from scipy import interpolate 

def Ryan_bin_spectrum_fast(wl_native, spectrum_native, wl_start, wl_end, R_bin):
    # Create binned wavelength grid at resolution R_bin
    delta_log_wl_bins = 1.0/R_bin
    N_wl_bins = (np.log(wl_end) - np.log(wl_start)) / delta_log_wl_bins
    N_wl_bins = np.around(N_wl_bins).astype(np.int64)
    log_wl_binned = np.linspace(np.log(wl_start), np.log(wl_end), N_wl_bins)
    wl_binned = np.exp(log_wl_binned)
    spectrum_binned = spectres.spectres(wl_binned, wl_native, spectrum_native)
    return wl_binned, spectrum_binned

def IntegrateF_nu_Or_F_lambda(wave_um, planet_flux, UseSpectres=True):

    WaveBinLimits = np.zeros((len(wave_um)+1))
    
    WaveBinLimits[0] = wave_um[0]
    WaveBinLimits[-1] = wave_um[-1]
    
    MiddlePoints = (wave_um[0:-1] + wave_um[1:])/2    
    
    WaveBinLimits[1:-1] = MiddlePoints    
    
    if not UseSpectres:
        f = interpolate.interp1d(wave_um, planet_flux)   
        FluxAtMiddlePoints = f(WaveBinLimits)
    
    if UseSpectres:
        FluxAtMiddlePoints = spectres.spectres(WaveBinLimits,wave_um,planet_flux)
    
    Integration = np.empty_like(planet_flux)
    
    for i in range(len(wave_um)):
        
        Integration[i] = np.trapz(FluxAtMiddlePoints[i:i+2], WaveBinLimits[i:i+2])   
        
    return Integration


# sampling wavelength grid uniformly in velocity space
def define_wavegrid(wavelmin,wavelmax,res):
   c=299792458.
   dx=np.log(1.+1./res)
   x=np.arange(np.log(wavelmin),np.log(wavelmax),dx)
   wavelength=np.exp(x)
   #waveno=1e4/wavelength
   return wavelength#,waveno

#resolving_power = 94600

starttime = time.time()

c_kms = const.c.to('km/s')

# SystemicVelocity = -20.567*u.km/u.s  ## kms for KELT-9b from NASA Exoplanet Archive
# PlanetRV = -200*u.km/u.s 
# BarycentricRVcorrection = 10*u.km/u.s 

# NettStellarRV = SystemicVelocity - BarycentricRVcorrection
# NettPlanetRV = PlanetRV + NettStellarRV

NettStellarRV = 0*u.km/u.s
NettPlanetRV = 0*u.km/u.s

NumParts = 20

PartsList = []

VStackOfSpec = np.empty((3))

#name = 'Fe'
#name = 'PseudoFeII'
#name = 'Al'
#name = 'Ca'
#name = 'CaII'
#name = 'CO_all_iso'
#name = 'H2O_main_iso'
#name = 'K'
#name = 'Mg'
#name = 'Na'
#name = 'SiO_main_iso'
#name = 'Ti'
#name = 'pseudoAlII'
#name = 'Cr'
#name = 'TiO_all_iso_Plez' 

#name = 'FeH_main_iso'
#name = 'VO'
#name = 'TiO_48_Exomol_McKemmish'
#name = 'TiO_48_Plez'
#name = 'FeII'
name = 'FeII_UsingFeI'



#ModelDescription = '0.25_+0.0_0.55'
#ModelDescription = '0.25_+1.0_0.55'
#ModelDescription = '0.25_+2.0_0.55'
#ModelDescription = '0.25_+2.3_0.55'

#ModelDescription = '0.50_-1.0_0.55'
ModelDescription = '0.50_+0.0_0.55'
#ModelDescription = '0.50_+1.0_0.55'
#ModelDescription = '0.50_+1.7_0.55'
#ModelDescription = '0.50_+2.0_0.55'
#ModelDescription = '0.50_+2.3_0.55'


#ModelDescription = '0.50_+0.0_0.55'
#ModelDescription = '0.50_+2.3_0.55'

for i in range(NumParts):
    
    print('Doing part %d of %d'%(i+1,NumParts))
    
    #FileName = 'KELT9bModelSpectra_Part%dof%d.txt.npy'%(i+1,NumParts)
#    FileName = '%s_parts/KELT9b_%s_Part%dof%d.txt'%(name,name,i+1,NumParts)
    #FileName = 'ModelDescription_%s/%s_parts/KELT9b_%s_Part%dof%d.txt'%(ModelDescription,name,name,i+1,NumParts)
    
    ### A quick hack for FeII
    FileName = 'ModelDescription_%s/%s_parts/KELT9b_%s_Part%dof%d.txt'%(ModelDescription,name,'FeII',i+1,NumParts)

    #FileName = 'ModelDescription_%s/%s_parts_AbunTimes100/KELT9b_%s_Part%dof%d.txt'%(ModelDescription,name,name,i+1,NumParts)
    
    
    SpecPart = np.loadtxt(FileName)
    
    #plt.plot(SpecPart[:,0],SpecPart[:,1])   
    
    VStackOfSpec = np.vstack((VStackOfSpec,SpecPart))   
    






VStackOfSpec = VStackOfSpec[1:,:]
    
SortedVStackOfSpec = VStackOfSpec[VStackOfSpec[:,0].argsort()]



WToInterpolateOverlaps = define_wavegrid(np.min(SortedVStackOfSpec[:,0]),np.max(SortedVStackOfSpec[:,0]),1e6)

VStackNoOverlaps = np.empty((len(WToInterpolateOverlaps),3))   

EmissionInterpObject = interpolate.interp1d(SortedVStackOfSpec[:,0],SortedVStackOfSpec[:,1])
TransInterpObject = interpolate.interp1d(SortedVStackOfSpec[:,0],SortedVStackOfSpec[:,2])

VStackNoOverlaps[:,0] = WToInterpolateOverlaps
VStackNoOverlaps[:,1] = EmissionInterpObject(WToInterpolateOverlaps)
VStackNoOverlaps[:,2] = TransInterpObject(WToInterpolateOverlaps)

SortedVStackOfSpec = VStackNoOverlaps

stationary_pl_spec_logspaced_wave = SortedVStackOfSpec[:,0]*u.um

pl_spec_logspaced_wave = stationary_pl_spec_logspaced_wave*(1+(NettPlanetRV/c_kms))

pl_spec_flux_per_Hz = SortedVStackOfSpec[:,1]*(u.erg/u.s/(u.cm**2)/u.Hz)

pl_spec_flux_per_um = pl_spec_flux_per_Hz.to(u.erg/u.s/(u.cm**2)/u.um, equivalencies=u.spectral_density(pl_spec_logspaced_wave))


pl_pRT_perHz_ForOutput = np.empty((len(stationary_pl_spec_logspaced_wave),2))
pl_pRT_perHz_ForOutput[:,0] = stationary_pl_spec_logspaced_wave.value
pl_pRT_perHz_ForOutput[:,1] = pl_spec_flux_per_Hz.value
#np.save('../ModelSpectra/ModelDescription_%s/KELT9b_%s_pRT_flux_per_Hz.npy'%(ModelDescription, name),pl_pRT_perHz_ForOutput)
#np.save('../ModelSpectra/ModelDescription_%s/KELT9b_%s_AbunTimes100_pRT_flux_per_Hz.npy'%(ModelDescription, name),pl_pRT_perHz_ForOutput)

SaveLocation = '../ModelSpectra/ModelDescription_%s'%(ModelDescription)

if not os.path.exists(SaveLocation):
    os.makedirs(SaveLocation)

np.save('%s/KELT9b_%s_%s_pRT_flux_per_Hz.npy'%(SaveLocation, name,ModelDescription),pl_pRT_perHz_ForOutput)


plt.plot(stationary_pl_spec_logspaced_wave.value,pl_spec_flux_per_Hz.value)



### For integrating the PHOENIX model spectrum to get flux at each wavelength point 

# stepsize = 0.01*u.AA


# #cwlims = ((5200,17100))*u.angstrom
# cwlims = ((SortedVStackOfSpec[0,0]*(u.um.to(u.AA)),SortedVStackOfSpec[-1][0]*(u.um.to(u.AA))))*u.angstrom

# w_constant_stepsize = np.arange(cwlims[0].value,cwlims[-1].value,stepsize.value)*u.AA

# bb_p = BlackBody(temperature=4048.99*u.K)
# flux_den_bbp = bb_p(w_constant_stepsize.to(u.um))*(np.pi*u.sr)  ### To convert intensity to flux, need to multiply by pi





# flux_den_bbp_per_um = flux_den_bbp.to(u.erg/u.s/(u.cm**2)/u.um, equivalencies=u.spectral_density(w_constant_stepsize.to(u.um)))

# Integrated_bbp = IntegrateF_nu_Or_F_lambda(w_constant_stepsize.to(u.um).value, flux_den_bbp_per_um.value, UseSpectres=True)*(u.erg/u.s/(u.cm**2))



# pl_em_spec_const_w_step = spectres.spectres(w_constant_stepsize.to(u.um).value, 
#                                             pl_spec_logspaced_wave.value,
#                                             pl_spec_flux_per_um.value)*(u.erg/u.s/(u.cm**2)/u.um)


# pl_em_spec_const_w_step_perHz = spectres.spectres(w_constant_stepsize.to(u.um).value, 
#                                             pl_spec_logspaced_wave.value,
#                                             pl_spec_flux_per_Hz.value)*(u.erg/u.s/(u.cm**2)/u.Hz)


# pl_em_spec_const_w_step_perHz_ForOutput = np.empty((len(pl_em_spec_const_w_step_perHz),2))
# pl_em_spec_const_w_step_perHz_ForOutput[:,0] = w_constant_stepsize.to(u.um).value
# pl_em_spec_const_w_step_perHz_ForOutput[:,1] = pl_em_spec_const_w_step_perHz.value
# #np.save('pRT_flux_per_Hz_constant_WaveStepSize.npy',pl_em_spec_const_w_step_perHz_ForOutput)


# plt.figure()
# plt.title('Flux densities of pRT and BB models')
# plt.plot(stationary_pl_spec_logspaced_wave,pl_spec_flux_per_Hz,label='pRT')
# plt.plot(w_constant_stepsize.to(u.um),flux_den_bbp,label='black body')
# plt.xlabel(r'wavelength ($\mu$m)')
# plt.ylabel(r'flux density (erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')
# plt.legend()
# plt.tight_layout()
# plt.savefig('ModelPlots/PlanetFluxDensityPerHz.png',dpi=400)


# plt.figure()
# plt.title('Flux densities of pRT and BB models')
# plt.plot(w_constant_stepsize.to(u.um),pl_em_spec_const_w_step,label='pRT')
# plt.plot(w_constant_stepsize.to(u.um),flux_den_bbp_per_um,label='black body')
# plt.xlabel(r'wavelength ($\mu$m)')
# plt.ylabel(r'flux density (erg s$^{-1}$ cm$^{-2}$ $\mu$m$^{-1}$)')
# plt.legend()
# plt.tight_layout()
# plt.savefig('ModelPlots/PlanetFluxDensityPerMicron.png',dpi=400)


# Integrated_pl_em_spec_const_w_step = IntegrateF_nu_Or_F_lambda(w_constant_stepsize.to(u.um).value, \
#                                                                pl_em_spec_const_w_step.value, \
#                                                             UseSpectres=True)*(u.erg/u.s/(u.cm**2))

# IntegratedPlanetSpec = np.empty((len(Integrated_pl_em_spec_const_w_step),2))
# IntegratedPlanetSpec[:,0] = w_constant_stepsize.to(u.um).value
# IntegratedPlanetSpec[:,1] = Integrated_pl_em_spec_const_w_step.value

# #np.save('IntegratedFluxPlanetSpec.npy',IntegratedPlanetSpec)

# plt.figure()
# plt.title('Integrated flux of pRT and BB models')
# plt.plot(w_constant_stepsize.to(u.um),Integrated_pl_em_spec_const_w_step,label='pRT')
# plt.plot(w_constant_stepsize.to(u.um),Integrated_bbp,label='black body')
# plt.xlabel(r'wavelength ($\mu$m)')
# plt.ylabel(r'Flux (erg s$^{-1}$ cm$^{-2}$)')
# plt.legend()
# plt.tight_layout()
# plt.savefig('ModelPlots/PlanetIntegratedFlux.png',dpi=400)



# ####################################
# ## For the PHOENIX model 


# st_wave_padding_A = 500

# #st_w_constant_stepsize = np.arange(cwlims[0].value-st_wave_padding_A,cwlims[-1].value+st_wave_padding_A,stepsize.value)*u.AA
# st_w_constant_stepsize = np.copy(w_constant_stepsize)

# bb_s = BlackBody(temperature=10170*u.K)
# flux_den_bbs = bb_s(st_w_constant_stepsize.to(u.um))*(np.pi*u.sr) ### To convert intensity to flux, need to multiply by pi
# flux_den_bbs_per_um = flux_den_bbs.to(u.erg/u.s/(u.cm**2)/u.um, equivalencies=u.spectral_density(st_w_constant_stepsize.to(u.um)))

# bbs_ForOutput = np.empty((len(st_w_constant_stepsize),2))
# bbs_ForOutput[:,0] = st_w_constant_stepsize.to(u.um).value
# bbs_ForOutput[:,1] = flux_den_bbs.value
# #np.save('FluxDensityPerHzStellarBlackBodySpec.npy',bbs_ForOutput)

# Integrated_bbs = IntegrateF_nu_Or_F_lambda(st_w_constant_stepsize.to(u.um).value, flux_den_bbs_per_um.value, UseSpectres=True)*(u.erg/u.s/(u.cm**2))

# IntegratedBBSpec = np.empty((len(Integrated_bbs),2))
# IntegratedBBSpec[:,0] = st_w_constant_stepsize.to(u.um).value
# IntegratedBBSpec[:,1] = Integrated_bbs.value

# #np.save('IntegratedFluxStellarBlackBodySpec.npy',IntegratedBBSpec)

# FluxDensityPerHzStellarBB = np.empty((len(flux_den_bbs),2))
# IntegratedBBSpec[:,0] = st_w_constant_stepsize.to(u.um).value
# IntegratedBBSpec[:,1] = flux_den_bbs.value

# np.save('IntegratedFluxStellarBlackBodySpec.npy',IntegratedBBSpec)


# modelfileflux = pyfits.open('../PHOENIX_model_spectrum/lte10200-4.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')
# modelfilew = pyfits.open('../PHOENIX_model_spectrum/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')

# mf_full = modelfileflux[0].data*(u.erg/u.s/(u.cm**2)/u.cm)
# mw_full = modelfilew[0].data*u.angstrom 
    

# wantedIndices = np.where((mw_full>cwlims[0])&(mw_full<cwlims[-1]))

# mf = mf_full[wantedIndices]
# mw = mw_full[wantedIndices]*(1+(NettStellarRV/c_kms))

# mf_per_Hz = mf.to(u.erg/u.s/(u.cm**2)/u.Hz, equivalencies=u.spectral_density(mw.to(u.um)))

# plt.figure()
# plt.title('Flux densities of PHOENIX and BB models')
# plt.plot(mw.to(u.um),mf_per_Hz,label='PHOENIX')
# plt.plot(st_w_constant_stepsize.to(u.um),flux_den_bbs,label='stellar BB')
# plt.xlabel(r'wavelength ($\mu$m)')
# plt.ylabel(r'Flux (erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')
# plt.legend()
# plt.tight_layout()
# plt.savefig('ModelPlots/StellarFluxDensityPerHz.png',dpi=400)




# mf_const_step_size = spectres.spectres(st_w_constant_stepsize.to(u.cm).value, mw.to(u.cm).value, mf.value,fill=np.nan)*(u.erg/u.s/(u.cm**2)/u.cm)

# mf_const_step_size_perHZ = spectres.spectres(st_w_constant_stepsize.to(u.cm).value, mw.to(u.cm).value, mf_per_Hz.value,fill=np.nan)*(u.erg/u.s/(u.cm**2)/u.Hz)


# #print('Integrating the log spaced waves')
# ####Integrated_mf_logwspacing = IntegrateF_nu_Or_F_lambda(mw.to(u.cm).value, mf.value, UseSpectres=True)*(u.erg/u.s/(u.cm**2))

# print('Integrating the const spaced waves')
# Integrated_mf_const_step_size = IntegrateF_nu_Or_F_lambda(st_w_constant_stepsize.to(u.cm).value, mf_const_step_size.value, UseSpectres=True)*(u.erg/u.s/(u.cm**2))

# ###############


# IntegratedStarSpec = np.empty((len(Integrated_mf_const_step_size),2))
# IntegratedStarSpec[:,0] = st_w_constant_stepsize.to(u.um).value
# IntegratedStarSpec[:,1] = Integrated_mf_const_step_size.value

# np.save('IntegratedFluxStellarSpec.npy',IntegratedStarSpec)

# plt.figure()
# plt.title('Integrated flux of PHOENIX and stellar BB models')
# plt.plot(st_w_constant_stepsize.to(u.um),Integrated_mf_const_step_size,label='PHOENIX')
# plt.plot(st_w_constant_stepsize.to(u.um),Integrated_bbs,label='stellar black body')
# plt.xlabel(r'wavelength ($\mu$m)')
# plt.ylabel(r'flux density (erg s$^{-1}$ cm$^{-2}$)')
# plt.legend()
# plt.tight_layout()
# plt.savefig('ModelPlots/StellarIntegratedFlux.png',dpi=400)

# plt.figure()
# plt.title('Black bodies flux density (per Hz)')
# plt.plot(st_w_constant_stepsize.to(u.um),flux_den_bbs,label='stellar BB')
# plt.plot(w_constant_stepsize.to(u.um),flux_den_bbp,label='planet BB')
# plt.ylabel(r'Flux density (erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')
# plt.xlabel(r'wavelength ($\mu$m)')
# plt.tight_layout()
# plt.legend()
# plt.savefig('ModelPlots/BlackBodyFluxDensityPerHzComparison.png',dpi=400)

# plt.figure()
# plt.title(r'Black bodies flux density (per $\mu$m)')
# plt.plot(st_w_constant_stepsize.to(u.um),flux_den_bbs_per_um,label='stellar BB')
# plt.plot(w_constant_stepsize.to(u.um),flux_den_bbp_per_um,label='planet BB')
# plt.ylabel(r'Flux density (erg s$^{-1}$ cm$^{-2}$ $\mu$m$^{-1}$)')
# plt.xlabel(r'wavelength ($\mu$m)')
# plt.tight_layout()
# plt.legend()
# plt.savefig('ModelPlots/BlackBodyFluxDensityPerMicronComparison.png',dpi=400)

# plt.figure()
# plt.title(r'Black bodies integrated flux')
# plt.plot(st_w_constant_stepsize.to(u.um),Integrated_bbs,label='stellar BB')
# plt.plot(w_constant_stepsize.to(u.um),Integrated_bbp,label='planet BB')
# plt.ylabel(r'Flux (erg s$^{-1}$ cm$^{-2}$)')
# plt.xlabel(r'wavelength ($\mu$m)')
# plt.tight_layout()
# plt.legend()
# plt.savefig('ModelPlots/BlackBodyFluxes.png',dpi=400)



# plt.figure()
# plt.title('Integrated flux of pRT and PHOENIX')
# plt.plot(st_w_constant_stepsize.to(u.um),Integrated_mf_const_step_size,label='PHOENIX')
# plt.plot(w_constant_stepsize.to(u.um),Integrated_pl_em_spec_const_w_step,label='pRT')
# plt.xlabel(r'wavelength ($\mu$m)')
# plt.ylabel(r'Flux (erg s$^{-1}$ cm$^{-2}$)')
# plt.legend()
# plt.tight_layout()
# plt.savefig('ModelPlots/StellarIntegratedFlux.png',dpi=400)


# ##### calculate ratios with the same wavelength scales 



# FluxDensityPerHz_pRT_over_PHOENIX = pl_em_spec_const_w_step_perHz/mf_const_step_size_perHZ
# FluxDensityPerHz_pRT_over_BB = pl_em_spec_const_w_step_perHz/flux_den_bbs
# FluxDensityPerHz_plBB_over_stBB = flux_den_bbp/flux_den_bbs


# plt.figure()
# plt.title('Ratio of planet and stellar flux densities')
# plt.plot(st_w_constant_stepsize.to(u.um),FluxDensityPerHz_pRT_over_PHOENIX,label='pRT over PHOENIX')
# plt.plot(w_constant_stepsize.to(u.um),FluxDensityPerHz_pRT_over_BB,label='pRT over stellar BB')
# plt.plot(w_constant_stepsize.to(u.um),FluxDensityPerHz_plBB_over_stBB,label='planet BB over stellar BB')
# plt.xlabel(r'wavelength ($\mu$m)')
# plt.ylabel(r'Flux density ratio')
# plt.legend()
# plt.tight_layout()
# plt.savefig('ModelPlots/RatioOfStellarAndPlanetFluxeDensities.png',dpi=400)


# Flux_pRT_over_PHOENIX = Integrated_pl_em_spec_const_w_step/Integrated_mf_const_step_size
# Flux_pRT_over_BB = Integrated_pl_em_spec_const_w_step/Integrated_bbs
# Flux_plBB_over_stBB = Integrated_bbp/Integrated_bbs


# plt.figure()
# plt.title('Ratio of planet and stellar fluxes')
# plt.plot(st_w_constant_stepsize.to(u.um),Flux_pRT_over_PHOENIX,'.',label='pRT over PHOENIX')
# plt.plot(w_constant_stepsize.to(u.um),Flux_pRT_over_BB,'.',label='pRT over stellar BB')
# plt.plot(w_constant_stepsize.to(u.um),Flux_plBB_over_stBB,'.',label='planet BB over stellar BB')
# plt.xlabel(r'wavelength ($\mu$m)')
# plt.ylabel(r'Flux ratio')
# plt.legend()
# plt.tight_layout()
# plt.savefig('ModelPlots/RatioOfStellarAndPlanetFluxes.png',dpi=400)





endtime = time.time()
TotalTime = endtime - starttime

print('Time required %.2f mins'%(TotalTime/60))