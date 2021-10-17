# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 02:02:43 2021

@author: ar-h1
"""


import numpy as np
#import matplotlib.pyplot as plt 
#import matplotlib.pyplot as plt 
from eniric import broaden 
import time
import os 

#d = np.load('../ModelSpectra/ModelDescription_0.50_+0.0_0.55/KELT9b_Fe_0.50_+0.0_0.55_pRT_flux_per_Hz.npy')

#name = 'Fe_0.50_+0.0_0.55'
#name = 'Al_0.50_+0.0_0.55'
#name = 'Al_0.50_+2.3_0.55'
#name = 'CaII_0.50_+0.0_0.55'
#name = 'CO_all_iso_0.50_+0.0_0.55'
#name = 'Cr_0.50_+0.0_0.55'
#name = 'FeH_main_iso_0.50_+0.0_0.55'
#name = 'K_0.50_+0.0_0.55'
#name = 'Mg_0.50_+0.0_0.55'
#name = 'Na_0.50_+0.0_0.55'
#name = 'Ti_0.50_+0.0_0.55'
#name = 'TiO_48_Exomol_McKemmish_0.50_+0.0_0.55'
#name = 'TiO_48_Plez_0.50_+0.0_0.55'

#name = 'TiO_all_iso_Plez_0.50_+0.0_0.55'

#name = 'VO_0.50_+0.0_0.55'

#name = 'FeII_0.50_+0.0_0.55'

#name = 'Ca_0.50_+0.0_0.55'

#name = 'Fe'

name = 'FeII_UsingFeI_0.50_+0.0_0.55'

ModelDescription = '0.50_+0.0_0.55'
#ModelDescription = '0.50_+2.3_0.55'

d = np.load('F:/KELT-9b_CARMENES_emission_data/ModelSpectra/ModelDescription_%s/KELT9b_%s_pRT_flux_per_Hz.npy'%(ModelDescription,name))

Savepath = 'F:/KELT-9b_CARMENES_emission_data/ModelSpectra/ModelDescription_%s/RotationalBroadening'%(ModelDescription)

# Savepath = 'F:/KELT-9b_CARMENES_emission_data/ScriptsForPetitRADTRANS/ModelsFromLoop/Fe'


# d = np.load('F:/KELT-9b_CARMENES_emission_data/ScriptsForPetitRADTRANS/ModelsFromLoop/%s/R1e6/K9b_%s_%s.npy'%(name,name,ModelDescription))



if not os.path.exists(Savepath):
    os.makedirs(Savepath)

# NumPointsToDo = 100000
# ext_w_pass_in = d[0:NumPointsToDo,0]
# ext_flux_pass_in = d[0:NumPointsToDo,1]

ext_w_pass_in = d[:,0]
ext_flux_pass_in = d[:,1]


output_w =  ext_w_pass_in[5:-5]

RotBroadened = np.zeros((len(output_w),2))

# model = d  

# phw = model[:,0]
# phf = model[:,1]
# c = 299792.458
# ModelInitialWavesVelocitySampling = c*np.diff(phw)/phw[1:]

# R = phw[1:]/np.diff(phw)


VsiniToPass = 6.63
epsilonToPass = 0.0 ## for no limb darkening 

RotStartTime = time.time()
rot = broaden.rotational_convolution(wavelength=output_w,extended_wav=ext_w_pass_in, extended_flux = ext_flux_pass_in,vsini=VsiniToPass,epsilon=epsilonToPass)
RotEndTime = time.time()

TotalTime = RotEndTime-RotStartTime

#print('%d, %f'%(NumPointsToDo,TotalTime))
print('%d, %f'%(len(ext_w_pass_in),TotalTime))

RotBroadened[:,0] = output_w
RotBroadened[:,1] = rot

np.save('%s/KELT9b_%s_Vrot%.2f_pRT_flux_per_Hz.npy'%(Savepath,name,VsiniToPass),RotBroadened)