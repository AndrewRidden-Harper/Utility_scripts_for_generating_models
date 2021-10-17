# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from petitRADTRANS import Radtrans
from petitRADTRANS import nat_cst as nc
import matplotlib.pyplot as plt 
import pandas as pd 
from astropy import units as u
from scipy import interpolate
import time
import os
from poor_mans_nonequ_chem import poor_mans_nonequ_chem as pm
from matplotlib.backends.backend_pdf import PdfPages





# atmosphere = Radtrans(line_species = [SpeciesOfInterest], \
#       rayleigh_species = ['H2', 'He'], \
#       continuum_opacities = ['H2-H2', 'H2-He'], \
#       wlen_bords_micron = [0.91, 2.48], \
#       mode = 'lbl')


#NumParts = 20
NumParts = 1
#FullWRange = 2.48 - 0.91
FullWRange = 1.07 - 0.66



 
 
 
WRangePerPart = FullWRange/NumParts

StaringWave = 0.66
##StaringWave = 0.9

BaseWaveLims = np.array([StaringWave,StaringWave+WRangePerPart])

WaveLimArray = np.zeros((NumParts,2))

WaveLimArray[0,:] = BaseWaveLims

print('BaseWaveLims')
print(BaseWaveLims)

for i in range(1,NumParts):
    
    WaveLimsPerPart = BaseWaveLims+(i)*WRangePerPart
    
    WaveLimArray[i,:] = WaveLimsPerPart
    
    print('WaveLimsPerPart %d: %s'%(i,WaveLimsPerPart))



# LineSpeciesList = ['CO_main_iso',
#                     'H2O_main_iso',
#                     'HCN_main_iso',
#                     'C2H2_main_iso',
#                     'CH4_main_iso',
#                     'PH3_main_iso',
#                     'CO2_main_iso',
#                     'NH3_main_iso',
#                     'H2S_main_iso',
#                     'VO',
#                     'TiO_48_Exomol_McKemmish',
#                     'Na',
#                     'K',
#                     'SiO_main_iso',
#                     'FeH_main_iso']
    
    
#CloudDeckPressure = 1.7e-5
CloudDeckPressure = None


#VMR_FeH = 10**(-5.25)
VMR_FeH = 10**(-7.56)
    
if CloudDeckPressure != None:
    SavePath = 'CloudP0_%.2e'%(CloudDeckPressure)
    SaveName = 'FeH_VMR%.3e_Cloud%.2e'%(VMR_FeH,CloudDeckPressure)
    
if CloudDeckPressure == None:    
    SavePath = 'NoClouds'
    SaveName = 'FeH_VMR%.3e_NoClouds'%(VMR_FeH)
    
    
if not os.path.exists(SavePath):
    os.makedirs(SavePath)
    

    
#LineSpeciesList = ['TiO_48_Exomol_McKemmish']
#LineSpeciesList = ['VO']
LineSpeciesList = ['FeH_main_iso']

pdf = PdfPages('%s/%s.pdf'%(SavePath,SaveName))

print('Constructing atmosphere')



    

for Part in range(NumParts):
    
    print('Doing part %d of %d'%(Part+1,NumParts))

    atmosphere = Radtrans(line_species = LineSpeciesList, \
          rayleigh_species = ['H2', 'He'], \
          continuum_opacities = ['H2-H2', 'H2-He'], \
          wlen_bords_micron = [WaveLimArray[Part,0], WaveLimArray[Part,1]], \
          mode = 'lbl')
    
    ### wlen_bords_micron = [0.91, 2.48], \
        
        
    ### wlen_bords_micron = [0.91, 1.7], \      
  
    
    ##Part = 2
    ### wlen_bords_micron = [1.7, 2.48], \
    
    
    print('Finished constructing atmosphere')
    
    #2.2, 2.4
    
    #wlen_bords_micron = [1.49, 2.46], \
        
    pressures = np.logspace(-10, 2, 130)
    atmosphere.setup_opa_structure(pressures)
    
    
    #####################3
    
    
    #R_pl = 1.838*nc.r_jup_mean
    #R_pl = 0.11643*nc.r_jup_mean
    R_pl = 1.3*nc.r_jup_mean

    
    #gravity = 1e1**2.45
    #gravity = 1640
    gravity = 1e1**2.4
    
    P0 = 0.01
    

    
    # kappa_IR = 0.01
    # gamma = 0.4
    # T_int = 200.
    # T_equ = 1500.
    
    #T_int = 90
    T_equ = 1400.0
    
    COs = 0.55 * np.ones_like(pressures)
    FeHs = 0. * np.ones_like(pressures)
    
    
    #temperature = nc.guillot_global(pressures, kappa_IR, gamma, gravity, T_int, T_equ)
    temperature = np.ones_like(pressures)*T_equ

    
    print('Getting chem eql abundances')
    
    mass_fractionsFromPM = pm.interpol_abundances(COs, \
                FeHs, \
                temperature, \
                pressures)
    
        
    print('Finished getting chem eql abundances')
        
    ###########################
    
    
    
    
    
    mass_fractions = {}
    mass_fractions['H2'] = mass_fractionsFromPM['H2']
    mass_fractions['He'] = mass_fractionsFromPM['He']
    
    #mass_fractions['CO_main_iso'] = mass_fractionsFromPM['CO']
    #mass_fractions['H2O_main_iso'] = mass_fractionsFromPM['H2O']
    
    #mass_fractions['HCN_main_iso'] = mass_fractionsFromPM['HCN']
    #mass_fractions['C2H2_main_iso'] = mass_fractionsFromPM['C2H2,acetylene']
    
    
    #mass_fractions['PH3_main_iso'] = mass_fractionsFromPM['PH3']
    #mass_fractions['CO2_main_iso'] = mass_fractionsFromPM['CO2']
    
    #mass_fractions['NH3_main_iso'] = mass_fractionsFromPM['NH3']
    
    
    #mass_fractions['H2S_main_iso'] = mass_fractionsFromPM['H2S']
    
    
    #mass_fractions['VO'] = mass_fractionsFromPM['VO']
    #mass_fractions['TiO_48_Exomol_McKemmish'] = mass_fractionsFromPM['TiO']
    
    
    # mass_fractions['Na'] = mass_fractionsFromPM['Na']
    # mass_fractions['K'] = mass_fractionsFromPM['K']
    
    
    # mass_fractions['SiO_main_iso'] = mass_fractionsFromPM['SiO']
    #mass_fractions['FeH_main_iso'] = mass_fractionsFromPM['FeH']
    
    
    
    
    
    
    
    ### MMW = mass_fractionsFromPM['MMW']
    MMW = 2.34*np.ones_like(pressures)
    
    
    ## VMR = mass_fractionsFromPM['CH4']*MMW/16.04
    #VMR = mass_fractionsFromPM['TiO']*MMW/63.866
    #VMR = mass_fractionsFromPM['VO']*MMW/66.9409
    #VMR = mass_fractionsFromPM['FeH']*MMW/56.853
    

    
    MolecWeight_FeH = 56.853
    
    mass_frac_FeH = VMR_FeH*MolecWeight_FeH/MMW
    
    mass_fractions['FeH_main_iso'] = mass_frac_FeH



    
    # RyanCH4VMR = 3e-3
    
    # RyanCH4MassFrac = RyanCH4VMR*16.04/MMW
    
    #mass_fractionsFromPM['CH4'] = RyanCH4MassFrac
    
    #mass_fractions['CH4_main_iso'] = RyanCH4MassFrac

    
    # with PdfPages('Andrew_Gliese486b_CH4.pdf') as pdf:
    #     plt.figure()
    #     plt.semilogy(temperature,pressures)
    #     plt.gca().invert_yaxis()
    #     plt.xticks(fontsize=20)
    #     plt.yticks(fontsize=20)
    #     plt.ylabel('Pressure (bar)',fontsize=20)
    #     plt.xlabel('Temperature (K)',fontsize=20)
    #     plt.title('Pressure-Temperature profile (Guillot)',fontsize=20)    
    #     pdf.savefig()
    
    #     plt.figure()
    #     plt.loglog(VMR,pressures)
    #     plt.gca().invert_yaxis()
    #     plt.ylabel('Pressure (bar)',fontsize=20)
    #     plt.xlabel('CH4 volume mixing ratio',fontsize=20)
    #     plt.xticks(fontsize=20)
    #     plt.yticks(fontsize=20)
    #     pdf.savefig()
        
    #     plt.figure()
    #     plt.loglog(mass_fractionsFromPM['CH4'],pressures)
    #     plt.gca().invert_yaxis()
    #     plt.ylabel('Pressure (bar)',fontsize=20)
    #     plt.xlabel('CH4 mass fraction',fontsize=20)
    #     plt.xticks(fontsize=20)
    #     plt.yticks(fontsize=20)
    #     pdf.savefig()
        
    #     plt.figure()
    #     plt.semilogy(MMW,pressures)
    #     plt.gca().invert_yaxis()
    #     plt.xticks(fontsize=20)
    #     plt.yticks(fontsize=20)
    #     plt.ylabel('Pressure (bar)',fontsize=20)
    #     plt.xlabel('MMW (amu)',fontsize=20)
    #     pdf.savefig()
        

    
    
    
    ##########################################
    
    if CloudDeckPressure != None:
        atmosphere.calc_transm(temperature, mass_fractions, gravity, MMW, R_pl=R_pl, P0_bar=P0,Pcloud=CloudDeckPressure)
    
    if CloudDeckPressure == None:
        atmosphere.calc_transm(temperature, mass_fractions, gravity, MMW, R_pl=R_pl, P0_bar=P0)

    
    plt.rcParams['figure.figsize'] = (10, 6)
    
    plt.plot(nc.c/atmosphere.freq/1e-4, atmosphere.transm_rad/nc.r_jup_mean)
    
    plt.xlabel('Wavelength (microns)')
    plt.ylabel(r'Transit radius ($\rm R_{Jup}$)')
    #plt.show()
    pdf.savefig()
    #plt.clf()
    
    
    #################
    
    

    
    #atmosphere.calc_flux(temperature, mass_fractions, gravity, MMW)
    
    # plt.figure()
    # plt.plot(nc.c/atmosphere.freq/1e-4, atmosphere.flux/1e-6)
    
    # #plt.xscale('log')
    # plt.xlabel('Wavelength (microns)')
    # plt.ylabel(r'Planet flux $F_\nu$ (10$^{-6}$ erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$)')
    # pdf.savefig()
    # plt.show()
    # plt.clf()
    
    results = np.zeros((len(nc.c/atmosphere.freq/1e-4),3))
    
    results[:,0] = nc.c/atmosphere.freq/1e-4
    results[:,1] = atmosphere.transm_rad/nc.r_jup_mean
    #results[:,2] = atmosphere.flux
    
    
    
    ################
    
    if not os.path.exists(SavePath):
        os.makedirs(SavePath)
    
    ##np.save('%s/Gl486b_ChemEql_Only%s.npy'%(SavePath,LineSpeciesList[0]),results)
    #np.save('%s/Gl486b_CH4_VMR_3e-3_New.npy'%(SavePath),results)
        
    np.save('%s/%s.npy'%(SavePath,SaveName),results)

plt.figure()
plt.plot(temperature,pressures)
plt.gca().invert_yaxis()
plt.yscale('log')
plt.ylabel('Pressure (bar)')
plt.xlabel('Temperature (K)')
pdf.savefig()

plt.figure()
plt.plot(MMW,pressures)
plt.gca().invert_yaxis()
plt.yscale('log')
plt.ylabel('Pressure (bar)')
plt.xlabel('MMW')
pdf.savefig()

plt.figure()
plt.plot(mass_frac_FeH*1e6,pressures)
plt.gca().invert_yaxis()
plt.yscale('log')
#plt.xscale('log')
plt.ylabel('Pressure (bar)')
plt.xlabel('%s mass fraction (ppm)'%(LineSpeciesList[0]))
pdf.savefig()
pdf.close()

    
    
    
    


