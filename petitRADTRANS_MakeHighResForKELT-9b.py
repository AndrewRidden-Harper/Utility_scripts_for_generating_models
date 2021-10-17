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


def IntegrateF_nu_Or_F_lambda(wave_um,planet_flux_from_file):

    f = interpolate.interp1d(wave_um, planet_flux_from_file)    
    
    WaveBinLimits = np.zeros((len(wave_um)+1))
    
    WaveBinLimits[0] = wave_um[0]
    WaveBinLimits[-1] = wave_um[-1]
    
    MiddlePoints = (wave_um[0:-1] + wave_um[1:])/2    
    
    WaveBinLimits[1:-1] = MiddlePoints    
    
    FluxAtMiddlePoints = f(WaveBinLimits)
    
    Integration = np.empty_like(planet_flux_from_file)
    
    for i in range(len(wave_um)):
        
        Integration[i] = np.trapz(FluxAtMiddlePoints[i:i+2], WaveBinLimits[i:i+2])        
    
    ##IntegratedJayeshEmissionModel = Integration*planet_flux_from_file  ## First wrong attempt 
    
    return np.abs(Integration)  ## integrating over Hz is negative because it goes from high to low freq.


def ConvertJayeshEmissionFluxToMatchpRT(planet_flux_from_file):


    fnu_wvl = (planet_flux_from_file)*1E7 # Flux converting from W/(m^2)/um to erg/(cm^2)/s.cm
    
    c_vel = 3*(10**10) ## velocity of light in cm/s
    
    wv_num_grid = np.arange(0,50010,10,dtype='float') ## Wavenumber grid 
    wv_num_grid[0] = 0.01 # modifying the 1st wavenumber to avoid divide by zero error
    wvl_grid = 1/wv_num_grid
    freq_grid = c_vel/wvl_grid 
    
    wvl_bin_size = []
    freq_bin_size = []
    for i in range(len(freq_grid) - 1):
     wvl_bin_size.append(wvl_grid[i] - wvl_grid[i+1])
     freq_bin_size.append(freq_grid[i+1] - freq_grid[i])
    
    wvl_bin_size = wvl_bin_size[::-1]
    freq_bin_size = freq_bin_size[::-1]
    
    fnu_freq = (fnu_wvl[:]*wvl_bin_size[:])/(freq_bin_size[:]) ## Flux in erg/(cm^2)/s/Hz

    return fnu_freq



#filenameToLoad = 'LogSpacePressureMin-6Max3Steps130Species_H2_CO_CH4_CO2_C2H2_H_OH_He_NH3_HITRAN_HCN_Na_K_TiO_SiO_H2S_FeH_PH3_HITRAN_H-_VO_H2O.txt'

TotalRunTimeStart = time.time()

print('loading the mass fractions')
#filenameToLoad = 'LogSpacePressureMin-6Max3Steps130Species_H2_CO_CH4_CO2_C2H2_H_OH_He_NH3_HITRAN_HCN_Na_K_TiO_SiO_H2S_Fe_FeH_PH3_HITRAN_H-_e-_VO_H2O.txt'
#filenameToLoad = 'LogSpacePressureMin-6Max3Steps130Species_H2_CO_CH4_CO2_C2H2_H_OH_He_NH3_HITRAN_HCN_Na_K_TiO_SiO_H2S_Fe_FeH_PH3_HITRAN_H-_e-_VO_H2O.txt'

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


#SpeciesName = 'Fe'
#SpeciesName = 'Ti'
#SpeciesName = 'Ca'
#SpeciesName = 'Cr'
#SpeciesName = 'Mg'

#SpeciesName = 'Al'
#SpeciesName = 'K'
#SpeciesName = 'Na'

#SpeciesName = 'CaII'
#SpeciesName = 'CO_all_iso'

#SpeciesName = 'TiO_all_iso_Plez' 
#SpeciesName = 'TiO_48_Exomol_McKemmish' 
#SpeciesName = 'TiO_48_Plez' 

#SpeciesName = 'FeH_main_iso'

#SpeciesName = 'VO'

SpeciesName = 'FeII'


filenameToLoad = 'Interp_K9b_PTs/OriginalBeforeLoop/LogSpacePressureMin-6Max3Steps130_T_P_abundances_highresnames_%s.txt'%(ModelDescription)
#filenameToLoad = 'Interp_K9b_PTs/OriginalBeforeLoop/ExoMol48TiOLogSpacePressureMin-6Max3Steps130_T_P_abundances_highresnames_%s.txt'%(ModelDescription)
#filenameToLoad = 'TiO_48_Plez_LogSpacePressureMin-6Max3Steps130_T_P_abundances_highresnames_%s.txt'%(ModelDescription)



print('finished loading the mass fractions')

InterpolatedInputs_df = pd.read_csv(filenameToLoad)

# pressure (bar), temperature (K), H2 mass fraction, CO mass fraction, H mass fraction, He mass fraction, Mg mass fraction, Fe mass fraction, H_ mass fraction

#['H2','CO','H','He','Mg','Fe','H_']
                                                          
                                    
# atmosphere = Radtrans(line_species = ['CO','H','He','Mg','Fe'], \
#       rayleigh_species = ['H2', 'He'], \
#       continuum_opacities = ['H2-H2', 'H2-He','H-'], \
#       wlen_bords_micron = [0.2, 1])
      #wlen_bords_micron = [0.2, 19])



# atmosphere = Radtrans(line_species = ['H2','CO','CH4','CO2','C2H2','H','OH','He',
#                                       'NH3_HITRAN','HCN','Na','K','TiO','SiO',
#                                       'H2S','FeH','PH3_HITRAN','VO','H2O'], \
#       rayleigh_species = ['H2', 'He'], \
#       continuum_opacities = ['H2-H2', 'H2-He','H-'], \
#       wlen_bords_micron = [0.2, 50])   ## 0.2, 15 works well 

StartTime = time.time()



EdgePadding = 0.05
#WavelengthRange = np.array([0.520-EdgePadding,1.710+EdgePadding]) ## CARMENES wavelength range 

#WavelengthRange = np.array([0.520-EdgePadding,0.595+EdgePadding])
#WavelengthRange = np.array([0.52-EdgePadding,0.52+0.238+EdgePadding])

#### Covers the entire CARMENES range and splits it into enough intervals 
FullWavelengthRange = np.array([0.595-EdgePadding,1.710+EdgePadding])
NIntervals = 20

# FullWavelengthRange = np.array([0.595,0.615])
# NIntervals = 1


IntervalSize = (FullWavelengthRange[1]-FullWavelengthRange[0])/NIntervals

print('Interval size: %.2f'%(IntervalSize))

for i in range(NIntervals):
    
    LoopStartTime = time.time()
    
    print('Doing part %d of %d'%((i+1),NIntervals))
    
    InputW1 = FullWavelengthRange[0]+i*IntervalSize - EdgePadding
    InputW2 = FullWavelengthRange[0]+(i+1)*IntervalSize + EdgePadding
    
    print('InputW1',InputW1)
    print('InputW2',InputW2)
    print('InputW range',InputW2-InputW1)
 
    print('Starting to load opacities')
    
    ### Radtrans object that works for low res    
    # atmosphere = Radtrans(line_species = ['CO','CH4','CO2','C2H2','OH','NH3_HITRAN','HCN','Na','K','TiO','SiO_main_iso','H2S','FeH','PH3_HITRAN','VO','H2O'], \
    #       rayleigh_species = ['H2', 'He'], \
    #       continuum_opacities = ['H2-H2', 'H2-He','H-'], \
    #       wlen_bords_micron = [0.2, 15])  # 0.3, 15 works well 
        
    ### Example for high res spectrum     
    # atmosphere = Radtrans(line_species = ['CO_all_iso','Fe','Mg','Ca','CaII',
    #       'H2O_main_iso','SiO_main_iso','Na', 'K','Al','Ti'], \
    #       rayleigh_species = ['H2', 'He','N2'], \
    #       continuum_opacities = ['H2-H2', 'H2-He','H-'], \
    #       wlen_bords_micron = [WavelengthRange[0], WavelengthRange[1]], \
    #       mode = 'lbl')
    
    
    
    #atmosphere = Radtrans(line_species = ['CaII'], \
    atmosphere = Radtrans(line_species = [SpeciesName], \
          rayleigh_species = ['H2', 'He','N2'], \
          continuum_opacities = ['H2-H2', 'H2-He','H-'], \
          wlen_bords_micron = [InputW1, InputW2], \
          mode = 'lbl')
        
    # print('finished loading opacities')
        
    # #Now, let’s define the pressures of the atmospheric layers. Note that the pressures must always be sorted in increasing order, and be equidistant in log-space. The pressure is in units of bar, although all other units in petitRADTRANS are in cgs. Typically include of the order of 100 layers in your computations:
        
    pressures = np.logspace(-6, 3, 130)
    
    # #Units in petitRADTRANS: all units in petitRADTRANS are in cgs, except for pressure, which is in bars, and the mean molecular weight (MMW), which is in units of atomic mass units.
    
    atmosphere.setup_opa_structure(pressures)
    
    # # ### Calculating a transmission spectrum
    
    temperature = InterpolatedInputs_df['temperature (K)'].to_numpy()
    
    mass_fractions = {}
    
    ### Mass fractions for low res example 
    mass_fractions['H2'] = InterpolatedInputs_df['H2_main_iso mass fraction'].to_numpy()
    #mass_fractions['H2_main_iso'] = InterpolatedInputs_df['H2_main_iso mass fraction'].to_numpy()
    mass_fractions['N2'] = InterpolatedInputs_df['N2 mass fraction'].to_numpy()
    #mass_fractions['CO_all_iso'] = InterpolatedInputs_df['CO_all_iso mass fraction'].to_numpy()
    # mass_fractions['CH4'] = InterpolatedInputs_df['CH4 mass fraction'].to_numpy()
    # mass_fractions['CO2'] = InterpolatedInputs_df['CO2 mass fraction'].to_numpy()
    # mass_fractions['C2H2'] = InterpolatedInputs_df['C2H2 mass fraction'].to_numpy()
    mass_fractions['H'] = InterpolatedInputs_df['H mass fraction'].to_numpy()
    # mass_fractions['OH'] = InterpolatedInputs_df['OH mass fraction'].to_numpy()
    mass_fractions['He'] = InterpolatedInputs_df['He mass fraction'].to_numpy()
    # mass_fractions['NH3_HITRAN'] = InterpolatedInputs_df['NH3_HITRAN mass fraction'].to_numpy()
    # mass_fractions['HCN'] = InterpolatedInputs_df['HCN mass fraction'].to_numpy()
    #mass_fractions['Na'] = InterpolatedInputs_df['Na mass fraction'].to_numpy()
    #mass_fractions['NaII'] = InterpolatedInputs_df['NaII mass fraction'].to_numpy()
    #mass_fractions['K'] = InterpolatedInputs_df['K mass fraction'].to_numpy()
    # mass_fractions['TiO'] = InterpolatedInputs_df['TiO mass fraction'].to_numpy()
    #mass_fractions['SiO_main_iso'] = InterpolatedInputs_df['SiO_main_iso mass fraction'].to_numpy()
    # mass_fractions['H2S'] = InterpolatedInputs_df['H2S mass fraction'].to_numpy()
    # mass_fractions['FeH'] = InterpolatedInputs_df['FeH mass fraction'].to_numpy()
    # mass_fractions['PH3_HITRAN'] = InterpolatedInputs_df['PH3_HITRAN mass fraction'].to_numpy()
    # mass_fractions['VO'] = InterpolatedInputs_df['VO mass fraction'].to_numpy()
    #mass_fractions['H2O_main_iso'] = InterpolatedInputs_df['H2O_main_iso mass fraction'].to_numpy()
    #mass_fractions['FeII'] = InterpolatedInputs_df['Fe mass fraction'].to_numpy()
    #mass_fractions['Fe'] = InterpolatedInputs_df['Fe mass fraction'].to_numpy()
    #mass_fractions['Mg'] = InterpolatedInputs_df['Mg mass fraction'].to_numpy()
    #mass_fractions['Ca'] = InterpolatedInputs_df['Ca mass fraction'].to_numpy()
    #mass_fractions['Al'] = InterpolatedInputs_df['Al mass fraction'].to_numpy()
    #mass_fractions['AlII'] = InterpolatedInputs_df['Al mass fraction'].to_numpy()
    #mass_fractions['Ti'] = InterpolatedInputs_df['Ti mass fraction'].to_numpy()
    #mass_fractions['CaII'] = InterpolatedInputs_df['CaII mass fraction'].to_numpy()
    #mass_fractions['Cr'] = InterpolatedInputs_df['Cr mass fraction'].to_numpy()
    mass_fractions['H-'] = InterpolatedInputs_df['H- mass fraction'].to_numpy()
    mass_fractions['e-'] = InterpolatedInputs_df['e- mass fraction'].to_numpy()
    mass_fractions[SpeciesName] = InterpolatedInputs_df['%s mass fraction'%(SpeciesName)].to_numpy()
    
    ### A hack to use Fe abundances for Fe II
    #mass_fractions['FeII'] = InterpolatedInputs_df['Fe mass fraction'].to_numpy()
    
    

    
    #np.savetxt('DumpMet2_FeAbunTimes100',InterpolatedInputs_df['Fe mass fraction'].to_numpy()*100)
    
    
    # # # # #Abundances in petitRADTRANS: abundances in pRT are in units of mass fractions, not number fractions (aka volume mixing ratio, VMR). You can convert between mass fractions and VMRs by using
    
    MMW = InterpolatedInputs_df['mixture mean molar mass'].to_numpy()
    
    ########################
    ### To make the plot of TiO abundance (in VMR)
    print()
    print(SpeciesName)
    #VMRTiO = mass_fractions[SpeciesName]*MMW/63.866
    VMRTiO = mass_fractions[SpeciesName]*MMW/22.989769
    

       
    lines = plt.loglog(pressures, VMRTiO,'.-')
      

    
    plt.ylabel('Abundance (VMR)')
    plt.xlabel('Pressure (bar)')
    ##xlim(1e-15,1)
    #tick_params(axis='x',labelsize='20')
    #tick_params(axis='y',labelsize='20')
    
    #ax = gca()
    #title('Chemical Abundances',fontsize=15,fontweight='bold')
    #legend(loc=0,fontsize=15).draw_frame(0)
    
    plt.gca().invert_yaxis()
    
    # plt.title('Goyal TiO VMR')
    # plt.savefig('Goyal_TiO_VMR.png',dpi=400)
    
    plt.title('Goyal Na VMR')
    plt.savefig('Goyal_Na_VMR.png',dpi=400)
    
    #raise Exception
    ### End making plot of TiO abundance (in VMR)
    
    R_pl = 1.891*nc.r_jup_mean
    gravity = 1e1**3.3002
    P0 = 0.01
    
    # # ## To generate transmission spectrum 
    atmosphere.calc_transm(temperature, mass_fractions, gravity, MMW, R_pl=R_pl, P0_bar=P0)
    
    TransitRadius_Rjup = atmosphere.transm_rad/nc.r_jup_mean
    
    
    
    # # ## To generate emission spectrum 
    atmosphere.calc_flux(temperature, mass_fractions, gravity, MMW)
    
    LoopEndTime = time.time()
    
    #ModelWaveMicron = (nc.c/atmosphere.freq/1e-4)*u.micron
    
    #ModelFreqHz = ModelWaveMicron.to(u.Hz, equivalencies=u.spectral()) 
    
    
    #ModelFluxPerHz = atmosphere.flux*(u.erg/((u.cm**2)*u.s*u.Hz))
    
    #ModelFlux = ModelFluxPerHz*ModelFreqHz
    
    #ModelFluxWm2 = ModelFlux.to(u.Watt/(u.m**2))
    
    #ModelFluxWperm2PerMircon = ModelFluxWm2/ModelWaveMicron
    
    
    pRTwave_um = nc.c/atmosphere.freq/1e-4
    
    
    # pRTspec = np.empty((len(pRTwave_um),2))
    # pRTspec[:,0] = pRTwave_um
    # pRTspec[:,1] = atmosphere.flux
    # np.savetxt('pRTspec.txt',pRTspec)
    
    spec_output = np.empty((len(pRTwave_um),3))
    spec_output[:,0] = pRTwave_um
    spec_output[:,1] = atmosphere.flux 
    spec_output[:,2] = TransitRadius_Rjup
    
    
    header = 'wavelength (micron),emission flux (erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}),transit radius (Rjup)'

    SaveLocation = 'ModelDescription_%s/%s_parts'%(ModelDescription,SpeciesName)
    #SaveLocation = 'ModelDescription_%s/%s_parts'%(ModelDescription,'FeII_UsingFeI')
    if not os.path.exists(SaveLocation):
        os.makedirs(SaveLocation)

    np.savetxt('%s/KELT9b_%s_Part%dof%d.txt'%(SaveLocation,SpeciesName,i+1,NIntervals),spec_output,header=header)
    #np.savetxt('KELT9bModelSpectra_Part%dof%d.txt'%((i+1),NIntervals),spec_output,header=header)
    ##np.save('KELT9bModelSpectra_Part%dof%d.txt'%((i+1),NIntervals),spec_output)

    
    plt.figure()
    
    plt.plot(pRTwave_um, TransitRadius_Rjup)
    
    plt.xlabel('Wavelength (microns)')
    plt.ylabel(r'Transit radius ($\rm R_{Jup}$)')
    
    
    plt.figure()
    plt.plot(pRTwave_um, atmosphere.flux/1e-6,label='petitRADTRANS emission spectrum\nbased on Jayesh\'s model atmosphere')
    plt.ylabel(r'Planet flux $F_\nu$ (10$^{-6}$ erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$)')
    plt.xlabel('Wavelength (microns)')
    plt.title('KELT-9b high res')
    
    print('Time taken to do %.2f um (%.2f to %.2f): %.2f mins'%(InputW2-InputW1,InputW1,InputW2,(LoopEndTime-LoopStartTime)/60))

TotalRunTimeEnd = time.time()
print('Total run time: %f mins'%((TotalRunTimeEnd-TotalRunTimeStart)/60))

# plt.legend()
# plt.savefig('EmissionSpectraComparison_%.2fto%.2fum.png'%(pRTwave_um[0],pRTwave_um[-1]),dpi=400)

