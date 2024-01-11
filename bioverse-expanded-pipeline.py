import glob
from bioverse.util import compute_t_ref
from bioverse.survey import TransitSurvey
from bioverse.generator import Generator
from atmos_to_PSG import atmosatm, psgspec
import pyatmos
import os
from tqdm import tqdm, trange
import numpy as np
import matplotlib.pyplot as plt

atmos = pyatmos.Simulation(code_path='~/atmos-master', DEBUG=True)
atmos.start()

curr_dir=os.getcwd()

O2_flux_list=[5E+10,9E+10,3.94985E+11,8E+11]


for i in range(len(O2_flux_list)):
    O2_flux=O2_flux_list[i]
    output_dir=curr_dir+'/sample_atmos_results/O2_test_{n}'.format(n=i)
    args = {
        'species_fluxes': {'O2' :O2_flux},
        'max_photochem_iterations' : 50000, 
        'max_clima_steps' : 10, 
        'output_directory' : output_dir}
    
    
    atmos.run(**args)

model_list=glob.glob('sample_atmos_results/O2_test*')

tel='JWST'


# generate sample of observed planets

generator=Generator('transit')
survey=TransitSurvey('default')
sample, detected, data = survey.quickrun(generator, t_total=10*365.25)
T_st_ref = survey.T_st_ref
R_st_ref = survey.R_st_ref
survey.measurements['has_H2O'].t_ref=7.5



# generate spectra and then calculate t_ref

for i in range(len(model_list)):
    model=model_list[i].split('/')[1]

    #create the null spectrum for t_ref
    
    null_newf=atmosatm(model_list[i],tel=tel,filebase=model,null_spec=True,removed_gas="O2",star='M')
    psgspec(model,null_newf,showplot=False,null_spec=True)
    null_rad=curr_dir+'/psg_output/%s_null_rad.txt' % model

    #create radiance spectrum

    newf=atmosatm(model_list[i],tel=tel,filebase=model,star='M')
    psgspec(model,newf,showplot=False)
    model_rad=curr_dir+'/psg_output/%s_rad.txt' % model
    t_ref = compute_t_ref(filenames=(model_rad,null_rad), t_exp=10, wl_min=0.4, wl_max=0.9,)
    print("Required exposure time for %s: {:.1f} hr".format(t_ref) % model)
    survey.measurements['has_O2'].t_ref=t_ref/24
    survey.save()
    
    
    t_exp, N_obs = survey.measurements['has_O2'].compute_exposure_time(data[detected['EEC']])

    fig, ax = plt.subplots(ncols=2, figsize=(16,8))

    bins = np.logspace(np.log10(0.01), np.log10(np.amax(t_exp)), 30)
    ax[0].hist(t_exp, bins=bins)
    ax[0].set_xscale('log')
    ax[0].set_xlabel('Exposure time (d)')
    ax[0].set_ylabel('Number of EECs')
    ax[0].axvline(1000/24, linestyle='dashed', lw=5, c='black')

    bins = np.logspace(0, 5, 30)
    ax[1].hist(N_obs, bins=bins)
    ax[1].set_xscale('log')
    ax[1].set_xlabel('Number of transit observations')
    ax[1].axvline(survey.N_obs_max, linestyle='dashed', lw=5, c='black')

    plt.subplots_adjust(wspace=0.3)

    plt.show()

    
    
