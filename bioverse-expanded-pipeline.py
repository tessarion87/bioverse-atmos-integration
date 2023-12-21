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

O2_flux_list=[8.7E+11,1.45E+12,9.64E+13]


# for i in range(len(O2_flux_list)):
    # O2_flux=O2_flux_list[i]
    # output_dir=curr_dir+'/sample_atmos_results/O2_test_{n}'.format(n=i)
    # args = {
        # 'species_fluxes': {'O2' :O2_flux},
        # 'max_photochem_iterations' : 50000, 
        # 'max_clima_steps' : 10, 
        # 'output_directory' : output_dir}
    
    
    # atmos.run(**args)

model_list=glob.glob('sample_atmos_results/O2_test*')

tel='JWST'


# generate sample of observed planets

generator=Generator('transit')
survey=TransitSurvey('default')
sample, detected, data = survey.quickrun(generator, t_total=10*365.25)

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
    
    bins = np.linspace(0, 150, 15)
    obs= ~np.isnan(detected['has_O2'])
    EEC = detected['EEC']
    plt.hist(data['d'][obs&EEC], density=True, histtype='step', lw=2, bins=bins, label='Observed')
    plt.hist(data['d'][~obs&EEC], density=True, histtype='step', lw=2, bins=bins, label='Not observed')
    plt.show()
