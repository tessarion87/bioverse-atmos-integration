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

atmos = pyatmos.Simulation(code_path='~/atmos-master',docker_image=None, DEBUG=True)
atmos.start()

curr_dir=os.getcwd()

H2_flux_list=[1E-3,1E-2,1E-1,1,1E2,1E4,1E6,1E8,1E10,1E12,1E13,1E14]


for i in range(len(H2_flux_list)):
    H2_flux=H2_flux_list[i]
    output_dir=curr_dir+'/sample_atmos_results/H2_test_{n}'.format(n=i)
    args = {
        'species_fluxes': {'H2' :H2_flux},
        'max_photochem_iterations' : 50000, 
        'max_clima_steps' : 10, 
        'output_directory' : output_dir}
    
    
    atmos.run(**args)

model_list=[]
for i in range(len(H2_flux_list)):
    model_name="sample_atmos_results/H2_test_"+str(i)
    model_list.append(model_name)


# generate sample of observed planets

generator=Generator('transit')
survey=TransitSurvey('default')
sample, detected, data = survey.quickrun(generator, t_total=10*365.25)

# generate spectra and then calculate t_ref
t_ref_list=[]
tel='JWST'

for i in range(len(model_list)):
    model=model_list[i].split('/')[1]

    #create the null spectrum for t_ref
    
    null_newf=atmosatm(model_list[i],tel=tel,filebase=model,null_spec=True,removed_gas="H2",star='G')
    psgspec(model,null_newf,showplot=False,null_spec=True)
    null_rad=curr_dir+'/psg_output/%s_null_rad.txt' % model

    #create radiance spectrum

    newf=atmosatm(model_list[i],tel=tel,filebase=model,star='G')
    psgspec(model,newf,showplot=False)
    model_rad=curr_dir+'/psg_output/%s_rad.txt' % model
    t_ref = compute_t_ref(filenames=(model_rad,null_rad), t_exp=10, wl_min=0.4, wl_max=0.9)
    t_ref_list.append(t_ref)
    print("Required exposure time for %s: {:.1f} hr".format(t_ref) % model)

plt.plot(H2_flux_list,t_ref_list)
plt.xlabel('H2 flux (molecules/cm$^2$/s)')
plt.ylabel('Required observation time (hrs)')
plt.scatter(H2_flux_list,t_ref_list)
plt.xscale('log')
plt.yscale('log')
plt.title(r'Observation time vs H$_2$ flux')
#plt.show()
plt.savefig('/home/tessa/Alien_Earths/bioverse-atmos-integration/figures/H2 vs obs time.jpg')


