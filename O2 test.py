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
import pandas as pd

# atmos = pyatmos.Simulation(code_path='~/atmos-master',docker_image=None, DEBUG=True)
# atmos.start()

curr_dir=os.getcwd()

O2_flux_list=[5.0E+8,5.0E+8,8.0E+8,1.0E+9,2.5E+9,3.5E+9,5.0E+9,6.0E+9,7.5E+9,8.5E+9,1.00E+10,5.0E+10,9.0E+10,1.50E+11,2.50E+11,
              3.0E+11,3.25E+11,3.5E+11,3.95E+11,4.2E+11,5.50E+11,1.0E+12,5.0E+12,9.0E+12,2.4E+13,5.0E+13,6.0E+13]

# for i in range(len(O2_flux_list)):
    # O2_flux=O2_flux_list[i]
    # output_dir=curr_dir+'/sample_atmos_results/O2_test_{n}'.format(n=i)
    # args = {
        # 'species_fluxes': {'O2' :O2_flux},
        # 'max_photochem_iterations' : 50000, 
        # 'max_clima_steps' : 10, 
        # 'output_directory' : output_dir}
    
    
    # atmos.run(**args)

model_list=[]
for i in range(len(O2_flux_list)):
    model_name="sample_atmos_results/O2_test_"+str(i)
    model_list.append(model_name)


# generate sample of observed planets

generator=Generator('transit')
survey=TransitSurvey('default')
sample, detected, data = survey.quickrun(generator, t_total=10*365.25)

# generate spectra and then calculate t_ref
t_ref_list_JWST=[]
tel='JWST'

for i in range(len(model_list)):
    model=model_list[i].split('/')[1]

    #create the null spectrum for t_ref
    
    null_newf=atmosatm(model_list[i],tel=tel,filebase=model,null_spec=True,removed_gas="O2",star='G')
    psgspec(model,null_newf,showplot=False,null_spec=True)
    null_rad=curr_dir+'/psg_output/%s_null_rad.txt' % model

    #create radiance spectrum

    newf=atmosatm(model_list[i],tel=tel,filebase=model,star='G')
    psgspec(model,newf,showplot=False)
    model_rad=curr_dir+'/psg_output/%s_rad.txt' % model
    t_ref = compute_t_ref(filenames=(model_rad,null_rad), t_exp=10, wl_min=0.4, wl_max=0.9)
    t_ref_list_JWST.append(t_ref)
    print("Required exposure time for %s: {:.1f} hr".format(t_ref) % model)

O2_flux_list.pop(0)
t_ref_list_JWST.pop(0) #since the first entry of this list always seem weirdly low


t_ref_list_nautilus=[]
tel='Nautilus'

for i in range(len(model_list)):
    model=model_list[i].split('/')[1]

    #create the null spectrum for t_ref
    
    null_newf=atmosatm(model_list[i],tel=tel,filebase=model,null_spec=True,removed_gas="O2",star='G')
    psgspec(model,null_newf,showplot=False,null_spec=True)
    null_rad=curr_dir+'/psg_output/%s_null_rad.txt' % model

    #create radiance spectrum

    newf=atmosatm(model_list[i],tel=tel,filebase=model,star='G')
    psgspec(model,newf,showplot=False)
    model_rad=curr_dir+'/psg_output/%s_rad.txt' % model
    t_ref = compute_t_ref(filenames=(model_rad,null_rad), t_exp=10, wl_min=0.4, wl_max=0.9)
    t_ref_list_nautilus.append(t_ref)
    print("Required exposure time for %s: {:.1f} hr".format(t_ref) % model)

t_ref_list_nautilus.pop(0) #since the first entry of this list always seem weirdly low


plt.plot(O2_flux_list,t_ref_list_JWST, label='JWST')
plt.scatter(O2_flux_list,t_ref_list_JWST)
plt.plot(O2_flux_list,t_ref_list_nautilus,label='Nautilus')
plt.scatter(O2_flux_list,t_ref_list_nautilus)
plt.xscale('log')
plt.yscale('log')
plt.vlines(2.39e13,0,1e10,linestyles='dashed',label='Modern Earth')
plt.ylabel('Required observation time (hrs)')
plt.xlabel(r'O$_2$ flux (molecules/cm$^2$/s)')
plt.suptitle(r'Observation time vs O$_2$ flux')
plt.xlim(1e10,1e14)
plt.legend()
#plt.show()
plt.savefig('/home/tessa/Alien Earths/bioverse-atmos-integration/figures/O2 vs obs time combined.jpg')

O2_table=pd.DataFrame()
O2_table.insert(0,'O2 surface flux',O2_flux_list)
O2_table.insert(1,'t_ref (JWST)',t_ref_list_JWST)
O2_table.insert(2,'t_ref (Nautilus)',t_ref_list_nautilus)
O2_table.to_csv('O2 t_ref lookup table.csv')
