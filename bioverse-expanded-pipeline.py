import glob
from bioverse.util import compute_t_ref
from atmos_to_PSG import atmosatm, psgspec
import pyatmos
import os
from tqdm import tqdm, trange

atmos = pyatmos.Simulation(code_path='/media/tessa/Storage/Dissertation_backup/atmos-master', DEBUG=True)
atmos.start()

curr_dir=os.getcwd()

biotic_flux_list=[9.64E+11]


for i in range(len(biotic_flux_list)):
	CH4_flux=biotic_flux_list[i]
	print(CH4_flux)
	output_dir=curr_dir+'/sample_atmos_results/integration_test_{n}'.format(n=i)
	args = {
		'species_fluxes': {'CH4' :CH4_flux},
		'max_photochem_iterations' : 50000, 
		'max_clima_steps' : 10, 
		'output_directory' : output_dir}
	
	
	atmos.run(**args)

model_list=glob.glob('sample_atmos_results/integration_test*')

tel='LUVOIR'
# generate spectra and then calculate t_ref
for i in range(len(model_list)):
	model=model_list[i].split('/')[1]
	
	#create the null spectrum for t_ref
	
	null_newf=atmosatm(model_list[i],tel=tel,filebase=model,null_spec=True,removed_gas="CH4")
	psgspec(model,null_newf,showplot=True,null_spec=True)
	null_rad=curr_dir+'/psg_output/%s_null_rad.txt' % model

	#create radiance spectrum

	newf=atmosatm(model_list[i],tel=tel,filebase=model)
	psgspec(model,newf,showplot=True)
	model_rad=curr_dir+'/psg_output/%s_rad.txt' % model
	t_ref = compute_t_ref(filenames=(model_rad,null_rad), t_exp=100, wl_min=0.4, wl_max=0.9,)
	print("Required exposure time for %s: {:.1f} hr".format(t_ref) % model)
	
