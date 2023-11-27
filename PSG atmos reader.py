import sys, os
import numpy as np
import matplotlib.pyplot as plt
import glob
import shutil
from tqdm import tqdm, trange


def atmosatm(model, location,filebase = ''):

	# Read the atmos profiles ------------------------------
	ptfile=model+'/profile.pt'
	hcfile=model+'/hcaer.out'
	data = np.genfromtxt(ptfile, skip_header=1)
	if data.shape[1]!=16:
		print("Invalid file dimensions: {}.".format(data.shape))
		exit()
	#Endif
	bin = 4
	mols = 'Altitude,H2O,CH4,C2H6,CO2,O2,O3,CO,H2CO,HNO3,NO2,SO2,N2O,N2'; ncols=16     
	data[:,0] = data[:,0]/1e5   # Convert altitude from cm to km
	dens = data[::bin,2]        # Density [molecules/cm3]
	nalt = int(data.shape[0]/bin)
	gprof = np.zeros([nalt,19])
	gprof[:,0:15] = data[::bin,[3,1,0,4,5,6,7,8,9,10,11,12,13,14,15]]
	gprof[:,15] = 1.0 - np.sum(gprof[:,3:],axis=1) # Compute N2 abundance
	gprof[:,17] = 1e-6
	mmol = [18.0,16.0,30.0,44.0,32.0,48.0,28.0,30.0,63.0,46.0,64.0,44.0,28.0]
	gatm = np.sum(gprof[1,3:16] * mmol)/np.sum(gprof[1,3:16])

	 # Read the haze profile --------------------------------
	if os.path.isfile(hcfile):
		data = np.genfromtxt(hcfile, skip_header=4)
		if data.shape[1]!=9:
			print("Invalid file dimensions: {}.".format(data.shape))
			exit()
		#Endif
		rho = 0.64                                               # Particle density [g/cm3]
		data[:,0] = data[:,0]/1e5                                # Convert altitude from cm to km
		naero = np.interp(gprof[:,2], data[:,0], data[:,1])      # Particle abundance [#/cm3]
		raero = np.interp(gprof[:,2], data[:,0], data[:,2])      # Particle radius [cm]
		vaero = np.pi*(4.0/3.0)*raero**3.0                       # Volume of particle [cm3]
		maero = naero*vaero*rho                                  # Haze mass density [g/cm3]
		xaero = maero/(gatm*dens/6.0221409e+23)                  # Haze abundance [g/g]
		gprof[:,16] = xaero
		gprof[:,17] = raero/1e2
		mols = '%s,Haze,Haze_size' % mols; ncols=19
	#Endif

	  # Write configuration file --------------------------------

	newf = []
	newf.append('<ATMOSPHERE-DESCRIPTION>ATMOS Photochemistry')  # Description establishing the source/reference for the vertical profile
	newf.append('<ATMOSPHERE-STRUCTURE>Equilibrium')             # The structure of the atmosphere, None / Equilibrium:'Hydrostatic equilibrium' / Coma:'Cometary expanding coma'
	newf.append('<ATMOSPHERE-PRESSURE>%.3f' % gprof[0,0])        # For equilibrium atmospheres, this field defines the surface pressure; while for cometary coma, this field indicates the gas production rate
	newf.append('<ATMOSPHERE-PUNIT>bar')                         # The unit of the ATMOSPHERE-PRESSURE field
	newf.append('<ATMOSPHERE-WEIGHT>%.3f' % gatm)                # Molecular weight of the atmosphere [g/mol] or expansion velocity [m/s] for expanding atmospheres
	newf.append('<ATMOSPHERE-NGAS>11')                           # Number of gas in Atmos profile.pt
	newf.append('<ATMOSPHERE-GAS>N2,H2O,CH4,CO2,O2,O3,CO,H2CO,NO2,SO2,N2O') # Name of the gases to include in the simulation,
	newf.append('<ATMOSPHERE-TYPE>HIT[22],HIT[1],HIT[6],HIT[2],HIT[7],HIT[3],HIT[5],HIT[20],HIT[10],HIT[9],HIT[4]') # Sub-type of the gases, e.g. 'HIT[1], HIT[2]'
	newf.append('<ATMOSPHERE-ABUN>1,1,1,1,1,1,1,1,1,1,1')      # Value of the scaler that would multiply the gas profile
	newf.append('<ATMOSPHERE-UNIT>scl,scl,scl,scl,scl,scl,scl,scl,scl,scl,scl') # Scaler unit

	# Add haze parameters

	if ncols>17:
		newf.append('<ATMOSPHERE-NAERO>1')                       # Number of aerosols to include in the simulation
		newf.append('<ATMOSPHERE-AEROS>Haze')                    # Name of the aerosols to include in the simulation
		newf.append('<ATMOSPHERE-ATYPE>KHARE_TITAN_THOLINS_HRI') # Sub-type of the aerosols
		newf.append('<ATMOSPHERE-AABUN>1')                       # Abundance of aerosols. The values can be assumed to be same across all altitudes/layers [%,ppm,ppb,ppt,Kg/m2], or as a multiplier [scaler] to the provided vertical profile
		newf.append('<ATMOSPHERE-AUNIT>scl')                     # Unit of the ATMOSPHERE-AABUN field, % / ppmv / ppbv / pptv / m2:'molecules/m2' / scl:'scaler of profile
		newf.append('<ATMOSPHERE-ASIZE>1')                       # Effective radius of the aerosol particles. The values can be assumed to be same across all layers [um, m, log(um)], or as a multiplier [scaler] to the provided size vertical profile
		newf.append('<ATMOSPHERE-ASUNI>scl')                     # Unit of the size of the aerosol particles
	else:
		newf.append('<ATMOSPHERE-NAERO>0')                       # Number of aerosols to include in the simulation
	#Endif

	# Add profile

	newf.append('<ATMOSPHERE-LAYERS-MOLECULES>%s' % mols)        # Molecules and aerosols quantified by the vertical profile
	newf.append('<ATMOSPHERE-LAYERS>%d' % nalt)                  # Number of layers of the atmospheric vertical profile
	for i in range(nalt):                                        # Values for that specific layer: Pressure[bar], Temperature[K], gases[mol/mol], aerosols [kg/kg] - Optional fields: Altitude[km], X_size[m, aerosol X particle size]
		str = '%.5e' % gprof[i,0]
		for j in range(1,ncols): str = '%s,%.5e' % (str,gprof[i,j])
		newf.append('<ATMOSPHERE-LAYER-%d>%s' % (i+1,str))
	#Endfor

	# Add surface parameters

	newf.append('<SURFACE-TEMPERATURE>%.3f' % gprof[0,1])        # Temperature of the surface [K]
	newf.append('<SURFACE-ALBEDO>0.25')                          # Albedo the surface [0:non-reflectance, 1:fully-reflective]
	newf.append('<SURFACE-EMISSIVITY>1.0')                       # Emissivity of the surface [0:non-emitting, 1:perfect-emitter]
	newf.append('<SURFACE-NSURF>0')                              # Number of components describing the surface properties [areal mixing]

	# Save atmospheric file
	if location=="home":
		dir_root='/media/tessa/Storage/'
	if location=="office":
		dir_root="/home/tessa/"
	if len(filebase):
		model=model.split('/')[1]
		model_file="%s_cfg.txt" % model
		config_name=dir_root+"Alien Earths/Bioverse atmosphere integration/psg_configs/"+model_file
		fw = open(config_name,'w')
		for line in newf: fw.write('%s\n' % line)
		fw.close()
	#Endif
	
	return newf

def psgspec(model, location,config=[], showplot=True):
	#psgurl = 'http://localhost:3000' # URL of the PSG server - For PSG/Docker
	#psgurl = 'http://localhost' # URL of the test PSG server
	psgurl = 'https://psg.gsfc.nasa.gov' # URL of the PSG server
	if location=="home":
		dir_root='/media/tessa/Storage/'
	if location=="office":
		dir_root="/home/tessa/"	# Save configuration and run via the API
	model_file="%s_cfg.txt" % model
	config_name=dir_root+"Alien Earths/Bioverse atmosphere integration/psg_configs/"+model_file
	target=dir_root+"Alien Earths/Bioverse atmosphere integration/"+model_file
	# Shuffle files around temporarily to keep everything tidy
	shutil.copyfile(config_name, target)
	# Actually query PSG and get results
	fw = open(model_file,'w')
	for line in newf: fw.write('%s\n' % line)
	fw.close()
	os.system('curl -s --data-urlencode file@%s_cfg.txt %s/api.php > %s_rad.txt' % (model, psgurl, model))
	dest_path=dir_root+"/Alien Earths/Bioverse atmosphere integration/psg_output/"+"%s_rad.txt" % model
	# More file shuffling
	shutil.copyfile("%s_rad.txt" % model,dest_path)
	os.remove("%s_cfg.txt" % model)
	os.remove("%s_rad.txt" % model)
	# Plot results

	if showplot:
		xunit = 'Wavelength [um]'; yunit = 'Contrast [ppm]'; cols = 'Total'
		fw = open(dest_path); lines = fw.readlines(); fw.close()
		for line in lines:
			if line[0]!='#': break
			if line[0:16]=='# Spectral unit:': xunit = line[17:-1]
			if line[0:16]=='# Radiance unit:': yunit = line[17:-1]
			if line[0:11]=='# Wave/freq': cols = line[12:-1]
		#Endfor
		data = np.genfromtxt(dest_path)
		cols = cols.split(); wnoise=0
location="office"
model_list=glob.glob('sample atmos results/*')
with trange(len(model_list)) as t:
	for i in t:
		model=model_list[i].split('/')[1]
		newf=atmosatm(model_list[i],location,filebase=model)
		psgspec(model,location)
