import sys, os
import numpy as np
import matplotlib.pyplot as plt
import shutil

def atmosatm(model,tel='',filebase = '',null_spec=False,removed_gas='',star=''):

	dir_root=os.getcwd()

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
	abun_list=[]
	object_name='<OBJECT-NAME>'+model.rsplit('/', 1)[-1]
	if null_spec==True:
		object_name=object_name+'_null'
	newf = []
	newf.append(object_name)
	if star !='':
		star_object='<OBJECT-STAR-TYPE>'+star
	newf.append('<ATMOSPHERE-DESCRIPTION>ATMOS Photochemistry')  # Description establishing the source/reference for the vertical profile
	newf.append('<ATMOSPHERE-STRUCTURE>Equilibrium')             # The structure of the atmosphere, None / Equilibrium:'Hydrostatic equilibrium' / Coma:'Cometary expanding coma'
	newf.append('<ATMOSPHERE-PRESSURE>%.3f' % gprof[0,0])        # For equilibrium atmospheres, this field defines the surface pressure; while for cometary coma, this field indicates the gas production rate
	newf.append('<ATMOSPHERE-PUNIT>bar')                         # The unit of the ATMOSPHERE-PRESSURE field
	newf.append('<ATMOSPHERE-WEIGHT>%.3f' % gatm)                # Molecular weight of the atmosphere [g/mol] or expansion velocity [m/s] for expanding atmospheres
	newf.append('<ATMOSPHERE-NGAS>11')                           # Number of gas in Atmos profile.pt
	gases='N2,H2O,CH4,CO2,O2,O3,CO,H2CO,NO2,SO2,N2O'
	gas_list=gases.split(',')
	for i in range(len(gas_list)):
		if null_spec==True and gas_list[i]==removed_gas:
			abun_list.append('0')
		else:
			abun_list.append('1')
	abuns=','.join(abun_list)
	newf.append('<ATMOSPHERE-GAS>'+gases) #Name of the gases to include in the simulation
	newf.append('<ATMOSPHERE-TYPE>HIT[22],HIT[1],HIT[6],HIT[2],HIT[7],HIT[3],HIT[5],HIT[20],HIT[10],HIT[9],HIT[4]') # Sub-type of the gases, e.g. 'HIT[1], HIT[2]'
	newf.append('<ATMOSPHERE-ABUN>'+abuns)      # Value of the scaler that would multiply the gas profile
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
	
	# Add telescope parameters
	if tel!='':
		if tel=='HabEx':
			tel_file='HabEx.dat'
		elif tel=='LUVOIR':
			tel_file='LUVOIR.dat'
		elif tel=='Nautilus':
			tel_file='Nautilus.dat'
		elif tel=='JWST':
			tel_file='JWST.dat'
		elif tel=='iSHELL':
			tel_file='iSHELL.dat'
		elif tel=='Keck_HIRES':
			tel_file='Keck_HIRES.dat'
		elif tel=='HST':
			tel_file='HST.dat'
		tel_path=dir_root+'/telescope_config/'+tel_file
		tel_cfg=open(tel_path,'r')
		lines=tel_cfg.readlines()
		for line in lines:
			newf.append(line)
	
	# Save atmospheric file

	if len(filebase):
		model=model.split('/')[1]
		if null_spec==True:
			model_file="%s_null_cfg.txt" % model
		else:
			model_file="%s_cfg.txt" % model
		config_name=dir_root+"/psg_configs/"+model_file
		fw = open(config_name,'w')
		for line in newf: fw.write('%s\n' % line)
		fw.close()
	#Endif
	
	return newf

def psgspec(model,newf, showplot=True,null_spec=False):
	psgurl = 'http://localhost:3000' # URL of the PSG server - For PSG/Docker
	#psgurl = 'http://localhost' # URL of the test PSG server
	#psgurl = 'https://psg.gsfc.nasa.gov' # URL of the PSG server
	dir_root=os.getcwd()
	# Save configuration and run via the API
	if null_spec==True:
		config_file='%s_null_cfg.txt' % model
	else:
		config_file="%s_cfg.txt" % model
	config_name=dir_root+"/psg_configs/"+config_file
	target=dir_root+'/'+config_file
	# Shuffle files around temporarily to keep everything tidy
	shutil.copyfile(config_name, target)
	# Actually query PSG and get results
	fw = open(config_name,'w')
	for line in newf: fw.write('%s\n' % line)
	fw.close()
	if null_spec==True:
		rad_file='%s_null_rad.txt' % model
	else:
		rad_file='%s_rad.txt' % model
	os.system('curl -s --data-urlencode file@%s %s/api.php > %s' % (config_file, psgurl, rad_file))
	dest_path=dir_root+"/psg_output/"+rad_file
	# More file shuffling
	shutil.copyfile(rad_file,dest_path)

	# Plot results
	if showplot:
		xunit = 'Wavelength [um]'; yunit = 'Contrast [ppm]'; cols = 'Total'
		fw = open(rad_file); lines = fw.readlines(); fw.close()
		for line in lines:
			if line[0]!='#': break
			if line[0:16]=='# Spectral unit:': xunit = line[17:-1]
			if line[0:16]=='# Radiance unit:': yunit = line[17:-1]
			if line[0:11]=='# Wave/freq': cols = line[12:-1]
		#Endfor
		
		data = np.genfromtxt(rad_file)
		cols = cols.split(); wnoise=0
		plt.figure(figsize=[10,5])
		for i in range(len(cols)):
			if cols[i]=='Noise': wnoise=i; continue
			plt.plot(data[:,0],data[:,i+1],label=cols[i])
		#Endfors
		plt.xscale('log')
		plt.xlabel(xunit)
		plt.ylabel(yunit)
		plt.ylim([max(data[:,1])/100.0,max(data[:,1])*1.2])
		plt.xlim([min(data[:,0]),max(data[:,0])])
		plt.legend()
		plt.title('ATMOS transmission spectrum generated by NASA/PSG')
		plt.tight_layout()
		if null_spec==True:
			png_file='%s_null_rad.png' % model
		else:
			png_file='%s_rad.png' % model
		figname=dir_root+"/figures/"+png_file
		plt.savefig(figname)
		plt.show()
	os.remove(config_file) # file clean-up
	os.remove(rad_file)
	

