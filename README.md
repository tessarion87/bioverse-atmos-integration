# Integration of atmospheric modeling into bioverse

![A diagram of the pipeline for integrating atmos into bioverse](https://github.com/tessarion87/bioverse-atmos-integration/blob/master/bioverse%20integration.png)

This is a set of several scripts and packages that will integrate atmospheric modeling (in the form of [atmos](https://github.com/VirtualPlanetaryLaboratory/atmos/tree/master) or [PyAtmos](https://github.com/PyAtmos/PyAtmos) into the bioverse pipeline. 

The primary addition is the package `atmos_to_PSG` developed by Estelle Janin and Tessa Fisher, which will take the profile.pt and hcaer.out files produced by an atmos model and use them to generate a configuration model for the [NASA Planetary Spectrum Generator](https://psg.gsfc.nasa.gov/index.php) and submit it. The resulting spectral data file can then used, for example, as input for bioverse's `t_ref` function.

Note: if you're using PyAtmos, insert the following code at line 401 of `simulation.py`:
```
     self._copy_container_file(self._atmos_directory+'/PHOTOCHEM/OUTPUT/hcaer.out', output_directory)
     self._copy_container_file(self._atmos_directory+'/PHOTOCHEM/OUTPUT/profile.pt', output_directory)
```

## Using atmos\_to\_PSG

`atmos_to_PSG` is made up of two functions, `atmosatm` and `psgspec`. The functions are designed to run in sequence one after another--first `atmosatm`, then `psgspec`.

### atmosatm (model, tel='',filebase='',null_spec=False,removed\_gas='') 
Creates a configuration file from the profile.pt and hcaer.out files. Parameters:

- _model_: the name of the folder containing the profile.pt and hcaer.out files
- _tel_: Optional. The telescope configuration used in generating the spectrum. Options: 'JWST', 'LUVOIR', 'HabEx', and 'Nautilus'
- _filebase_ : Optional. The prefix used on the name of the resulting configuration file
- _null\_spec_ : Determines whether or not the resulting configuration has the atmospheric species of interest removed, for use with bioverse's `t_ref` function. Requires a value for removed_gas (see below)
- _removed\_gas_: The atmospheric species of interest that's removed for use with `t_ref`. Options: 'H2O', 'CH4', 'C2H6', 'CO2', 'O2', 'O3', 'CO', 'H2CO', 'HNO3', 'NO2', 'SO2', 'N2O', and 'N2'
    
### psgspec (model,newf,showplot=True,null_spec=False)
- _model_: the name desired for output file(s)
- _newf_: the configuration file generated by atmosatm
- _showplot_: Determines whether or not psgspec generates and displays a spectral radiance plot
- _null\_spec_: Determines whether or not the configuration file for the null spectrum is used, for use with bioverse's `t_ref` function.

Required packages:
- NumPy

To do:
- Include noise parameters in configuration file
- Add more telescopes as options
