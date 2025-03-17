# Input files

CABLE can be used in many configurations, in part determined by the input files supplied to the code.
Table 1 lists the various input files used in offline CABLE. A description of each file can be accessed through the left navigation bar.

## Table 1: CABLE input files for the offline case

|   Input file         	 | Description |
|------------------------|-------------|
| cable.nml            	 | main configuration file for CABLE |
| pft_params.nml       	 | default parameter values for each PFT |
| cable_soilparm.nml   	 | default parameter values for each soil type |
| pftlookup.csv        	 | default parameter values for CASA-CNP |          
| Meteorological forcing | atmospheric forcing data for CABLE |
| Surface forcing        | information about the surface characteristics |
| Restart                | information from a previous CABLE run to restart a simulation | 

## Example configurations

CABLE is used as a standalone model or in the ACCESS Earth System model. You can find example configurations for each of its applications at:

- CABLE standalone: files under CABLE/src/offline
- ACCESS-ESM1.5: ACCESS-NRI supported configurations at [access-esm1.5-configs repository][access-esm1.5-configs]

[access-esm1.5-configs]: https://github.com/ACCESS-NRI/ACCESS-esm1.5-configs/