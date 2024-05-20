# pft_params.nml options

In CABLE the characteristics of the vegetation cover are parametrised and generally remain constant throughout the simulation. These **vegetation parameters** in general vary by **vegetation type**.
CABLE supports tiling of vegetation, i.e. multiple vegetation types with specified area fractions for a single location/grid-cell (also referred to as patches in CABLE). 
Altogether CABLE considers the 17 distinct types of vegetation cover listed in Table 1.
Eleven of those are plant functional types (PFT), two are generally unoccupied, and the remaining four are not vegetated. 

## Table 1: CABLE vegetation types

| Number |        Type          |
|--------|----------------------|
|    1 	 | Evergreen Needleleaf |
|    2 	 | Evergreen Broadleaf  |
|    3 	 | Deciduous Needleleaf |
|    4 	 | Deciduous Broadleaf  |
|    5 	 | Shrub 		|
|    6 	 | C3 Grassland 	|
|    7 	 | C4 Grassland 	|
|    8 	 | Tundra 		|
|    9 	 | C3 Cropland 	        |
|    10	 |  C4 Cropland 	|
|    11	 |  Wetland 		|
|    12	 |  empty 		|
|    13	 |  empty 		|
|    14	 |  Barren 		|
|    15	 |  Urban 		|
|    16	 |  Lakes 		|
|    17	 |  Ice 		|

The CABLE vegetation parameters are listed in Table 2 along with their value ranges, units, and descriptions. All parameters in Table 2 are of type real and dimension 17, corresponding to the number of vegetation types in Table 1.
The ranges listed are “physically possible” and mostly correspond to the optional checks in the CABLE offline code.
## Table 2: CABLE vegetation parameters

| Name (vegin%)         | Range            | Units                                                                                               | Description |
|-----------------------|------------------|-----------------------------------------------------------------------------------------------------|-------------|                                                                      
| 	a1gs            | >=0.0        	   | \( (-) \) 												 | a1 parameter in stomatal conductance model. Represents the sensitivity of stomatal conductance to the assimilation rate |
| 	alpha           | >=0.0        	   | \( (mol (electrons) \cdot mol^{-1} (photons) (C3) \cdot mol (CO2) \cdot mol^{-1} (photons) (C4)) \) | Initial slope of J-Q response curve |
| 	canst1          | 0.05 – 0.15  	   | \( (mm \cdot LAI^{-1}) \) 										 | Maximum intercepted water by canopy |                                   
| 	cfrd            | >=0.0        	   | \( (-) \) 												 | Ratio of day respiration to vcmax |                                                      
| 	clitt           | >=0.0        	   | \( (tC \cdot ha^{-1}) \)                                                                            | Leaf litter (alters resistance to soil evaporation) |
| 	conkc0 	        | >=0.0        	   | \( (bar) \) 											 | Michaelis-Menton constant for carboxylase |
| 	conko0 	        | >=0.0        	   | \( (bar) \) 											 | Michaelis-Menton constant for oxygenase   |                                           
| 	convex 	        | >=0.0        	   | \( (-) \) 												 | Convexity of J-Q response curve |                                                       
| 	cplant1         | >=0.0        	   | \( (gC \cdot m^{-2}) \) 										 | Plant carbon in 1st vegetation carbon store |
| 	cplant2         | >=0.0        	   | \( (gC \cdot m^{-2}) \) 										 | Plant carbon in 2nd vegetation carbon store |
| 	cplant3         | >=0.0        	   | \( (gC \cdot m^{-2}) \) 										 | Plant carbon in 3rd vegetation carbon store |                             
| 	csoil1 	        | >=0.0        	   | \( (gC \cdot m^{-2}) \) 										 | Soil carbon in 1st soil carbon store |
| 	csoil2 	        | >=0.0        	   | \( (gC \cdot m^{-2}) \) 										 | Soil carbon in 2nd soil carbon store |                                    
| 	d0gs            | >=0.0        	   | \( (kPa) \) 											 | d0 in stomatal conductance model |                                                    
| 	ejmax           | 0 – 3.0E-4   	   | \( (mol \cdot m^{-2} \cdot s^{-1}) \) 								 | Maximum potential electron transp rate top leaf |           
| 	ekc 	        | >=0.0        	   | \( (J \cdot mol^{-1}) \) 										 | Activation energy for carboxylase |
| 	eko 	        | >=0.0        	   | \( (J \cdot mol^{-1}) \) 										 | Activation energy for oxygenase   |                                      
| 	extkn 	        | 0.0 – 10.0   	   | \( (-) \) 												 | Extinction coeficient for vertical profile of N |                                       
| 	frac4 	        | 0.0 – 1.0    	   | \( (-) \) 												 | Fraction of c4 plants |                                                                 
| 	froot1 	        | 0.0 – 1.0    	   | \( (-) \) 												 | Fraction of root in soil layer 1 |
| 	froot2 	        | 0.0 – 1.0    	   | \( (-) \) 												 | Fraction of root in soil layer 2 |
| 	froot3 	        | 0.0 – 1.0    	   | \( (-) \) 												 | Fraction of root in soil layer 3 |
| 	froot4 	        | 0.0 – 1.0    	   | \( (-) \) 												 | Fraction of root in soil layer 4 |
| 	froot5 	        | 0.0 – 1.0    	   | \( (-) \) 												 | Fraction of root in soil layer 5 |
| 	froot6 	        | 0.0 – 1.0    	   | \( (-) \) 												 | Fraction of root in soil layer 6 |                                                        
| 	g0     	        | -0.5 – 0.5   	   | \( (mol \cdot m^{-2} \cdot s^{-1}) \) 								 | Belinda's stomatal model intercept. Residual stomatal conductance as net assimilation rate reaches zero |
| 	g1     	        | 0.0 – 20.0   	   | \( (kPa) \) 											 | Belinda's stomatal model slope. Sensitivity of stomatal conductance to the assimilation rate |
| 	gswmin 	        | >=0.0        	   | \( (mol \cdot m^{-2} \cdot s^{-1}) \) 								 | Minimal stomatal conductance |                              
| 	hc     	        | 0.0 – 100.0  	   | \( (m) \) 												 | Roughness height of canopy (veg - snow) |                                               
| 	lai    	        | 0.0 – 8.0    	   | \( (-) or (m^{2} \cdot m^{-2} ) \) 								 | Leaf area index of each plant functional type |                
| 	length 	        | >=0.0        	   | \( (m) \) 												 | Leaf length |                                                                           
| 	ratecp1         | 0.01 – 3.0   	   | \( (year^{-1}) \) 											 | Plant carbon pool rate constant in 1st vegetation carbon store  |
| 	ratecp2         | 0.01 – 3.0   	   | \( (year^{-1}) \) 											 | Plant carbon pool rate constant in 2nd vegetation carbon store  |
| 	ratecp3         | 0.01 – 3.0   	   | \( (year^{-1}) \) 											 | Plant carbon pool rate constant in 3rd vegetation carbon store  |
| 	ratecs1         | 0.01 – 3.0   	   | \( (year^{-1}) \) 											 | Soil carbon pool rate constant in 1st soil carbon store 	 |
| 	ratecs2         | 0.01 – 3.0   	   | \( (year^{-1}) \) 											 | Soil carbon pool rate constant in 2nd soil carbon store 	 |                
| 	refl1 	        | 0.0 – 0.5    	   | \( (-) \) 												 | Leaf reflectance in 1st radiation band |
| 	refl2 	        | 0.0 – 0.5    	   | \( (-) \) 												 | Leaf reflectance in 2nd radiation band |
| 	refl3 	        | 0.0 – 0.5    	   | \( (-) \) 												 | Leaf reflectance in 3rd radiation band |                                                
| 	rootbeta        | 0.0 – 1.0    	   | \( (-) \) 												 | Beta parameter to calculate froot (Jackson et al. 1996) |                               
| 	rp20            | 0.0 – 10.0   	   | \( (-) \) 												 | Plant respiration coefficient at 20 C |                                                  
| 	rpcoef          | 0.05 – 1.5   	   | \( (^{\circ}C^{-1}) \) 										 | Temperature coefficient non-leaf plant respiration |                       
| 	rs20   	        | 0.0 – 10.0   	   | \( (mol \cdot m^{-2} \cdot s^{-1} ) \) 								 | Soil respiration at 20 C |                                 
| 	shelrb 	        | 1.0 – 3.0    	   | \( (-) \) 												 | Sheltering factor |                                                                     
| 	taul1 	        | 0.0 – 0.3    	   | \( (-) \) 												 | Leaf transmittance in 1st radiation band |
| 	taul2 	        | 0.0 – 0.3    	   | \( (-) \) 												 | Leaf transmittance in 2nd radiation band |
| 	taul3 	        | 0.0 – 0.3    	   | \( (-) \) 												 | Leaf transmittance in 3rd radiation band |                                              
| 	tmaxvj 	        | -15.0 – 30.0 	   | \( (^{\circ}C) \) 											 | Maximum temperature of the start of photosynthesis |
| 	tminvj 	        | -20.0 – 15.0 	   | \( (^{\circ}C) \) 											 | Minimum temperature of the start of photosynthesis |                            
| 	vbeta 	        | -999999 – 999999 | \( (-) \) 												 | Stomatal sensitivity to soil water |                                                    
| 	vcmax 	        | 5.0E-6 – 1.5E-4  | \( ( mol \cdot m^{-2} \cdot s^{-1} ) \) 								 | Maximum RuBP carboxylation rate top leaf |                
| 	vegcf 	        | 0.0 – 100.0      | \( (-) \) 												 | Scaling on soil respiration (place-holder scheme) |                                     
| 	width 	        | >=0.0 	   | \( (m) \) 												 | Leaf width |
| 	xfang 	        | -1.0 – 0.5 	   | \( (-) \) 												 | Leaf angle |                                                                            
| 	zr              | >=0.0 	   | \( (cm) \) 											 | Maximum rooting depth |                                                                



The CABLE default values of the vegetation parameters in Table 2 for the vegetation types in Table 1 are listed in Table 3.

## Table 3: Default vegetation parameter values

|                | Evergreen Needleleaf | Evergreen Broadleaf | Deciduous Needleleaf | Deciduous Broadleaf | Shrub | C3 Grassland | C4 Grassland | Tundra | C3 Cropland | C4 Cropland | Wetland | empty | empty | Barren | Urban | Lakes | Ice |
|----------------|----------------------|---------------------|----------------------|---------------------|-------|--------------|--------------|--------|-------------|-------------|---------|-------|-------|--------|-------|-------|-----|
| a1gs   	 | 9  |   9  |   9  |   9  |   9  |   9  |   4  |   9  |   9  |   4  |   9  |   9  |   9  |   9  |   9  |   9  |   9 |
| alpha  	 | 0.2  |   0.2  |   0.2  |   0.2  |   0.2  |   0.2  |   0.05  |   0.2  |   0.2  |   0.05  |   0.2  |   0.2  |   0.2  |   0.2  |   0.2  |   0.2  |   0.2 |
| canst1 	 | 0.1  |   0.1  |   0.1  |   0.1  |   0.1  |   0.1  |   0.1  |   0.1  |   0.1  |   0.1  |   0.1  |   0.1  |   0.1  |   0.1  |   0.1  |   0.1  |   0.1 |
| cfrd   	 | 0.015  |   0.015  |   0.015  |   0.015  |   0.015  |   0.025  |   0.015  |   0.015  |   0.025  |   0.015  |   0.015  |   0.015  |   0.015  |   0.015  |   0.015  |   0.015  |   0.015 |
| clitt  	 | 20  |   6  |   10  |   13  |   2  |   2  |   0.3  |   0.3  |   0  |   0  |   2  |   2  |   0  |   0  |   0  |   0  |   0 |
| conkc0 	 | 0.000302  |   0.000302  |   0.000302  |   0.000302  |   0.000302  |   0.000302  |   0.000302  |   0.000302  |   0.000302  |   0.000302  |   0.000302  |   0.000302  |   0.000302  |   0.000302  |   0.000302  |   0.000302  |   0.000302 |
| conko0 	 | 0.256  |   0.256  |   0.256  |   0.256  |   0.256  |   0.256  |   0.256  |   0.256  |   0.256  |   0.256  |   0.256  |   0.256  |   0.256  |   0.256  |   0.256  |   0.256  |   0.256 |
| convex 	 | 0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.8  |   0.01  |   0.01  |   0.8  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01 |
| cplant1  	 | 200  |   300  |   200  |   300  |   159  |   250  |   250  |   250  |   150  |   150  |   250  |   1  |   0.1  |   0  |   1  |   1  |   0 |
| cplant2  	 | 10217  |   16833  |   5967  |   12000  |   5000  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0 |
| cplant3  	 | 876  |   1443  |   511  |   1029  |   500  |   500  |   500  |   500  |   607  |   607  |   500  |   1  |   0.1  |   0  |   1  |   1  |   0 |
| csoil1   	 | 184  |   303  |   107  |   216  |   100  |   275  |   275  |   275  |   149  |   149  |   275  |   1  |   0.1  |   1  |   1  |   1  |   1 |
| csoil2   	 | 367  |   606  |   214  |   432  |   250  |   314  |   314  |   314  |   300  |   300  |   314  |   1  |   0.1  |   1  |   1  |   1  |   1 |
| d0gs     	 | 1500  |   1500  |   1500  |   1500  |   1500  |   1500  |   1500  |   1500  |   1500  |   1500  |   1500  |   1500  |   1500  |   1500  |   1500  |   1500  |   1500 |
| ejmax    	 | 0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0 |
| ekc      	 | 59430  |   59430  |   59430  |   59430  |   59430  |   59430  |   59430  |   59430  |   59430  |   59430  |   59430  |   59430  |   59430  |   59430  |   59430  |   59430  |   59430 |
| eko      	 | 36000  |   36000  |   36000  |   36000  |   36000  |   36000  |   36000  |   36000  |   36000  |   36000  |   36000  |   36000  |   36000  |   36000  |   36000  |   36000  |   36000 |
| extkn    	 | 0.001  |   0.001  |   0.001  |   0.001  |   0.001  |   0.001  |   0.001  |   0.001  |   0.001  |   0.001  |   0.001  |   0.001  |   0.001  |   0.001  |   0.001  |   0.001  |   0.001 |
| frac4    	 | 0  |   0  |   0  |   0  |   0  |   0  |   1  |   0  |   0  |   1  |   0  |   0  |   0  |   0  |   0  |   0  |   0|
| froot1 	 | 0.05  |   0.2  |   0.2  |   0.2  |   0.2  |   0.15  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0 |
| froot2 	 | 0.05  |   0.2  |   0.2  |   0.2  |   0.2  |   0.15  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0 |
| froot3 	 | 0.05  |   0.2  |   0.2  |   0.2  |   0.2  |   0.15  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0 |
| froot4 	 | 0.05  |   0.2  |   0.2  |   0.2  |   0.2  |   0.15  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0 |
| froot5 	 | 0.05  |   0.2  |   0.2  |   0.2  |   0.2  |   0.15  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0 |
| froot6 	 | 0.05  |   0.2  |   0.2  |   0.2  |   0.2  |   0.15  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0 |
| g0 	 	 | 0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0 |
| g1 	 	 | 2.346064  |   4.114762  |   2.346064  |   4.447321  |   4.694803  |   5.2485  |   1.616178  |   2.222156  |   5.789377  |   1.616178  |   5.2485  |   5.2485  |   0  |   5.2485  |   5.2485  |   5.2485  |   5.2485 |
| gswmin   	 | 0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.04  |   0.01  |   0.01  |   0.04  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01 |
| hc       	 | 17  |   35  |   15.5  |   20  |   0.6  |   0.567  |   0.567  |   0.567  |   0.55  |   0.55  |   0.567  |   0.2  |   6.017  |   0.2  |   0.2  |   0.2  |   0.2 |
| lai      	 | 4  |   5  |   0  |   0  |   0  |   0.2  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0  |   0 |
| length   	 | 0.055  |   0.1  |   0.04  |   0.15  |   0.1  |   0.3  |   0.3  |   0.3  |   0.3  |   0.3  |   0.3  |   0.03  |   0.242  |   0.03  |   0.03  |   0.03  |   0.03 |
| ratecp1  	 | 1  |   1  |   1  |   1  |   1  |   1  |   1  |   1  |   1  |   1  |   1  |   1  |   1  |   1  |   1  |   1  |   1 |
| ratecp2  	 | 0.03  |   0.03  |   0.03  |   0.03  |   0.03  |   0.03  |   0.03  |   0.03  |   0.03  |   0.03  |   0.03  |   0.03  |   0.03  |   0.03  |   0.03  |   0.03  |   0.03 |
| ratecp3  	 | 0.14  |   0.14  |   0.14  |   0.14  |   0.14  |   0.14  |   0.14  |   0.14  |   0.14  |   0.14  |   0.14  |   0.14  |   0.14  |   0.14  |   0.14  |   0.14  |   0.14 |
| ratecs1  	 | 2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2 |
| ratecs2  	 | 0.5  |   0.5  |   0.5  |   0.5  |   0.5  |   0.5  |   0.5  |   0.5  |   0.5  |   0.5  |   0.5  |   0.5  |   0.5  |   0.5  |   0.5  |   0.5  |   0.5 |
| refl1 	 | 0.09  |   0.09  |   0.075  |   0.09  |   0.09  |   0.11  |   0.11  |   0.075  |   0.11  |   0.11  |   0.108  |   0.055  |   0.091  |   0.238  |   0.143  |   0.143  |   0.159 |
| refl2 	 | 0.3  |   0.29  |   0.3  |   0.29  |   0.3  |   0.34  |   0.34  |   0.32  |   0.34  |   0.34  |   0.343  |   0.19  |   0.31  |   0.457  |   0.275  |   0.275  |   0.305 |
| refl3 	 | 0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01 |
| rootbeta 	 | 0.943  |   0.962  |   0.966  |   0.961  |   0.964  |   0.943  |   0.943  |   0.943  |   0.961  |   0.961  |   0.943  |   0.975  |   0.961  |   0.961  |   0.961  |   0.961  |   0.961 |
| rp20     	 | 3  |   0.6  |   3  |   2.2  |   1  |   1.5  |   2.8  |   2.5  |   1.5  |   1  |   1.5  |   1  |   1  |   1  |   1  |   1  |   1 |
| rpcoef   	 | 0.0832  |   0.0832  |   0.0832  |   0.0832  |   0.0832  |   0.0832  |   0.0832  |   0.0832  |   0.0832  |   0.0832  |   0.0832  |   0.0832  |   0.0832  |   0.0832  |   0.0832  |   0.0832  |   0.0832 |
| rs20     	 | 1  |   1  |   1  |   1  |   1  |   1  |   1  |   1  |   1  |   1  |   1  |   0  |   1  |   0  |   0  |   0  |   0 |
| shelrb   	 | 2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2  |   2 |
| taul1 	 | 0.09  |   0.09  |   0.075  |   0.09  |   0.09  |   0.11  |   0.11  |   0.075  |   0.11  |   0.11  |   0.075  |   0.023  |   0.059  |   0.039  |   0.023  |   0.023  |   0.026 |
| taul2 	 | 0.3  |   0.29  |   0.3  |   0.29  |   0.3  |   0.34  |   0.34  |   0.32  |   0.34  |   0.34  |   0.146  |   0.198  |   0.163  |   0.189  |   0.113  |   0.113  |   0.113 |
| taul3 	 | 0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01 |
| tmaxvj 	 | -10  |   -10  |   10  |   15  |   -10  |   -10  |   -10  |   -10  |   -10  |   -10  |   -10  |   -10  |   -10  |   -10  |   -10  |   -10  |   -10 |
| tminvj 	 | -15  |   -15  |   5  |   5  |   -15  |   -15  |   -15  |   -15  |   -15  |   -15  |   -15  |   -15  |   -15  |   -15  |   -15  |   -15  |   -15 |
| vbeta 	 | 2  |   2  |   2  |   2  |   4  |   4  |   4  |   4  |   2  |   2  |   4  |   4  |   2  |   4  |   4  |   4  |   4 |
| vcmax 	 | 4E-05  |   5.5E-05  |   4E-05  |   6E-05  |   4E-05  |   6E-05  |   1E-05  |   4E-05  |   8E-05  |   8E-05  |   6E-05  |   1.7E-05  |   1E-06  |   1.7E-05  |   1.7E-05  |   1.7E-05  |   1.7E-05 |
| vegcf 	 | 9  |   14  |   9  |   8  |   5  |   7  |   7  |   5  |   7  |   1  |   7  |   1  |   1  |   1  |   1  |   1  |   1 |
| width 	 | 0.001  |   0.05  |   0.001  |   0.08  |   0.005  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.01  |   0.003  |   0.015  |   0.001  |   0.001  |   0.001  |   0.001 |
| xfang 	 | 0.01  |   0.1  |   0.01  |   0.25  |   0.01  |   -0.3  |   -0.3  |   -0.3  |   -0.3  |   -0.3  |   -0.3  |   0.1  |   0  |   0  |   0  |   0  |   0|
| zr             | 1.8  |   3  |   2  |   2  |   2.5  |   0.5  |   0.5  |   0.5  |   0.5  |   0.5  |   1.8  |   3.1  |   3  |   1  |   1  |   1  |   1 |



The CABLE distribution provides the default vegetation parameter values from Table 3 in the namelist file pft_params.nml, including the vegetation types from Table 1 in the top part.
The chosen parameter values in offline cases can be checked against the pre-defined realistic parameter value ranges using the CABLE namelist variable check%ranges=.true. in cable.nml. 
## Example pft_params.nml file

!!! Note "Namelist file format explanation"

    The asterisk symbol, in the namelist file, means the following entry is repeated multiple times. It does not mean arithmetic multiplication.

``` fortran-free-form
&cable_pftparm
vegin%desc='Evergreen Needleleaf', 'Evergreen Broadleaf', 'Deciduous Needleleaf', 'Deciduous Broadleaf',
                    'Shrub', 'C3 Grassland', 'C4 Grassland', 'Tundra', 'C3 Cropland', 'C4 Cropland',
                    'Wetland', 'empty', 'empty', 'Barren', 'Urban', 'Lakes', 'Ice',
vegin%a1gs=6*9.000000,4.000000,9.000000,9.000000,4.000000,7*9.000000,
vegin%alpha=6*0.200000,0.050000,0.200000,0.200000,0.050000,7*0.200000,
vegin%canst1=17*0.100000,
vegin%cfrd=6*0.015000,0.025000,0.015000,0.015000,0.025000,7*0.015000,
vegin%clitt=20.000000,6.000000,10.000000,13.000000,2.000000,2.000000,0.300000,0.300000,0.000000,0.000000,2.000000,2.000000,5*0.000000,
vegin%conkc0=17*0.000302,
vegin%conko0=17*0.256000,
vegin%convex=6*0.010000,0.800000,0.010000,0.010000,0.800000,7*0.010000,
vegin%cplant1=200.000000,300.000000,200.000000,300.000000,159.000000,250.000000,250.000000,250.000000,150.000000,150.000000,250.000000,1.000000,0.100000,0.000000,1.000000,1.000000,0.000000,
vegin%cplant2=10217.000000,16833.000000,5967.000000,12000.000000,5000.000000,12*0.000000,
vegin%cplant3=876.000000,1443.000000,511.000000,1029.000000,500.000000,500.000000,500.000000,500.000000,607.000000,607.000000,500.000000,1.000000,0.100000,0.000000,1.000000,1.000000,0.000000,
vegin%csoil1=184.000000,303.000000,107.000000,216.000000,100.000000,275.000000,275.000000,275.000000,149.000000,149.000000,275.000000,1.000000,0.100000,1.000000,1.000000,1.000000,1.000000,
vegin%csoil2=367.000000,606.000000,214.000000,432.000000,250.000000,314.000000,314.000000,314.000000,300.000000,300.000000,314.000000,1.000000,0.100000,1.000000,1.000000,1.000000,1.000000,
vegin%d0gs=17*1500.000000,
vegin%ejmax=17*0.000000,
vegin%ekc=17*59430.000000,
vegin%eko=17*36000.000000,
vegin%extkn=17*0.001000,
vegin%frac4=6*0.000000,1.000000,0.000000,0.000000,1.000000,7*0.000000,
vegin%froot1=0.050000,0.200000,0.200000,0.200000,0.200000,0.150000,11*0.000000,
vegin%froot2=0.050000,0.200000,0.200000,0.200000,0.200000,0.150000,11*0.000000,
vegin%froot3=0.050000,0.200000,0.200000,0.200000,0.200000,0.150000,11*0.000000,
vegin%froot4=0.050000,0.200000,0.200000,0.200000,0.200000,0.150000,11*0.000000,
vegin%froot5=0.050000,0.200000,0.200000,0.200000,0.200000,0.150000,11*0.000000,
vegin%froot6=0.050000,0.200000,0.200000,0.200000,0.200000,0.150000,11*0.000000,
vegin%g0=17*0.000000,
vegin%g1=2.346064,4.114762,2.346064,4.447321,4.694803,5.248500,1.616178,2.222156,5.789377,1.616178,5.248500,5.248500,0.000000,5.248500,5.248500,5.248500,5.248500,
vegin%gswmin=6*0.010000,0.040000,0.010000,0.010000,0.040000,7*0.010000,
vegin%hc=17.000000,35.000000,15.500000,20.000000,0.600000,0.567000,0.567000,0.567000,0.550000,0.550000,0.567000,0.200000,6.017000,0.200000,0.200000,0.200000,0.200000,
vegin%lai=4.0,5.0,0.0,0.0,0.0,0.2,11*0.0,
vegin%length=0.055000,0.100000,0.040000,0.150000,0.100000,6*0.300000,0.030000,0.242000,0.030000,0.030000,0.030000,0.030000,
vegin%ratecp1=17*1.000000,
vegin%ratecp2=17*0.030000,
vegin%ratecp3=17*0.140000,
vegin%ratecs1=17*2.000000,
vegin%ratecs2=17*0.500000,
vegin%refl1=0.090000,0.090000,0.075000,0.090000,0.090000,0.110000,0.110000,0.075000,0.110000,0.110000,0.108000,0.055000,0.091000,0.238000,0.143000,0.143000,0.159000,
vegin%refl2=0.300000,0.290000,0.300000,0.290000,0.300000,0.340000,0.340000,0.320000,0.340000,0.340000,0.343000,0.190000,0.310000,0.457000,0.275000,0.275000,0.305000,
vegin%refl3=17*0.010000,
vegin%taul1=0.090000,0.090000,0.075000,0.090000,0.090000,0.110000,0.110000,0.075000,0.110000,0.110000,0.075000,0.023000,0.059000,0.039000,0.023000,0.023000,0.026000,
vegin%taul2=0.300000,0.290000,0.300000,0.290000,0.300000,0.340000,0.340000,0.320000,0.340000,0.340000,0.146000,0.198000,0.163000,0.189000,0.113000,0.113000,0.113000,
vegin%taul3=17*0.010000,
vegin%rootbeta=0.943000,0.962000,0.966000,0.961000,0.964000,0.943000,0.943000,0.943000,0.961000,0.961000,0.943000,0.975000,5*0.961000,
vegin%rp20=3.000000,0.600000,3.000000,2.200000,1.000000,1.500000,2.800000,2.500000,1.500000,1.000000,1.500000,6*1.000000,
vegin%rpcoef=17*0.083200,
vegin%rs20=11*1.000000,0.000000,1.000000,0.000000,0.000000,0.000000,0.000000,
vegin%shelrb=17*2.000000,
vegin%tmaxvj=-10.000000,-10.000000,10.000000,15.000000,13*-10.000000,
vegin%tminvj=-15.000000,-15.000000,5.000000,5.000000,13*-15.000000,
vegin%vbeta=2.000000,2.000000,2.000000,2.000000,4.000000,4.000000,4.000000,4.000000,2.000000,2.000000,4.000000,4.000000,2.000000,4.000000,4.000000,4.000000,4.000000,
vegin%vcmax=0.000040,0.000055,0.000040,0.000060,0.000040,0.000060,0.000010,0.000040,0.000080,0.000080,0.000060,0.000017,0.000001,0.000017,0.000017,0.000017,0.000017,
vegin%vegcf=9.000000,14.000000,9.000000,8.000000,5.000000,7.000000,7.000000,5.000000,7.000000,1.000000,7.000000,6*1.000000,
vegin%wai=1.000000,1.000000,1.000000,1.000000,13*0.000000,
vegin%width=0.001000,0.050000,0.001000,0.080000,0.005000,6*0.010000,0.003000,0.015000,0.001000,0.001000,0.001000,0.001000,
vegin%xalbnir=17*1.000000,
vegin%xfang=0.010000,0.100000,0.010000,0.250000,0.010000,6*-0.300000,0.100000,5*0.000000,
vegin%zr=1.800000,3.000000,2.000000,2.000000,2.500000,5*0.500000,1.800000,3.100000,3.000000,1.000000,1.000000,1.000000,1.000000,
/
```

## Unused parameters

| Name | Type | Available | Default | Description |
|-|-|-|-|-|
| vegin%cplant | real | >=0.0 | ? | Plant carbon in 1st vegetation carbon store, now split into `cplant1-3` \( (gC \cdot m^{-2}) \) |
| vegin%csoil | real | >=0.0 | ? | Soil carbon in 1st soil carbon store, now split into `csoil1-2` \( (gC \cdot m^{-2}) \) |
| vegin%dleaf | real | 0.005 – 0.4 | ? | Characteristic length of leaf (not used) \( (m) \) |
| vegin%froot | real  | 0.0 – 1.0 | ? | Root fraction, now split into `froot1-6` (not used) \( (-) \) |
| vegin%ratecp | real | 0.01 – 3.0 | ? | Plant carbon pool rate constant in 1st vegetation carbon store, now split into `ratecp1-3` \( (year^{-1}) \) |
| vegin%ratecs | real | 0.01 – 3.0 | ? | Soil carbon pool rate constant in 1st soil carbon store, now split into `ratecs1-3` \( (year^{-1}) \) |
| vegin%refl | real | 0.0 – 0.5 | ? | Leaf reflectance in 1st radiation band, now split into `refl1-3` \( (-) \) |
| vegin%taul1 | real | 0.0 – 0.3 | ? | Leaf transmittance in 1st radiation band, now split into `tau1-3` \( (-) \) |
| vegin%wai | real | 0.0 – 5.0 | 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Wood area index (stem+branches+twigs) (not used) \( (-) \) |
| vegin%xalbnir | real | 0.0 – 1.5 | 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 | ? (not used) |

