# pft_params.nml options

| Namelist variable| Type | Available values | Default values | Description |
|------------------|------|------------------|----------------|-------------|
| vegin%a1gs | real | | 9,9,9,9,9,9,4,9,9,4,9,9,9,9,9,9,9| a1 parameter in stomatal conductance model |
| vegin%alpha | real | | 0.2,0.2,0.2,0.2,0.2,0.2,0.05,0.2,0.2,0.05,0.2,0.2,0.2,0.2,0.2,0.2,0.2| Initial slope of J-Q response curve |
| vegin%canst1 | real | | 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1| Maximum intercepted water by canopy /( (mm \cdot LAI^{-1}) /) |
| vegin%cfrd | real | | 0.015,0.015,0.015,0.015,0.015,0.025,0.015,0.015,0.025,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015| Ratio of day respiration to vcmax |
| vegin%clitt | real | | 20,6,10,13,2,2,0.3,0.3,0,0,2,2,0,0,0,0,0|  |
| vegin%conkc0 | real | | 0.000302,0.000302,0.000302,0.000302,0.000302,0.000302,0.000302,0.000302,0.000302,0.000302,0.000302,0.000302,0.000302,0.000302,0.000302,0.000302,0.000302| Michaelis-menton constant for carboxylase |
| vegin%conko0 | real | | 0.256,0.256,0.256,0.256,0.256,0.256,0.256,0.256,0.256,0.256,0.256,0.256,0.256,0.256,0.256,0.256,0.256| Michaelis-menton constant for oxygenase |
| vegin%convex | real | | 0.01,0.01,0.01,0.01,0.01,0.01,0.8,0.01,0.01,0.8,0.01,0.01,0.01,0.01,0.01,0.01,0.01| Convexity of J-Q response curve |
| vegin%cplant1 | real | | 200,300,200,300,159,250,250,250,150,150,250,1,0.1,0,1,1,0| |
| vegin%cplant2 | real | | 10217,16833,5967,12000,5000,0,0,0,0,0,0,0,0,0,0,0,0| |
| vegin%cplant3 | real | | 876,1443,511,1029,500,500,500,500,607,607,500,1,0.1,0,1,1,0| |
| vegin%csoil1 | real | | 184,303,107,216,100,275,275,275,149,149,275,1,0.1,1,1,1,1| |
| vegin%csoil2 | real | | 367,606,214,432,250,314,314,314,300,300,314,1,0.1,1,1,1,1| d0 in stomatal conductance model |
| vegin%d0gs | real | | 1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500| |
| vegin%desc | character | | 'Evergreen Needleleaf','Evergreen Broadleaf','Deciduous Needleleaf','Deciduous Broadleaf','Shrub','C3 Grassland','C4 Grassland','Tundra','C3 Cropland','C4 Cropland','Wetland','empty','empty','Barren','Urban','Lakes','Ice'| |
| vegin%ejmax | real | | 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0|Maximum potential electron transp rate top leaf /( (mol \cdot m^{-2} \cdot s^{-1}) /) |
| vegin%ekc | real | | 59430,59430,59430,59430,59430,59430,59430,59430,59430,59430,59430,59430,59430,59430,59430,59430,59430| Activation energy for caroxylagse |
| vegin%eko | real | | 36000,36000,36000,36000,36000,36000,36000,36000,36000,36000,36000,36000,36000,36000,36000,36000,36000| Acvtivation enegery for oxygenase |
| vegin%extkn | real | | 0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001| Extinction coeficient for vertical |
| vegin%frac4 | real | | 0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0| Fraction of c4 plants |
| vegin%froot1 | real | | 0.05,0.2,0.2,0.2,0.2,0.15,0,0,0,0,0,0,0,0,0,0,0 | Fraction of root in soil layer 1 |
| vegin%froot2 | real | | 0.05,0.2,0.2,0.2,0.2,0.15,0,0,0,0,0,0,0,0,0,0,0 | Fraction of root in soil layer 2 |
| vegin%froot3 | real | | 0.05,0.2,0.2,0.2,0.2,0.15,0,0,0,0,0,0,0,0,0,0,0 | Fraction of root in soil layer 3 |
| vegin%froot4 | real | | 0.05,0.2,0.2,0.2,0.2,0.15,0,0,0,0,0,0,0,0,0,0,0 | Fraction of root in soil layer 4 |
| vegin%froot5 | real | | 0.05,0.2,0.2,0.2,0.2,0.15,0,0,0,0,0,0,0,0,0,0,0 | Fraction of root in soil layer 5 |
| vegin%froot6 | real | | 0.05,0.2,0.2,0.2,0.2,0.15,0,0,0,0,0,0,0,0,0,0,0 | Fraction of root in soil layer 6 |
| vegin%g0 | real | | 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0| Belinda's stomatal model intercept |
| vegin%g1 | real | | 2.346064,4.114762,2.346064,4.447321,4.694803,5.2485,1.616178,2.222156,5.789377,1.616178,5.2485,5.2485,0,5.2485,5.2485,5.2485,5.2485| Belinda's stomatal model slope |
| vegin%gswmin | real | | 0.01,0.01,0.01,0.01,0.01,0.01,0.04,0.01,0.01,0.04,0.01,0.01,0.01,0.01,0.01,0.01,0.01| Minimal stomatal conductance |
| vegin%hc | real | | 17,35,15.5,20,0.6,0.567,0.567,0.567,0.55,0.55,0.567,0.2,6.017,0.2,0.2,0.2,0.2| Roughness height of canopy (veg - snow) |
| vegin%lai | real | | 4,5,0,0,0,0.2,0,0,0,0,0,0,0,0,0,0,0| |
| vegin%length | real | | 0.055,0.1,0.04,0.15,0.1,0.3,0.3,0.3,0.3,0.3,0.3,0.03,0.242,0.03,0.03,0.03,0.03| |
| vegin%ratecp1 | real | | 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1| |
| vegin%ratecp2 | real | | 0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03| |
| vegin%ratecp3 | real | | 0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14| |
| vegin%ratecs1 | real | | 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2| |
| vegin%ratecs2 | real | | 0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5| |
| vegin%refl1 | real | | 0.09,0.09,0.075,0.09,0.09,0.11,0.11,0.075,0.11,0.11,0.108,0.055,0.091,0.238,0.143,0.143,0.159| |
| vegin%refl2 | real | | 0.3,0.29,0.3,0.29,0.3,0.34,0.34,0.32,0.34,0.34,0.343,0.19,0.31,0.457,0.275,0.275,0.305| |
| vegin%refl3 | real | | 0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01| |
| vegin%rootbeta | real | | 0.943,0.962,0.966,0.961,0.964,0.943,0.943,0.943,0.961,0.961,0.943,0.975,0.961,0.961,0.961,0.961,0.961| |
| vegin%rp20 | real | | 3,0.6,3,2.2,1,1.5,2.8,2.5,1.5,1,1.5,1,1,1,1,1,1| Plant respiration coefficient at 20 C |
| vegin%rpcoef | real | | 0.0832,0.0832,0.0832,0.0832,0.0832,0.0832,0.0832,0.0832,0.0832,0.0832,0.0832,0.0832,0.0832,0.0832,0.0832,0.0832,0.0832| Temperature coef nonleaf plant respiration (1/C) |
| vegin%rs20 | real | | 1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0| Soil respiration at 20 C \( (mol \cdot m^{-2} \cdot s^{-1}] ) /) |
| vegin%shelrb | real | | 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2| Sheltering factor (dimensionless) |
| vegin%taul1 | real | | 0.09,0.09,0.075,0.09,0.09,0.11,0.11,0.075,0.11,0.11,0.075,0.023,0.059,0.039,0.023,0.023,0.026| |
| vegin%taul2 | real | | 0.3,0.29,0.3,0.29,0.3,0.34,0.34,0.32,0.34,0.34,0.146,0.198,0.163,0.189,0.113,0.113,0.113| |
| vegin%taul3 | real | | 0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01| |
| vegin%tmaxvj | real | | -10,-10,10,15,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10| Maximum temperature of the start of photosynthesis |
| vegin%tminvj | real | | -15,-15,5,5,-15,-15,-15,-15,-15,-15,-15,-15,-15,-15,-15,-15,-15| Minimum temperature of the start of photosynthesis |
| vegin%vbeta | real | | 2,2,2,2,4,4,4,4,2,2,4,4,2,4,4,4,4| Maximum RuBP carboxylation rate top leaf /( ( mol \cdot m^{-2} \cdot s^{-1} ) /) |
| vegin%vcmax | real | | 4E-05,5.5E-05,4E-05,6E-05,4E-05,6E-05,1E-05,4E-05,8E-05,8E-05,6E-05,1.7E-05,1E-06,1.7E-05,1.7E-05,1.7E-05,1.7E-05| |
| vegin%vegcf | real | | 9,14,9,8,5,7,7,5,7,1,7,1,1,1,1,1,1| Wood area index (stem+branches+twigs) |
| vegin%wai | real | | 1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0| |
| vegin%width | real | | 0.001,0.05,0.001,0.08,0.005,0.01,0.01,0.01,0.01,0.01,0.01,0.003,0.015,0.001,0.001,0.001,0.001| |
| vegin%xalbnir | real | | 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1| |
| vegin%xfang | real | | 0.01,0.1,0.01,0.25,0.01,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,0.1,0,0,0,0,0| Leaf angle |
| vegin%zr | real | | 1.8,3,2,2,2.5,0.5,0.5,0.5,0.5,0.5,1.8,3.1,3,1,1,1,1 |
| vegin%meth | integer |  |  | Method for calculation of canopy fluxes and temperature |
| vegin%iveg | integer |  |  | Vegetation type |
| vegin%iLU | integer |  |  | Land use type |
| vegin%dleaf | real |  |  | Chararacteristc legnth of leaf (m) |
| vegin%vlai  | real |  |  | Leaf area index |
| vegin%toptvj | real |  |  | Optimum temperature of the start of photosynthesis |
| vegin%vlaimax | real |  |  |  |
| vegin%deciduous | logical  |  |  | Flag used for phenology fix |
| vegin%froot | real  |  |  |  |
