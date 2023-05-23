# pft_params.nml options

## Example pft_params.nml file

!!! Note

    The asterisk symbol does not mean arithmetic multiplication in the namelist file.
    It means that the following value is repeated multiple times.

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

## Descriptive table of parameters

| Namelist variable| Type | Available values | Default values | Description |
|------------------|------|------------------|----------------|-------------|
| vegin%a1gs | real | >=0.0 | 9, 9, 9, 9, 9, 9, 4, 9, 9, 4, 9, 9, 9, 9, 9, 9, 9 | a1 parameter in stomatal conductance model. Represents the sensitivity of stomatal conductance to the assimilation rate \( (-) \) |
| vegin%alpha | real | >=0.0 | 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.05, 0.2, 0.2, 0.05, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 | Initial slope of J-Q response curve \( (mol (electrons) \cdot mol^{-1} (photons) (C3) \cdot mol (CO2) \cdot mol^{-1} (photons) (C4)) \) |
| vegin%canst1 | real | 0.05 – 0.15 | 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 | Maximum intercepted water by canopy \( (mm \cdot LAI^{-1}) \) |
| vegin%cfrd | real | >=0.0 | 0.015, 0.015, 0.015, 0.015, 0.015, 0.025, 0.015, 0.015, 0.025, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015 | Ratio of day respiration to vcmax \( (-) \) |
| vegin%clitt | real | >=0.0 | 20, 6, 10, 13, 2, 2, 0.3, 0.3, 0, 0, 2, 2, 0, 0, 0, 0, 0 | Leaf litter (alters resistance to soil evaporation) \( (tC \cdot ha^{-1}) \) |
| vegin%conkc0 | real | >=0.0 | 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302 | Michaelis-Menton constant for carboxylase \( (bar) \) |
| vegin%conko0 | real | >=0.0 | 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256 | Michaelis-Menton constant for oxygenase \( (bar) \) |
| vegin%convex | real | >=0.0 | 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.8, 0.01, 0.01, 0.8, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 | Convexity of J-Q response curve \( (-) \) |
| vegin%cplant1 | real | >=0.0 | 200, 300, 200, 300, 159, 250, 250, 250, 150, 150, 250, 1, 0.1, 0, 1, 1, 0 | Plant carbon in 1st vegetation carbon store \( (gC \cdot m^{-2}) \) |
| vegin%cplant2 | real | >=0.0 | 10217, 16833, 5967, 12000, 5000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Plant carbon in 2nd vegetation carbon store \( (gC \cdot m^{-2}) \) |
| vegin%cplant3 | real | >=0.0 | 876, 1443, 511, 1029, 500, 500, 500, 500, 607, 607, 500, 1, 0.1, 0, 1, 1, 0 | Plant carbon in 3rd vegetation carbon store \( (gC \cdot m^{-2}) \) |
| vegin%csoil1 | real | >=0.0 | 184, 303, 107, 216, 100, 275, 275, 275, 149, 149, 275, 1, 0.1, 1, 1, 1, 1 | Soil carbon in 1st soil carbon store \( (gC \cdot m^{-2}) \) |
| vegin%csoil2 | real | >=0.0 | 367, 606, 214, 432, 250, 314, 314, 314, 300, 300, 314, 1, 0.1, 1, 1, 1, 1 | Soil carbon in 2nd soil carbon store \( (gC \cdot m^{-2}) \) |
| vegin%d0gs | real | >=0.0 | 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500 | d0 in stomatal conductance model \( (kPa) \) |
| vegin%desc | character | | 'Evergreen Needleleaf', 'Evergreen Broadleaf', 'Deciduous Needleleaf', 'Deciduous Broadleaf', 'Shrub', 'C3 Grassland', 'C4 Grassland',  'Tundra', 'C3 Cropland', 'C4 Cropland', 'Wetland', 'empty', 'empty', 'Barren', 'Urban', 'Lakes', 'Ice' | Description of plant functional type |
| vegin%ejmax | real | 0 – 3.0E-4 | 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Maximum potential electron transp rate top leaf \( (mol \cdot m^{-2} \cdot s^{-1}) \) |
| vegin%ekc | real | >=0.0 | 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430 | Activation energy for carboxylase \( (J \cdot mol^{-1}) \) |
| vegin%eko | real | >=0.0 | 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000 | Activation energy for oxygenase \( (J \cdot mol^{-1}) \) |
| vegin%extkn | real | 0.0 – 10.0 | 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001 | Extinction coeficient for vertical profile of N \( (-) \) |
| vegin%frac4 | real | 0.0 – 1.0 | 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0| Fraction of c4 plants \( (-) \) |
| vegin%froot1 | real | 0.0 – 1.0 | 0.05, 0.2, 0.2, 0.2, 0.2, 0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Fraction of root in soil layer 1 \( (-) \) |
| vegin%froot2 | real | 0.0 – 1.0 | 0.05, 0.2, 0.2, 0.2, 0.2, 0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Fraction of root in soil layer 2 \( (-) \) |
| vegin%froot3 | real | 0.0 – 1.0 | 0.05, 0.2, 0.2, 0.2, 0.2, 0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Fraction of root in soil layer 3 \( (-) \) |
| vegin%froot4 | real | 0.0 – 1.0 | 0.05, 0.2, 0.2, 0.2, 0.2, 0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Fraction of root in soil layer 4 \( (-) \) |
| vegin%froot5 | real | 0.0 – 1.0 | 0.05, 0.2, 0.2, 0.2, 0.2, 0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Fraction of root in soil layer 5 \( (-) \) |
| vegin%froot6 | real | 0.0 – 1.0 | 0.05, 0.2, 0.2, 0.2, 0.2, 0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Fraction of root in soil layer 6 \( (-) \) |
| vegin%g0 | real | -0.5 – 0.5 | 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Belinda's stomatal model intercept. Residual stomatal conductance as net assimilation rate reaches zero \( (mol \cdot m^{-2} \cdot s^{-1}) \) |
| vegin%g1 | real | 0.0 – 20.0 | 2.346064, 4.114762, 2.346064, 4.447321, 4.694803, 5.2485, 1.616178, 2.222156, 5.789377, 1.616178, 5.2485, 5.2485, 0, 5.2485, 5.2485, 5.2485, 5.2485 | Belinda's stomatal model slope. Sensitivity of stomatal conductance to the assimilation rate \( (kPa) \) |
| vegin%gswmin | real | >=0.0 | 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.04, 0.01, 0.01, 0.04, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 | Minimal stomatal conductance \( (mol \cdot m^{-2} \cdot s^{-1}) \) |
| vegin%hc | real | 0.0 – 100.0 | 17, 35, 15.5, 20, 0.6, 0.567, 0.567, 0.567, 0.55, 0.55, 0.567, 0.2, 6.017, 0.2, 0.2, 0.2, 0.2 | Roughness height of canopy (veg - snow) \( (m) \) |
| vegin%lai | real | 0.0 – 8.0 | 4, 5, 0, 0, 0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Leaf area index of each plant functional type \( (-) or (m^{2} \cdot m^{-2} ) \) |
| vegin%length | real | >=0.0 | 0.055, 0.1, 0.04, 0.15, 0.1, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.03, 0.242, 0.03, 0.03, 0.03, 0.03 | Leaf length \( (m) \) |
| vegin%ratecp1 | real | 0.01 – 3.0 | 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 | Plant carbon pool rate constant in 1st vegetation carbon store \( (year^{-1}) \) |
| vegin%ratecp2 | real | 0.01 – 3.0 | 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03 | Plant carbon pool rate constant in 2nd vegetation carbon store  \( (year^{-1}) \) |
| vegin%ratecp3 | real | 0.01 – 3.0 | 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14 | Plant carbon pool rate constant in 3rd vegetation carbon store  \( (year^{-1}) \) |
| vegin%ratecs1 | real | 0.01 – 3.0 | 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 | Soil carbon pool rate constant in 1st soil carbon store \( (year^{-1}) \) |
| vegin%ratecs2 | real | 0.01 – 3.0 | 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 | Soil carbon pool rate constant in 2nd soil carbon store \( (year^{-1}) \) |
| vegin%refl1 | real | 0.0 – 0.5 | 0.09, 0.09, 0.075, 0.09, 0.09, 0.11, 0.11, 0.075, 0.11, 0.11, 0.108, 0.055, 0.091, 0.238, 0.143, 0.143, 0.159 | Leaf reflectance in 1st radiation band \( (-) \) |
| vegin%refl2 | real | 0.0 – 0.5 | 0.3, 0.29, 0.3, 0.29, 0.3, 0.34, 0.34, 0.32, 0.34, 0.34, 0.343, 0.19, 0.31, 0.457, 0.275, 0.275, 0.305 | Leaf reflectance in 2nd radiation band \( (-) \) |
| vegin%refl3 | real | 0.0 – 0.5 | 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 | Leaf reflectance in 3rd radiation band \( (-) \) |
| vegin%rootbeta | real | 0.0 – 1.0 | 0.943, 0.962, 0.966, 0.961, 0.964, 0.943, 0.943, 0.943, 0.961, 0.961, 0.943, 0.975, 0.961, 0.961, 0.961, 0.961, 0.961 | Beta parameter to calculate froot (Jackson et al. 1996) \( (-) \) |
| vegin%rp20 | real | 0.0 – 10.0 | 3, 0.6, 3, 2.2, 1, 1.5, 2.8, 2.5, 1.5, 1, 1.5, 1, 1, 1, 1, 1, 1 | Plant respiration coefficient at 20 C \( (-) \) |
| vegin%rpcoef | real | 0.05 – 1.5 | 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832 | Temperature coefficient non-leaf plant respiration \( (^{\circ}C^{-1}) \) |
| vegin%rs20 | real | 0.0 – 10.0 | 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0 | Soil respiration at 20 C \( (mol \cdot m^{-2} \cdot s^{-1} ) \) |
| vegin%shelrb | real | 1.0 – 3.0 | 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 | Sheltering factor \( (-) \) |
| vegin%taul1 | real | 0.0 – 0.3 | 0.09, 0.09, 0.075, 0.09, 0.09, 0.11, 0.11, 0.075, 0.11, 0.11, 0.075, 0.023, 0.059, 0.039, 0.023, 0.023, 0.026 | Leaf transmittance in 1st radiation band \( (-) \) |
| vegin%taul2 | real | 0.0 – 0.3 | 0.3, 0.29, 0.3, 0.29, 0.3, 0.34, 0.34, 0.32, 0.34, 0.34, 0.146, 0.198, 0.163, 0.189, 0.113, 0.113, 0.113 | Leaf transmittance in 2nd radiation band \( (-) \) |
| vegin%taul3 | real | 0.0 – 0.3 | 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 | Leaf transmittance in 3rd radiation band \( (-) \) |
| vegin%tmaxvj | real | -15.0 – 30.0 | -10, -10, 10, 15, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10 | Maximum temperature of the start of photosynthesis \( (^{\circ}C) \) |
| vegin%tminvj | real | -20.0 – 15.0 | -15, -15, 5, 5, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15 | Minimum temperature of the start of photosynthesis \( (^{\circ}C) \) |
| vegin%vbeta | real | -999999 – 999999 | 2, 2, 2, 2, 4, 4, 4, 4, 2, 2, 4, 4, 2, 4, 4, 4, 4 | Stomatal sensitivity to soil water \( (-) \) |
| vegin%vcmax | real | 5.0E-6 – 1.5E-4 | 4E-05, 5.5E-05, 4E-05, 6E-05, 4E-05, 6E-05, 1E-05, 4E-05, 8E-05, 8E-05, 6E-05, 1.7E-05, 1E-06, 1.7E-05, 1.7E-05, 1.7E-05, 1.7E-05 | Maximum RuBP carboxylation rate top leaf \( ( mol \cdot m^{-2} \cdot s^{-1} ) \) |
| vegin%vegcf | real | 0.0 – 100.0 | 9, 14, 9, 8, 5, 7, 7, 5, 7, 1, 7, 1, 1, 1, 1, 1, 1 | Scaling on soil respiration (place-holder scheme) \( (-) \) |
| vegin%width | real | >=0.0 | 0.001, 0.05, 0.001, 0.08, 0.005, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.003, 0.015, 0.001, 0.001, 0.001, 0.001 | Leaf width \( (m) \) |
| vegin%xfang | real | -1.0 – 0.5 | 0.01, 0.1, 0.01, 0.25, 0.01, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, 0.1, 0, 0, 0, 0, 0| Leaf angle \( (-) \) |
| vegin%zr | real | >=0.0 | 1.8, 3, 2, 2, 2.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.8, 3.1, 3, 1, 1, 1, 1 | Maximum rooting depth \( (cm) \) |

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

