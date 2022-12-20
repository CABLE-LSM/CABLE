# pft_params.nml options

| Namelist variable| Type | Available values | Default values | Description |
|------------------|------|------------------|----------------|-------------|
| vegin%a1gs | real | >=0.0 | 9, 9, 9, 9, 9, 9, 4, 9, 9, 4, 9, 9, 9, 9, 9, 9, 9 | a1 parameter in stomatal conductance model. Represents the sensitivity of stomatal conductance to the assimilation rate \( (-) \) |
| vegin%alpha | real | >=0.0 | 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.05, 0.2, 0.2, 0.05, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 | Initial slope of J-Q response curve \( (mol (electrons) \cdot mol^{-1} (photons) (C3) \cdot mol (CO2) \cdot mol^{-1} (photons) (C4)) \) |
| vegin%canst1 | real | 0.05 – 0.15 | 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 | Maximum intercepted water by canopy \( (mm \cdot LAI^{-1}) \) |
| vegin%cfrd | real | >=0.0 | 0.015, 0.015, 0.015, 0.015, 0.015, 0.025, 0.015, 0.015, 0.025, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015 | Ratio of day respiration to vcmax \( (-) \) |
| vegin%clitt | real | >=0.0 | 20, 6, 10, 13, 2, 2, 0.3, 0.3, 0, 0, 2, 2, 0, 0, 0, 0, 0 | Leaf litter (alters resistance to soil evaporation) \( (tC \cdot ha^{-1}) \) |
| vegin%conkc0 | real | >=0.0 | 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302, 0.000302 | Michaelis-menton constant for carboxylase \( (bar) \) |
| vegin%conko0 | real | >=0.0 | 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256, 0.256 | Michaelis-menton constant for oxygenase \( (bar) \) |
| vegin%convex | real | >=0.0 | 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.8, 0.01, 0.01, 0.8, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 | Convexity of J-Q response curve \( (-) )/ |
| vegin%cplant1 | real | >=0.0 | 200, 300, 200, 300, 159, 250, 250, 250, 150, 150, 250, 1, 0.1, 0, 1, 1, 0 | Plant carbon in 1st vegetation carbon store \( (gC \cdot m^{-2}) \) |
| vegin%cplant2 | real | >=0.0 | 10217, 16833, 5967, 12000, 5000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Plant carbon in 2nd vegetation carbon store \( (gC \cdot m^{-2}) \) |
| vegin%cplant3 | real | >=0.0 | 876, 1443, 511, 1029, 500, 500, 500, 500, 607, 607, 500, 1, 0.1, 0, 1, 1, 0 | Plant carbon in 3rd vegetation carbon store \( (gC \cdot m^{-2}) \) |
| vegin%csoil1 | real | >=0.0 | 184, 303, 107, 216, 100, 275, 275, 275, 149, 149, 275, 1, 0.1, 1, 1, 1, 1 | Soil carbon in 1st soil carbon store \( (gC \cdot m^{-2}) \) |
| vegin%csoil2 | real | >=0.0 | 367, 606, 214, 432, 250, 314, 314, 314, 300, 300, 314, 1, 0.1, 1, 1, 1, 1 | Soil carbon in 2nd soil carbon store \( (gC \cdot m^{-2}) \) |
| vegin%d0gs | real | >=0.0 | 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500 | d0 in stomatal conductance model \( (kPa) \) |
| vegin%desc | character | | 'Evergreen Needleleaf', 'Evergreen Broadleaf', 'Deciduous Needleleaf', 'Deciduous Broadleaf', 'Shrub', 'C3 Grassland', 'C4 Grassland',  'Tundra', 'C3 Cropland', 'C4 Cropland', 'Wetland', 'empty', 'empty', 'Barren', 'Urban', 'Lakes', 'Ice' | Description of plant functional type |
| vegin%ejmax | real | >=0.0 | 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Maximum potential electron transp rate top leaf \( (mol \cdot m^{-2} \cdot s^{-1}) \) |
| vegin%ekc | real | >=0.0 | 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430, 59430 | Activation energy for caroxylagse \( (J \cdot mol^{-1}) \) |
| vegin%eko | real | >=0.0 | 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000 | Acvtivation enegery for oxygenase \( (J \cdot mol^{-1}) \) |
| vegin%extkn | real | 0.0 – 10.0 | 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001 | Extinction coeficient for vertical profile of N \( (-) \) |
| vegin%frac4 | real | 0.0 – 1.0 | 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0| Fraction of c4 plants \( (-) \) |
| vegin%froot1 | real | 0.0 – 1.0 | 0.05, 0.2, 0.2, 0.2, 0.2, 0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Fraction of root in soil layer 1 \( (-) \) |
| vegin%froot2 | real | 0.0 – 1.0 | 0.05, 0.2, 0.2, 0.2, 0.2, 0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Fraction of root in soil layer 2 \( (-) \) |
| vegin%froot3 | real | 0.0 – 1.0 | 0.05, 0.2, 0.2, 0.2, 0.2, 0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Fraction of root in soil layer 3 \( (-) \) |
| vegin%froot4 | real | 0.0 – 1.0 | 0.05, 0.2, 0.2, 0.2, 0.2, 0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Fraction of root in soil layer 4 \( (-) \) |
| vegin%froot5 | real | 0.0 – 1.0 | 0.05, 0.2, 0.2, 0.2, 0.2, 0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Fraction of root in soil layer 5 \( (-) \) |
| vegin%froot6 | real | 0.0 – 1.0 | 0.05, 0.2, 0.2, 0.2, 0.2, 0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Fraction of root in soil layer 6 \( (-) \) |
| vegin%g0 | real | >=0.0 | 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Belinda's stomatal model intercept. Residual stomatal conductance as net assimilation rate reaches zero \( (mol \cdot m^{-2} \cdot s^{-1}) \) |
| vegin%g1 | real | >=0.0 | 2.346064, 4.114762, 2.346064, 4.447321, 4.694803, 5.2485, 1.616178, 2.222156, 5.789377, 1.616178, 5.2485, 5.2485, 0, 5.2485, 5.2485, 5.2485, 5.2485 | Belinda's stomatal model slope. Sensitivity of stomatal conductance to the assimilation rate \( (kPa) \) |
| vegin%gswmin | real | >=0.0 | 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.04, 0.01, 0.01, 0.04, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 | Minimal stomatal conductance \( (mol \cdot m^{-2} \cdot s^{-1}) \) |
| vegin%hc | real | >=0.0 | 17, 35, 15.5, 20, 0.6, 0.567, 0.567, 0.567, 0.55, 0.55, 0.567, 0.2, 6.017, 0.2, 0.2, 0.2, 0.2 | Roughness height of canopy (veg - snow) \( (m) \) |
| vegin%lai | real | >=0.0 | 4, 5, 0, 0, 0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Leaf area index of each plant functional type \( (-) or (m^{2} \cdot m^{-2} ) \) |
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
| vegin%rpcoef | real | 0.05 – 1.5 | 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832, 0.0832 | Temperature coeficient non-leaf plant respiration \( (C^{-1}) \) |
| vegin%rs20 | real | 0.0 – 10.0 | 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0 | Soil respiration at 20 C \( (mol \cdot m^{-2} \cdot s^{-1}] ) \) |
| vegin%shelrb | real | >=0.0 | 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 | Sheltering factor \( (-) \) |
| vegin%taul1 | real | 0.0 – 0.3 | 0.09, 0.09, 0.075, 0.09, 0.09, 0.11, 0.11, 0.075, 0.11, 0.11, 0.075, 0.023, 0.059, 0.039, 0.023, 0.023, 0.026 | Leaf transmittance in 1st radiation band \( (-) \) |
| vegin%taul2 | real | 0.0 – 0.3 | 0.3, 0.29, 0.3, 0.29, 0.3, 0.34, 0.34, 0.32, 0.34, 0.34, 0.146, 0.198, 0.163, 0.189, 0.113, 0.113, 0.113 | Leaf transmittance in 2nd radiation band \( (-) \) |
| vegin%taul3 | real | 0.0 – 0.3 | 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 | Leaf transmittance in 3rd radiation band \( (-) \) |
| vegin%tmaxvj | real | >=0.0 | -10, -10, 10, 15, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10 | Maximum temperature of the start of photosynthesis \( (^{\circ}C) \) |
| vegin%tminvj | real | >=0.0 | -15, -15, 5, 5, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15 | Minimum temperature of the start of photosynthesis \( (^{\circ}C) \) |
| vegin%vbeta | real | >=0.0 | 2, 2, 2, 2, 4, 4, 4, 4, 2, 2, 4, 4, 2, 4, 4, 4, 4 | Stomatal sensitivity to soil water \( (-) \) |
| vegin%vcmax | real | >=0.0 | 4E-05, 5.5E-05, 4E-05, 6E-05, 4E-05, 6E-05, 1E-05, 4E-05, 8E-05, 8E-05, 6E-05, 1.7E-05, 1E-06, 1.7E-05, 1.7E-05, 1.7E-05, 1.7E-05 | Maximum RuBP carboxylation rate top leaf \( ( mol \cdot m^{-2} \cdot s^{-1} ) \) |
| vegin%vegcf | real | 0.0 – 10.0 | 9, 14, 9, 8, 5, 7, 7, 5, 7, 1, 7, 1, 1, 1, 1, 1, 1 | Scaling on soil respiration (place-holder scheme) \( (-) \) |
| vegin%width | real | >=0.0 | 0.001, 0.05, 0.001, 0.08, 0.005, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.003, 0.015, 0.001, 0.001, 0.001, 0.001 | Leaf width \( (m) \) |
| vegin%xfang | real | -1.0 – 0.5 | 0.01, 0.1, 0.01, 0.25, 0.01, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, 0.1, 0, 0, 0, 0, 0| Leaf angle \( (-) \) |
| vegin%zr | real | >=0.0 | 1.8, 3, 2, 2, 2.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.8, 3.1, 3, 1, 1, 1, 1 | Maximum rooting depth \( (cm) \) |

## Old and/or unused parameters

| Name | Type | Available | Default | Description |
|-|-|-|-|-|
| deciduous | logical  |  |  | Flag used for phenology fix |
| dleaf | real |  |  | Chararacteristc legnth of leaf \( (m) \) |
| froot | real  |  |  |  |
| iLU | integer |  |  | Land use type |
| isoil | integer |  | 1 – 30 | Soil type |
| iveg | integer |  | 1 – 30 | Vegetation type |
| meth | integer |  |  | Method for calculation of canopy fluxes and temperature |
| toptvj | real |  |  | Optimum temperature of the start of photosynthesis \( ^{/circ}C \) |
| vegin%wai | real | 0.0 – 1.0 | 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 | Wood area index (stem+branches+twigs) (not used) \( (-) \) |
| vegin%xalbnir | real | >=0.0 | 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 | Not used |
| vlai  | real |  |  | Leaf area index |
| vlaimax | real |  |  |  |

