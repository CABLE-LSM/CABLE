# List of inputs

## Vegetation parameters file def_veg_params.txt

| Variable      | CABLE variable | Comment                                                                                      |
|---------------|----------------|----------------------------------------------------------------------------------------------|
| pft_number    | jveg           |                                                                                              |
| kind          | vegtypetmp     |                                                                                              |
| name          | vegnametmp     |                                                                                              |
| canopy_height | vegin%hc       | Canopy height                                                                                |
| leaf_angle    | vegin%xfang    | Leaf angle. Not used.                                                                        |
| leaf_width    | vegin%width    | Leaf width. Used for dlength: the characteristic length of leaf                              |
| leaf_length   | vegin%length   | Leaf width. Used for dlength: the characteristic length of leaf                              |
| c4_fraction   | vegin%frac4    | Fraction of vegetation that are c4 plants.                                                   |
| rholeaf_vis   | vegin%refl     | leaf reflectivity                                                                            |
| rholeaf_nir   | vegin%refl     | leaf reflectivity                                                                            |
| rholeaf_therm | vegin%refl     | leaf reflectivity                                                                            |
| rhowood_vis   | x              | Not used                                                                                     |
| rhowood_nir   | x              | Not used                                                                                     |
| rhowood_therm | x              | Not used                                                                                     |
| tauleaf_vis   | vegin%taul     | Leaf transmisivity                                                                           |
| tauleaf_nir   | vegin%taul     | Leaf transmisivity                                                                           |
| tauleaf_therm | vegin%taul     | Leaf transmisivity                                                                           |
| tauwood_vis   | x              | Not used                                                                                     |
| tauwood_nir   | x              | Not used                                                                                     |
| tauwood_therm | x              | Not used                                                                                     |
| rhosoil_vis   | x              | Not used                                                                                     |
| rhosoil_nir   | x              | Not used                                                                                     |
| rhosoil_therm | x              | Not used                                                                                     |
| xalbnir       | xalbnir        | Modifier for surface albedo in the near infra-red band.                                      |
| laimax        | x              | Not used                                                                                     |
| woodai        | vegin%wai      | Not used                                                                                     |
| canst1        | vegin%canst1   | BATS-type canopy saturation proportional to LAI                                              |
| shelrb        | vegin%shelrb   | Sheltering factor (dimensionless) Something to do with vegetation boundary-layer conductance |
| vegcf         | vegin%vegcf    | Kdcorbin 08/10. Used in the calculation of plant and soil respiration                        |
| extkn         | vegin%extkn    | Extinction coefficient for vertical                                                          |
| vcmax         | vegin%vcmax    | Max RuBP carboxylation rate top leaf (mol/m2/s)                                              |
| rp20          | vegin%rp20     | Plant respiration coefficient at 20 degrees C                                                |
| rpcoeff       | vegin%rpcoef   | Temperature coef nonleaf plant respiration (1/C)                                             |
| rs20          | vegin%rs20     | Soil respiration at 20 C [mol m-2 s-1]                                                       |
| tvjmin        | Vegin%tminvj   | Min temperature of the start of photosynthesis                                               |
| tvjmax        | vegin%tmaxvj   | Max temperature of the start of photosynthesis                                               |
| vbeta         | vegin%vbeta    | Stomatal sensitivity to soil water                                                           |
| rootbeta      | vegin%rootbeta | Estimaties vertical root mass distribution (froot)                                           |
| rootdepth     | x              | Not used                                                                                     |
| clitt         | x              | Not used                                                                                     |
| leaf_pool     | vegin%cplant   | Plant carbon (g C/m2)                                                                        |
| wood_pool     | vegin%cplant   |                                                                                              |
| root_pool     | vegin%cplant   |                                                                                              |
| soilfast_pool | vegin%csoil    | Soil carbon (g C/m2)                                                                         |
| soilslow_pool | vegin%csoil    |                                                                                              |
| leaf_rate     | vegin%ratecp   | Plant carbon rate constant (1/year)                                                          |
| wood_rate     | vegin%ratecp   |                                                                                              |
| root_rate     | vegin%ratecp   |                                                                                              |
| soilfast_rate | vegin%ratecs   | Soil carbon rate constant (1/year)                                                           |
| soilslow_rate | vegin%ratecs   |                                                                                              |
| a1            | x              |                                                                                              |
| d0            | x              |                                                                                              |
| alpha         | x              |                                                                                              |
| convex        | x              |                                                                                              |
| cfrd          | x              |                                                                                              |
| gswmin        | x              |                                                                                              |
| conkc0        | x              |                                                                                              |
| conko0        | x              |                                                                                              |
| ekc           | x              |                                                                                              |
| eko           | x              |                                                                                              |
| g0c3or4       | x              |                                                                                              |
| g1c3or4       | x              |                                                                                              |
