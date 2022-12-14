# pft_params.nml options

|Namelist variable| Type | Available values | Default values | Description |
|-|-|-|-|-|
|          vegin%meth,    | integer |  |  | method for calculation of canopy fluxes and temp. |
|          vegin%iveg , | integer |  |  | vegetation type |
|          vegin%iLU | integer |  |  | land use type |
|          vegin%canst1,  | real |  |  | max intercepted water by canopy (mm/LAI) |
|          vegin%dleaf,   | real |  |  | chararacteristc legnth of leaf (m) |
|          vegin%ejmax,   | real |  |  | max pot. electron transp rate top leaf(mol/m2/s) |
|          vegin%frac4,   | real |  |  | fraction of c4 plants |
|          vegin%hc,      | real |  |  | roughness height of canopy (veg - snow) |
|          vegin%vlai,    | real |  |  | leaf area index |
|          vegin%xalbnir, | real |  |  |  |
|          vegin%rp20,    | real |  |  | plant respiration coefficient at 20 C |
|          vegin%rpcoef,  | real |  |  | temperature coef nonleaf plant respiration (1/C) |
|          vegin%rs20,    | real |  |  | soil respiration at 20 C [mol m-2 s-1] |
|          vegin%shelrb,  | real |  |  | sheltering factor (dimensionless) |
|          vegin%vegcf,   | real |  |  | kdcorbin, 08/10 |
|          vegin%tminvj,  | real |  |  | min temperature of the start of photosynthesis |
|          vegin%toptvj,  | real |  |  | opt temperature of the start of photosynthesis |
|          vegin%tmaxvj,  | real |  |  | max temperature of the start of photosynthesis |
|          vegin%vbeta,   | real |  |  | |
|          vegin%vcmax,   | real |  |  | max RuBP carboxylation rate top leaf (mol/m2/s) |
|          vegin%xfang,   | real |  |  | leaf angle PARAMETER |
|          vegin%extkn,   | real |  |  | extinction coef for vertical |
|          vegin%vlaimax, | real |  |  | extinction coef for vertical |
|          vegin%wai,     | real |  |  | wood area index (stem+branches+twigs) |
|          vegin%a1gs,    | real |  |  | a1 parameter in stomatal conductance model |
|          vegin%d0gs,    | real |  |  | d0 in stomatal conductance model |
|          vegin%alpha,   | real |  |  | initial slope of J-Q response curve |
|          vegin%convex,  | real |  |  | convexity of J-Q response curve |
|          vegin%cfrd,    | real |  |  | ratio of day respiration to vcmax |
|          vegin%gswmin,  | real |  |  | minimal stomatal conductance |
|          vegin%conkc0,  | real |  |  | Michaelis-menton constant for carboxylase |
|          vegin%conko0,  | real |  |  | Michaelis-menton constant for oxygenase |
|          vegin%ekc,     | real |  |  | activation energy for caroxylagse |
|          vegin%eko,     | real |  |  | acvtivation enegery for oxygenase |
|          vegin%g0,      | real |  |  | Belinda's stomatal model intercept, Ticket #56. |
|          vegin%g1       | real |  |  | Belinda's stomatal model slope, Ticket #56. |
|          vegin%deciduous | logical  |   |   |  flag used for phenology fix |
|          vegin%refl | real  ||||
|          vegin%taul | real  ||||
|          vegin%froot | real  ||| fraction of root in each soil layer |
