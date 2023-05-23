# pftlookup.csv

This file is required by CASA-CNP and contains biome (plant function types) specific parameters
on soil, vegetation carbon and nutrients dynamics.

!!! Note "Parameter names"

     The Variable column shows the parameter names as they appear in the pftlookup.csv table. The CABLE variable is the name of the variable in the CABLE source code. All alternative parameter names in the Description column are from Appendix A of Wang et al. (2010).

| Variable       | <div style="width:250px">CABLE variable</div> | Description           |
|----------------|-----------------------------------------------|-----------------------|
| Kroot          | `casabiome%kroot(nv)`                         | Unused \( (m^{-1}) \) |
| rootdepth      | `casabiome%rootdepth(nv)`                     | Unused \( (m) \) |
| kuptake        | `casabiome%kuptake(nv)`                       | Unused \( (-) \) |
| krootlen       | `casabiome%krootlen(nv)`                      | Unused \( (m \cdot gC^{-1}) \) |
| KminN          | `casabiome%kminN(nv)`                         | $k_{n,up}$ an empirical parameter relating plant nitrogen uptake to soil mineral N amount. Used in `casa_nuptake()` during computation of N uptake by plants (`xnuptake`) and allocation of uptaken N to plants. \( (gN \cdot m^{-2}) \) |
| Kuplabp        | `casabiome%Kuplabp(nv)`                       | $k_{p,up}$ an empirical parameter relating plant P uptake to labile P pool size in the soil. Used in `casa_puptake()` during calculation of P uptake by plants and the allocation to leaf wood and root. \( (gP \cdot m^{-2}) \)|
| Fracherb       | `xfherbivore(nv)`                             | Fraction of net foliar production consumed by herbivores. Fraction of the turnover rate that affects senescence. `Fracherb=1` then senescence is only related to cold and dry stress and there is no turnover factor, if `Fracherb=0` senescence has an additional turnover rate as set by leaf age. Used to calculate `casabiome%plantrate(nv,leaf)` used in calculation of plant pool senescence rate. I think herbivores were removed when coupled to CABLE. \( (-) \) |
| leaf age       | `leafage(nv)`                                 | $1/\mu_{leaf}$ mean residence time of plant tissue. Used to `calculate casabiome%plantrate(nv,leaf)` used in calculation of plant pool senescence rate. \( (year) \)|
| wood age       | `woodage(nv)`                                 | $1/\mu_{wood}$ mean residence time of plant tissue. Used to `calculate casabiome%plantrate(nv,wood)` used in calculation of plant pool senescence rate. \( (year) \)|
| froot age      | `frootage(nv)`                                | $1/\mu_{wood}$ mean residence time of plant tissue. Used to calculate `casabiome%plantrate(nv,froot)` used in calculation of plant pool senescence rate. \( (year) \)  |
| met age        | `metage(nv)`                                  | $1/\mu_{metabolic}$ mean residence time of metabolic litter. Used to calculate `casabiome%litterrate(nv,metb)` used in calculation of decomposition rate of litter pool. \( (year) \) |
| str age        | `strage(nv)`                                  | $1/\mu_{structural}$ mean residence time of structural litter. Used to calculate `casabiome%litterrate(nv,str)` used in calculation of decomposition rate of litter pool. \( (year) \)|
| cwd age        | `cwdage(nv)`                                  | $1/\mu_{cwd}$ mean residence time of coarse woody debris litter. Used to calculate `casabiome%litterrate(nv,cwd)` used in calculation of decomposition rate of litter pool. \( (year) \) |
| mic age        | `micage(nv)`                                  | $1/\mu_{microbial}$ mean residence time of microbial biomass in soil. Used to calculate `casabiome%soilrate(nv,mic)` used in calculation of decomposition rate of soil pool. \( (year) \)|
| slow age       | `slowage(nv)`                                 | $1/\mu_{slow}$ mean residence time of slow pool in soil. Used to calculate `casabiome%soilrate(nv,slow)` used in calculation of decomposition rate of soil pool. \( (year) \)|
| pass age       | `passage(nv)`                                 | $1/\mu_{passive}$ mean residence time of passive pool in soil. Used to calculate `casabiome%soilrate(nv,pass)` used in calculation of decomposition rate of soil pool. \( (year) \)|
| klabile        | `clabileage(nv)`                              | $1/\mu_{labile}$ mean residence time of labile pool. Used to calculate `kclabrate(nv)`, used to calculate maintenance respiration of woody and fine roots in `casa_rplant()`. \( (year) \)  |
| SLA            | `slawright(nv)` & `casabiome%sla(nv)`         | Specific leaf area. Used to calculate `glai` during conversion of leaf carbon pool to leaf area index. \( (m^2 \cdot gC^{-1}) \)|
| Calloc_leaf    | `casabiome%fracnpptoP(nv,leaf)`               | $a_{leaf}$ NPP allocation coefficient during steady leaf growth. Used in `casa_allocation()`. \( (-) \) |
| Calloc_wood    | `casabiome%fracnpptoP(nv,wood)`               | $a_{wood}$ NPP allocation coefficient during steady leaf growth. Used in `casa_allocation()`. \( (-) \) |
| Calloc_froot   | `casabiome%fracnpptoP(nv,froot)`              | $a_{root}$ NPP allocation coefficient during steady leaf growth. Used in `casa_allocation()`. \( (-) \) |
| rmleaf         | `casabiome%rmplant(nv,leaf)`                  | Tissue respiration rate at 10°C. Used to calculate maintenance respiration of leaf. \( (year^{-1}) \)|
| rmwood         | `casabiome%rmplant(nv,wood)`                  | Tissue respiration rate at 10°C. Used to calculate maintenance respiration of wood. \( (year^{-1}) \)|
| rmfroot        | `casabiome%rmplant(nv,froot)`                 | Tissue respiration rate at 10°C. Used to calculate maintenance respiration of fine roots. \( (year^{-1}) \)|
| rmclabile      | Not loaded                                    | Unused \( (year^{-1}) \) |
| C:N leaf       | `ratioCNplant(nv,leaf)`                       | $n_{leaf}$ Carbon:Nitrogen ratio of leaf biomass. Used to calculate `casabiome%ratioPcplantmin/max(nv,leaf)` used in `casa_Prequire()`. \( (gC \cdot gN^{-1}) \) |
| C:N wood       | `ratioCNplant(nv,wood)`                       | $n_{wood}$ Carbon:Nitrogen ratio of wood biomass. Used to calculate `casabiome%ratioPcplantmin/max(nv,wood)` used in `casa_Prequire()`. \( (gC \cdot gN^{-1}) \) |
| C:N froot      | `ratioCNplant(nv,froot)`                      | $n_{root}$ Carbon:Nitrogen ratio of fine root biomass. Used to calculate `casabiome%ratioPcplantmin/max(nv,froot)` used in `casa_Prequire()`. \( (gC \cdot gN^{-1}) \) |
| Ntrans_leaf    | `casabiome%ftransNPtoL(nv,leaf)`              | Fraction of nutrient fluxes to litter pools. \( (-) \)  |
| Ntrans_wood    | `casabiome%ftransNPtoL(nv,wood)`              | Fraction of nutrient fluxes to litter pools. \( (-) \)  |
| Ntrans_frt     | `casabiome%ftransNPtoL(nv,froot)`             | Fraction of nutrient fluxes to litter pools. \( (-) \)  |
| lignin leaf    | `casabiome%fracligninplant(nv,leaf)`          | Fraction of lignin in plant leaf. \( (g lignin \cdot gC^{-1}) \) |
| lignin CWD     | `casabiome%fracligninplant(nv,wood)`          | Fraction of lignin in plant wood. \( (g lignin \cdot gC^{-1}) \) |
| lignin froot   | `casabiome%fracligninplant(nv,froot)`         | Fraction of lignin in plant fine roots. \( (g lignin \cdot gC^{-1}) \) |
| C:N mic        | `ratioCNsoil(nv,mic)`                         | $n_{mic}$ Carbon:Nitrogen ratio of microbial biomass. Used to calculate `casapool%ratioNCsoil(npt,:)`. \( (gC \cdot gN^{-1}) \) |
| C:N slow       | `ratioCNsoil(nv,slow)`                        | $n_{slow}$ Carbon:Nitrogen ratio of slow soil pool. Used to calculate `casapool%ratioNCsoil(npt,:)`. \( (gC \cdot gN^{-1}) \) |
| C:N pass       | `ratioCNsoil(nv,pass)`                        | $n_{pass}$ Carbon:Nitrogen ratio of passive soil pool. Used to calculate `casapool%ratioNCsoil(npt,:)`. \( (gC \cdot gN^{-1}) \) |
| C:Nmicmin      | `ratioCNsoilmin(nv,mic)`                      | Minimum Carbon:Nitrogen ratio in microbial pool. \( (gC \cdot gN^{-1}) \) |
| C:Nslowmin     | `ratioCNsoilmin(nv,slow)`                     | Minimum Carbon:Nitrogen ratio in slow pool. \( (gC \cdot gN^{-1}) \) |
| C:Npassmin     | `ratioCNsoilmin(nv,pass)`                     | Minimum Carbon:Nitrogen ratio in passive pool. \( (gC \cdot gN^{-1}) \) |
| C:Nmaxmic      | `ratioCNsoilmax(nv,mic)`                      | Maximum Carbon:Nitrogen ratio in microbial pool. \( (gC \cdot gN^{-1}) \) |
| C:Nslowmax     | `ratioCNsoilmax(nv,slow)`                     | Maximum Carbon:Nitrogen ratio in slow pool. \( (gC \cdot gN^{-1}) \) |
| C:Npassmax     | `ratioCNsoilmax(nv,pass)`                     | Maximum Carbon:Nitrogen ratio in passive pool. \( (gC \cdot gN^{-1}) \) |
| Laimax         | `casabiome%glaimax(nv)`                       | Maximum leaf area index. \( (-) \) |
| Laimin         | `casabiome%glaimin(nv)`                       | Minimum leaf area index. \( (-) \) |
| Leaf C         | `cleaf(nv)`                                   | Unused \( (gC \cdot m^{-2}) \) |
| Wood C         | `cwood(nv)`                                   | Unused \( (gC \cdot m^{-2}) \) |
| Froot C        | `cfroot(nv)`                                  | Unused \( (gC \cdot m^{-2}) \) |
| met C          | `cmet(nv)`                                    | Unused \( (gC \cdot m^{-2}) \) |
| str C          | `cstr(nv)`                                    | Unused \( (gC \cdot m^{-2}) \) |
| CWD C          | `ccwd(nv)`                                    | Unused \( (gC \cdot m^{-2}) \) |
| mic C          | `cmic(nv)`                                    | Unused \( (gC \cdot m^{-2}) \) |
| slow C         | `cslow(nv)`                                   | Unused \( (gC \cdot m^{-2}) \) |
| pass C         | `cpass(nv)`                                   | Unused \( (gC \cdot m^{-2}) \) |
| Tkshed         | `phen%TKshed(nv)`                             | Used in `casa_xrateplant()` for cold and drought stress on death rate of leaf. \( (K) \) |
| xkleafcoldmax  | `xxkleafcoldmax(nv)`                          | Used to calculate `casabiome%xkleafcoldmax(nv)` used in `casa_xrateplant()` for cold and drought stress on death rate of leaf. \( (year^{-1}) \)  |
| xkleafcoldexp  | `casabiome%xkleafcoldexp(nv)`                 | exponent in calculation of cold and drought leaf death rate. \( (-) \) |
| xkleafdrymax   | `xxkleafdrymax(nv)`                           | Used to calculate `casabiome%xkleafdrymax(nv)` used in `casa_xrateplant()` for cold and drought stress on death rate of leaf. \( (year^{-1}) \)|
| xkleafdryexp   | `casabiome%xkleafdryexp(nv)`                  | Exponent in calculation of cold and drought leaf death rate. \( (-) \) |
| Tkchill        | Not loaded                                    | \( (K) \) |
| Tkwarm         | Not loaded                                    | \( (K) \) |
| GDD2stdy       | Not loaded                                    |  |
| nd2onset       | Not loaded                                    | \( (day) \) |
| nd2grow        | Not loaded                                    | \( (day) \) |
| nd2dorm        | Not loaded                                    | \( (day) \) |
| phena          | Not loaded                                    | \( (-) \) |
| phenb          | Not loaded                                    | \( (-) \) |
| phenc          | Not loaded                                    | \( (-) \) |
| N/Cleafmi      | `casabiome%ratioNCplantmin(nv,leaf)`          | Minimum ratio of nitrogen to carbon in leaf pool \( (gN \cdot gC^{-1}) \) |
| N/Cleafmx      | `casabiome%ratioNCplantmax(nv,leaf)`          | Maximum ratio of nitrogen to carbon in leaf pool \( (gN \cdot gC^{-1}) \) |
| N/Cwdmin       | `casabiome%ratioNCplantmin(nv,wood)`          | Minimum ratio of nitrogen to carbon in wood pool \( (gN \cdot gC^{-1}) \) |
| N/Cwdmax       | `casabiome%ratioNCplantmax(nv,wood)`          | Maximum ratio of nitrogen to carbon in wood pool \( (gN \cdot gC^{-1}) \) |
| N/Cfrtmin      | `casabiome%ratioNCplantmin(nv,froot)`         | Minimum ratio of nitrogen to carbon in root pool \( (gN \cdot gC^{-1}) \) |
| N/Cfrtmax      | `casabiome%ratioNCplantmax(nv,froot)`         | Maximum ratio of nitrogen to carbon in root pool \( (gN \cdot gC^{-1}) \) |
| xNminloss      | `casaflux%fNminloss(nv)`                      | Nitrogen mineralisation loss factor \( (-) \) |
| xNleach        | `xfNminleach(nv)`                             | Nitrogen leach loss factor \( (year^{-1}) \) |
| nfixrate       | `xnfixrate(nv)`                               | Nitrogen fixation rate (unused) \( (gN \cdot m^{-2} \cdot year^{-1}) \) |
| Nleaf          | `nleaf(nv)`                                   | Leaf pool nitrogen \( (gN \cdot m^{-2}) \) |
| Nwood          | `nwood(nv)`                                   | Wood pool nitrogen \( (gN \cdot m^{-2}) \) |
| Nfrt           | `nfroot(nv)`                                  | Root pool nitrogen \( (gN \cdot m^{-2}) \) |
| N met          | `nmet(nv)`                                    | Metabolic pool nitrogen \( (gN \cdot m^{-2}) \) |
| N Str          | `nstr(nv)`                                    | Structural pool nitrogen \( (gN \cdot m^{-2}) \) |
| Ncwd           | `ncwd(nv)`                                    | Coarse woody debris pool nitrogen \( (gN \cdot m^{-2}) \) |
| N mic          | `nmic(nv)`                                    | Microbial pool nitrogen \( (gN \cdot m^{-2}) \) |
| N slow         | `nslow(nv)`                                   | Xlow pool nitrogen \( (gN \cdot m^{-2}) \) |
| N Pass         | `npass(nv)`                                   | Passive pool nitrogen \( (gN \cdot m^{-2}) \) |
| N Nmin         | `xnsoilmin(nv)`                               | Minimum soil nitrogen \( (gN \cdot m^{-2}) \) |
| N/Pleafmin     | `xratioNPleafmin(nv)`                         | Minimum nitrogen to phosphorus ratio in leaf pool \( (gN \cdot gP^{-2}) \) |
| N/Pleafmx      | `xratioNPleafmax(nv)`                         | Maximum nitrogen to phosphorus ratio in leaf pool \( (gN \cdot gP^{-2}) \) |
| N/Pwdmin       | `xrationpwoodmin(nv)`                         | Minimum nitrogen to phosphorus ratio in woody pool \( (gN \cdot gP^{-2}) \) |
| N/Pwdmax       | `xrationpwoodmax(nv)`                         | Maximum nitrogen to phosphorus ratio in woody pool \( (gN \cdot gP^{-2}) \) |
| N/Pfrtmin      | `xratioNPfrootmin(nv)`                        | Minimum nitrogen to phosphorus ratio in root pool \( (gN \cdot gP^{-2}) \) |
| N/Pfrtmax      | `xratioNPfrootmax(nv)`                        | Maximum nitrogen to phosphorus ratio in root pool \( (gN \cdot gP^{-2}) \) |
| fpptoL(leaf)   | `casabiome%ftransPPtoL(nv,leaf)`              | Flux factor of leaf phosphorus to litter pools \( (-) \) |
| fpptoL(wd)     | `casabiome%ftransPPtoL(nv,wood)`              | Flux factor of wood phosphorus to litter pools \( (-) \) |
| fpptoL(frt)    | `casabiome%ftransPPtoL(nv,froot)`             | Flux factor of root phosphorus to litter pools \( (-) \) |
| xkmlabp        | `xkmlabp(iso)`                                | Phosphorus absorption \( (gP \cdot m^{-2}) \)  |
| xpsorbmax      | `xpsorbmax(iso)`                              | Maximum phosphorus absorption \( (gP \cdot m^{-2}) \) |
| xfpleach       | `xfPleach(iso)`                               | Phosphorus leaching \( (-) \) |
| N:Psoil (mic)  | `ratioNPsoil(iso,mic)`                        | Nitrogen to phosphorus ratio in microbial soil pool \( (gN \cdot gP^{-2}) \) |
| N:Psoil (slow) | `ratioNPsoil(iso,slow)`                       | Nitrogen to phosphorus ratio in slow soil pool \( (gN \cdot gP^{-2}) \) |
| N:Psoil (pass) | `ratioNPsoil(iso,pass)`                       | Nitrogen to phosphorus ratio in passive soil pool \( (gN \cdot gP^{-2}) \) |
| Pleaf          | `xpleaf(nv)`                                  | Phosphorus in leaf pool \( (gP \cdot m^{-2}) \) |
| Pwood          | `xpwood(nv)`                                  | Phosphorus in wood pool \( (gP \cdot m^{-2}) \) |
| Pfroot         | `xpfroot(nv)`                                 | Phosphorus in root pool \( (gP \cdot m^{-2}) \) |
| Pmet           | `xpmet(nv)`                                   | Phosphorus in metabolic pool \( (gP \cdot m^{-2}) \) |
| Pstr           | `xpstr(nv)`                                   | Phosphorus in structural pool \( (gP \cdot m^{-2}) \) |
| Pcwd           | `xpcwd(nv)`                                   | Phosphorus in coarse woody debris pool \( (gP \cdot m^{-2}) \) |
| Pmic           | `xpmic(nv)`                                   | Phosphorus in microbial pool \( (gP \cdot m^{-2}) \) |
| Pslow          | `xpslow(nv)`                                  | Phosphorus in slow pool \( (gP \cdot m^{-2}) \) |
| Ppass          | `xppass(nv)`                                  | Phosphorus in passive pool \( (gP \cdot m^{-2}) \) |
| Plab           | `xplab(nv)`                                   | Phosphorus in labile pool \( (gP \cdot m^{-2}) \) |
| Psorb          | `xpsorb(nv)`                                  | Phosphorus in sorbed pool \( (gP \cdot m^{-2}) \) |
| Pocc           | `xpocc(nv)`                                   | Phosphorus in occluded pool \( (gP \cdot m^{-2}) \) |

1. Wang, Y. P., Law, R. M. & Pak, B. A global model of carbon, nitrogen and phosphorus cycles for the terrestrial biosphere. Biogeosciences 7, 2261–2282 (2010).

