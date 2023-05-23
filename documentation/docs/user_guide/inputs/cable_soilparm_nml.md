# cable_soilparm.nml options

!!! Note

    The asterisk symbol does not mean arithmetic multiplication in the namelist file.
    It means that the following sequence is repeated multiple times.

!!! Note

    In JULES, unintialised values are assigned as a missing data indicator (MDI)

```fortran
&cable_soilparm

soilin%desc='Coarse sand/Loamy sand',
 'Medium clay loam/silty clay loam/silt loam',
 'Fine clay',
 'Coarse-medium sandy loam/loam',
 'Coarse-fine sandy clay',
 'Medium-fine silty clay',
 'Coarse-medium-fine sandy clay loam',
 'Organic peat',
 'Permanent ice',
soilin%bch=4.2,7.1,11.4,5.15,10.4,10.4,7.12,5.83,7.1,
soilin%clay=0.09,0.3,0.67,0.2,0.42,0.48,0.27,0.17,0.3,
soilin%css=7*850,1920,2100,
soilin%hyds=0.000166,0.000004,0.000001,0.000021,0.000002,0.000001,0.000006,0.0008,0.000001,
soilin%rhosoil=1600,1600,1381,1373,1476,1521,1373,1537,910,
soilin%sand=0.83,0.37,0.16,0.6,0.52,0.27,0.58,0.13,0.37,
soilin%sfc=0.143,0.301,0.367,0.218,0.31,0.37,0.255,0.45,0.301,
soilin%silt=0.08,0.33,0.17,0.2,0.06,0.25,0.15,0.7,0.33,
soilin%ssat=0.398,0.479,0.482,0.443,0.426,0.482,0.42,0.451,0.479,
soilin%sucs=-0.106,-0.591,-0.405,-0.348,-0.153,-0.49,-0.299,-0.356,-0.153,
soilin%swilt=0.072,0.216,0.286,0.135,0.219,0.283,0.175,0.395,0.216,
/
```

| Namelist variable | Type              | Available values | Default values | Description                                                                              |
|-------------------|-------------------|------------------|----------------|------------------------------------------------------------------------------------------|
| soilin%desc       | character(*)      |                  | uninitialised  | Description of soil types                                                                |
| soilin%bch        | real(n_soiltypes) | >=0.0            | uninitialised  | Parameter b in Campbell equation for the pore size distribution index (1,2) \( ( - ) \). |
| soilin%silt       | real(n_soiltypes) | >=0.0            | uninitialised  | Fraction of soil which is silt \( ( - ) \).                                              |
| soilin%clay       | real(n_soiltypes) | >=0.0            | uninitialised  | Fraction of soil which is clay \( ( - ) \).                                              |
| soilin%sand       | real(n_soiltypes) | >=0.0            | uninitialised  | Fraction of soil which is sand \( ( - ) \).                                              |
| soilin%swilt      | real(n_soiltypes) | >=0.0            | uninitialised  | Volume of H2O at wilting \( (m^{3} \cdot m^{-3}) \).                                     |
| soilin%sfc        | real(n_soiltypes) | >=0.0            | uninitialised  | Volume of H2O at field capacity \( (m^{3} \cdot m^{-3}) \).                              |
| soilin%ssat       | real(n_soiltypes) | >=0.0            | uninitialised  | Volume of H2O at saturation \( (m^{3} \cdot m^{-3}) \).                                  |
| soiln%hyds        | real(n_soiltypes) | >=0.0            | uninitialised  | Hydraulic conductivity at saturation \( (m^{-1}) \).                                     |
| soilin%sucs       | real(n_soiltypes) | <=0.0            | uninitialised  | Suction at saturation \( (m) \).                                                         |
| soilin%rhosoil    | real(n_soiltypes) | >=0.0            | uninitialised  | Soil bulk density \( (kg \cdot m^{-3}) \).                                               |
| soilin%css        | real(n_soiltypes) | >=0.0            | uninitialised  | Soil specific heat capacity \( (J \cdot kg^{-1} \cdot K^{-1}) \).                        |


1) Brooks, R. H., and Corey, A. T. 1964. Hydraulic Properties of Porous Media. Fort Collins: Colorado State University.

2) Campbell, G. S. 1974. A simple method for determining unsaturated conductivity from moisture retention data. Soil Science, 117(6), 311-314. [doi:10.1097/00010694-197406000-00001
