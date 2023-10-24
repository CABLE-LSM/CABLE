
## Introduction

CABLE is a community software project. 
For the benefit of the whole community, developers are required to adhere to certain coding standards that better enable portability, efficient maintenance, future development, performance and readability. 

Submitted developments that are inconsistent with the standards discussed below, and/or poorly commented/documented code will not be accepted without revision. 

The JULES and CABLE (JAC) project is a long term initiative to include CABLE as an inter-operable Land-Surface Model (LSM) within the Unified Model (UM). For several reasons, CABLE is coupled to the UM via its native LSM, JULES. From a coding standards perspective, both the UM and JULES developments are strictly goverened by the UK Met Office and thus CABLE's coding standards are somewhat constrained the [JULES coding standards](http://jules-lsm.github.io/coding_standards/).

Currently not all of the CABLE code is compliant with the standards described here. Work is in progress to apply these standards across the existing code base, however these standards are expected to be met by any new development.

Whilst syntax/text formatting style is relatively trivial, it must still be adhered to. Some more significant general rules follow.

## General Rules

### 1. IO must be restricted to the top level

In terms of the offline application of the model this means IO should be restricted to the `offline/` directory. There should be no IO code elsewhere in the CABLE directory. 

### 2. Comment liberally

Detrimental to any software is the lack of meaningful and/or outdated comments. 
Liberal inline commenting is expected, as well as use of Markdown for FORD generated documentation. See [cable_roughness.F90](https://github.com/CABLE-LSM/CABLE/blob/main/src/science/roughness/cable_roughness.F90) for an example, which is then rendered [here](https://cable.readthedocs.io/en/latest/api/module/cable_roughness_module.html)

### 3. Do not use WHERE loops

`WHERE` loops are sometimes used instead of `DO` loops with an `IF`. Although they might be easier to write and may improve readability, `WHERE` loops are difficult to debug and hamper portability. They are therefore discouraged. 

### 4. Variable names

Developers are encouraged to name variables meaningfully, and remember the number of characters for variables is not anymore limited. A variable name should help to understand what is going on when reading the code without fully relying on comment lines. 

Maintain globally consistent names of variables across files, subroutines and modules. 

### 5. File names

All files in CABLE source code should use the prefix `cbl_` in their name.
For example: `cbl_egName.F90`. 

### 6. MODULE names

Everything should be contained within modules.
Modules should be named using the filename followed by a `_mod` suffix. For example, the file `cbl_egName.F90` should contain the module declaration `MODULE cbl_egName_mod`

### 7. SUBROUTINE names

Subroutines should be of the form "\*_cbl()". Following from the above file `egName_cbl.F90` `SUBROUTINE egName_cbl()`. This is not only helpful for consistency but necessary within JULES/UM applications where it distinguishes the CABLE version of an analogous JULES subroutine. 

### 8. Miscellaneous

 * It will not be permitted to `USE` data through `MODULE`s. Data will be be passed through argument lists. 
 * The last character of a line **must** not be placed beyond column 80. 
 * `SAVE` attributes are not permitted.
 * Intrinsic Fortran subroutines, functions, etc should be capitalized
 * Subroutines start on the margin. This includes MODULE, CONTAINS, SUBROUTINE, USE statements, and declarations. Indentation of two spaces occurs in loops/layers. This is vital in readily identifying the extent of loops/layers. i.e. Under which condition(s) a section of code is relevantin an IF layer, or at what index is being traversed in a DO loop. 

 * Specify `IMPLICIT NONE` in all modules to prevent type mismatch errors 
 * Dummy arguments must be declared with the `INTENT()` attribute
 * Declare in the `MODULE` preamble which elements are `PUBLIC`
 * Import specific components using `ONLY` with `USE` statements to avoid inadvertent inheritance and
   increase transparency of underlying data flow 
 * Label each side of a loop, especially if it extends
   more than a page.
 * All `IF` constructs must be of the form:
```fortran
IF() THEN 
  Even if they only contain a simple, single line statement following the IF. 
END IF
```
 * In summary, be pedantic, be clear

## Final Note about JAC

Scientific models have naturally evolved over time to include increasing complexity. Coupled with the demand for higher resolution, larger ensembles etc, this has required ever increasing computational performance.
For several decades, this icreasing performance demand has been met by technological advances in hardware.   
This era has come to an end.
The Next Generation Modelling System (NGMS) is a major program at the UKMO to address this issue.
The NGMS is an initiative to optimize the software of their entire suite of models.
The specific coding standards required by the NGMS have not yet been strictly defined. 

In general, models (including CABLE) have evolved over decades, involving dozens of developers. A key feature emerging from NGMS is the tight reign on memory. As listed above, it will not be permitted to `USE` data through `MODULE`s. Data will be be passed through argument lists. There are further rules emerging in this respect which relax the rules slightly to allow time independent scalars to be `USE`d. Nevertheless, it is straightforward enough to pass these as well. 

## Specific Develoments
Code developments broadly fall into two categories: bug fixes and new developments. Both can involve a wide range of complexity. However, at the extremes a bug fix may be as simple as moving a bracket. Alternatively, the bug fix might require you to rewrite a section, or even several sections of code. 
A new development might be as simple as an alternative section or calculation, or it might be an entirely new model that can be plugged in to CABLE, either as an alternative model, or an extra feature to be used in the standard model.

In all cases, an issue should be raised. The accompanying fix contained should then be implemented in a unique branch, corresponding to the raised issue. A bugfix might simply involve reference to this issue. For a more complicated, added feature, it is a requirement that this addition can be isolated and switched off via a configuration switch, so that the model can be also run as if there had been no alteration. There is of course every possibility that this feature will be accepted as a standard part of the model. Although this is not the expectation initially and the developer should be mindful of this fundamental requirement. 

We proceed assuming that the "issue" raised requires adding a new and significant feature to CABLE, requiring the addition of new files/modules to CABLE. Lesser modifications are covered under this assumption anyway.

For discussion purposes we assume that you are including a new model to do some new stuff. This can conveniently be added in a single file. 

## Code Template 

Assuming the *new* model is called from eg_driver.F90.

```fortran
PROGRAM eg_driver()

USE new_mod_cbl, ONLY : new

IMPLICIT NONE
Header Declarations
...
CALL eg_subr( Arguments )
...

END PROGRAM eg_driver
```
Where,

```fortran
MODULE cbl_eg_mod
! Module description in markdown - e.g. [cable_roughness.F90](https://github.com/CABLE-LSM/CABLE/blob/main/src/science/roughness/cable_roughness.F90) for an example, which is then rendered [here](https://cable.readthedocs.io/en/latest/api/module/cable_roughness_module.html)

PUBLIC :: eg_subr 
PRIVATE

CONTAINS

SUBROUTINE eg_subr ( Arguments )
...
END SUBROUTINE eg_subr 

END MODULE cbl_eg_mod
```

Where,

```fortran
SUBROUTINE eg_subr( eg_arg, ilen, jlen, sometype )
! Subroutine description in markdown - see??

USE cbl_some_type_mod, ONLY : some_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: ilen              ! example arg length
INTEGER, INTENT(IN) :: jlen              ! example arg length
REAL, INTENT(OUT) :: eg_arg(ilen,jlen)   ! example arg [kg//m^2/s]

!local variables
REAL :: eg_local(ilen,jlen)              ! example local [kg//m^2/s]

ilen_loop: DO i=1,ilen

  jlen_loop: DO j=1,jlen
  
    condition1: IF( condition1 ) THEN
  
      condition1.1: IF( condition1.1 ) THEN  
        
        eg_arg = some_type%some_type_member1 * some_type%some_type_member2     &
                * some_type%some_type_member3 
      
      END IF condition1.1
      
    ELSE
        
      condition1.2: IF( condition1.2 ) THEN 
      
        eg_arg = some_type%some_type_member1 * some_type%some_type_member2     &
                * some_type%some_type_member4 
      
      END IF condition1.2
                          
    END IF condition1
  
  END DO jlen_loop

END DO ilen_loop

RETURN
END SUBROUTINE eg_subr 
```
