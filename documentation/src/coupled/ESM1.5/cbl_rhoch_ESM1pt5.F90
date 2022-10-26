MODULE cbl_rhoch_ESM1pt5_module

  USE cable_def_types_mod, ONLY: mp, nrb
  IMPLICIT NONE

  PUBLIC rhoch_gl, c1_gl, xk_gl

  REAL,ALLOCATABLE, SAVE :: c1_gl(:,:)
  REAL,ALLOCATABLE, SAVE :: rhoch_gl(:,:)
  REAL,ALLOCATABLE, SAVE :: xk_gl(:,:)

END MODULE cbl_rhoch_ESM1pt5_module
