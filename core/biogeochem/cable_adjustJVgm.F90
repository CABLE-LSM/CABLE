!
! ==============================================================================
! Purpose: adjustment of Ci-based Jmax and Vcmax to their Cc-based values
!          (accounting for a finite mesophyll conductance) using a nonlinear
!          curve fitting routine as described in Knauer et al. 2019 GCB.
!
! Called from: SUBROUTINE bgcdriver in casa_cable.F90
!
! History: Juergen Knauer July/August 2019
! ==============================================================================
!
MODULE cable_adjust_JV_gm_module

  USE cable_def_types_mod, ONLY: dp => r_2
  USE cable_def_types_mod, ONLY: mp, veg_parameter_type
  USE cable_data_module,   ONLY: icanopy_type, point2constants
  USE minpack

  TYPE(icanopy_type) :: C

  INTEGER, PARAMETER :: nrci=3000
  REAL(dp) :: gmmax25, Vcmax25Ci, Jmax25Ci, Vcmax25Cc, Jmax25Cc
  REAL(dp) :: Rd
  REAL(dp) :: Kc_ci, Ko_ci, gammastar_ci, Km_ci
  REAL(dp) :: Kc_cc, Ko_cc, gammastar_cc, Km_cc
  
CONTAINS 

  SUBROUTINE adjust_JV_gm(veg)  

    IMPLICIT NONE

    TYPE (veg_parameter_type), INTENT(INOUT)  :: veg  ! vegetation parameters

    ! local variables
    LOGICAL  :: Cc_based_OK, sw ! sw = stability switch
    INTEGER  :: p,i,k,z
    INTEGER  :: kmax=20  ! maximum nr of iterations (inner loop)
    INTEGER  :: zmax=8   ! maximum nr of iterations (outer loop)
    INTEGER  :: lAn, cntr
    REAL(dp) :: vstart, v
    REAL(dp) :: Vcmax25Cct1  ! Vcmax25Cc of previous iteration
    REAL(dp) :: Vcmax_diff
    REAL(dp) :: maxdiff=0.002e-6_dp
    REAL(dp), DIMENSION(nrci)  :: An1, Ci1
    REAL(dp), DIMENSION(:), ALLOCATABLE :: An, Ci, Cc, An_Cc
    
    ! MINPACK params
    INTEGER, PARAMETER      :: N=2     ! Number of variables
    REAL(dp), DIMENSION(N)  :: X
    REAL(dp), ALLOCATABLE   :: FVEC(:)
    INTEGER                 :: info
    REAL(dp)                :: tol=0.00001_dp

    
    ! assign local ptrs to constants defined in cable_data_module
    CALL point2constants(C)
    
    DO p=1,mp

      Ci1          = (/(real(i,dp),i=1,nrci,1)/) / 2.0_dp * 1.0e-6_dp  ! 1-1500 umol mol-1
      Rd           = real(veg%cfrd(p) * veg%vcmax(p),dp)
      gmmax25      = real(veg%gmmax(p),dp)
      Vcmax25Ci    = real(veg%vcmax(p),dp)
      Jmax25Ci     = real(veg%ejmax(p),dp)
      Kc_ci        = real(C%conkc0,dp)
      Ko_ci        = real(C%conko0,dp)
      gammastar_ci = real(C%gam0,dp)
      Kc_cc        = real(C%conkc0cc,dp)
      Ko_cc        = real(C%conko0cc,dp)
      gammastar_cc = real(C%gam0cc,dp)

      Km_ci        = Kc_ci * (1.0_dp + 0.21_dp / Ko_ci)
      Km_cc        = Kc_cc * (1.0_dp + 0.21_dp / Ko_cc)

      IF (veg%frac4(p) .lt. 0.001_dp) THEN ! not C4

        !! 1) Calculate An-Ci curve
        CALL PHOTOSYN25(Ci1,nrci,Vcmax25Ci,Jmax25Ci,Rd,Km_ci,gammastar_ci,An1)

        !! 2) Exclude negative parts of the An-Ci curve
        lAn = count(An1 .GT. 0.0_dp)

        ALLOCATE(An(lAn))
        ALLOCATE(An_Cc(lAn))
        ALLOCATE(Ci(lAn))
        ALLOCATE(Cc(lAn))
        ALLOCATE(FVEC(lAn))

        cntr = 0
        DO i = 1, 3000
          IF (An1(i) .GT. 0.0_dp) THEN
            cntr = cntr + 1
            An(cntr) = An1(i)
            Ci(cntr) = Ci1(i)
          END IF
        END DO

        Cc_based_OK = .FALSE.
        z = 0
        !! 3) calculate Cc based on gm and An
        DO WHILE(.NOT. Cc_based_OK .AND. z < zmax) ! if it iterates more than once, check gm and Vcmax, Jmax

           z = z + 1
           k = 0
           X(:) = [Vcmax25Ci,Jmax25Ci]
           sw = .FALSE.
           vstart = 1.0_dp
           Vcmax25Cct1 = Vcmax25Ci
           Vcmax_diff = 1.0e-6_dp
           An_Cc = An

           DO WHILE (Vcmax_diff > maxdiff .AND. k < kmax)

              k = k + 1
              Cc = Ci - An_Cc / gmmax25

              CALL LMDIF1(PHOTOSYN25_f,lAn,N,X,FVEC,tol,info,An_Cc,Cc,Rd,Km_cc,gammastar_cc)
              Vcmax25Cc = X(1)
              Jmax25Cc  = X(2)

              Vcmax_diff = ABS(Vcmax25Cc - Vcmax25Cct1)
              Vcmax25Cct1 = Vcmax25Cc

              CALL PHOTOSYN25(Cc,lAn,Vcmax25Cc,Jmax25Cc,Rd,Km_cc,gammastar_cc,An_Cc)

              ! safety switch ensuring stability
              IF (MINVAL(An_Cc) < 0.0_dp .AND. (.NOT. sw)) THEN
                 sw = .TRUE.
                 v  = vstart
              ENDIF

              IF (sw) THEN
                 v = MAX(v - (vstart/(0.8_dp*kmax)),0.0_dp)
                 An_Cc = v * An + (1.0_dp-v) * An_Cc
              ENDIF
           END DO   
           !! Avoid unrealistic Vcmax and Jmax values
           IF (Vcmax25Cc < 0.9_dp*Vcmax25Ci .OR. Vcmax25Cc > 2.5_dp*Vcmax25Ci &
               .OR. Jmax25Cc < 0.9_dp*Jmax25Ci .OR. Jmax25Cc > 1.5_dp*Jmax25Ci) THEN
              gmmax25 = 1.2_dp * gmmax25  ! If no solution, try again with higher gmmax25          
           ELSE       
              Cc_based_OK = .TRUE.
              veg%vcmaxcc(p) = real(Vcmax25Cc)
              veg%ejmaxcc(p) = real(Jmax25Cc)
           ENDIF

        END DO
     
        DEALLOCATE(An)
        DEALLOCATE(An_Cc)
        DEALLOCATE(Ci)
        DEALLOCATE(Cc)
        DEALLOCATE(FVEC)
        
       
      ELSE  ! For C4 plants same as for Ci for now...

        veg%vcmaxcc(p) = real(Vcmax25Ci)
        veg%ejmaxcc(p) = real(Jmax25Ci)
         
      ENDIF ! C4 flag
    END DO ! tile loop
          
  END SUBROUTINE adjust_JV_gm


  ! Function to use within LMDIF1
  SUBROUTINE PHOTOSYN25_f(M,N,X,FVEC,IFLAG,Anx,Cix,Rd,Km,gammastar)

    INTEGER,                INTENT(IN)    :: M,N,IFLAG
    REAL(dp), DIMENSION(N), INTENT(INOUT) :: X
    REAL(dp), DIMENSION(M), INTENT(OUT)   :: FVEC
    REAL(dp), DIMENSION(M), INTENT(IN)    :: Anx,Cix
    REAL(dp),               INTENT(IN)    :: Rd,Km,gammastar
    ! local
    REAL(dp), DIMENSION(M) :: Ac, Aj

    
    Ac = (X(1) * (Cix - gammastar) / (Cix + Km))
    Aj = (X(2) * (Cix - gammastar) / 4.0 / (Cix + 2.0 * gammastar))
    
    ! avoid discontinuity (e.g. Duursma 2015, PLOS ONE)
    FVEC = Anx  - ( (Ac + Aj - SQRT((Ac + Aj)**2 - 4.0*0.99999999_dp*Ac*Aj)) / &
                    (2.0*0.99999999_dp) - Rd  &
                  )
  
  END SUBROUTINE PHOTOSYN25_f


  
  ! Function to calculate An-Ci curve under standard conditions
  SUBROUTINE PHOTOSYN25(Ciz,nrci,Vcmax25,Jmax25,Rd,Km,gammastar,Anz)

    INTEGER,                   INTENT(IN)  :: nrci
    REAL(dp), DIMENSION(nrci), INTENT(IN)  :: Ciz
    REAL(dp),                  INTENT(IN)  :: Vcmax25,Jmax25,Rd,Km,gammastar
    REAL(dp), DIMENSION(nrci), INTENT(OUT) :: Anz
    ! local
    REAL(dp), DIMENSION(nrci)  :: Wc,We


    ! Rubisco-limited
    Wc =  Vcmax25 * (Ciz - gammastar) / (Ciz + Km)

    ! RuBP regeneration-limited
    We =  Jmax25 * (Ciz - gammastar) / 4.0 / (Ciz + 2.0 * gammastar)

    ! Net photosynthesis
    Anz = Min(Wc,We) - Rd

  END SUBROUTINE PHOTOSYN25

  
END MODULE cable_adjust_JV_gm_module
