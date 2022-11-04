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

  use cable_def_types_mod, only: dp => r_2
  use cable_def_types_mod, only: mp, veg_parameter_type
  use cable_data_module,   only: icanopy_type, point2constants
  use cable_abort_module,  only: nc_abort
  use cable_canopy_module, only: light_inhibition
  use netcdf
  use minpack

  type(icanopy_type) :: C

  integer, parameter :: nrci=3000
  integer, parameter :: nrcic4=1200
  real(dp) :: gmmax25, Vcmax25Ci, Jmax25Ci, Vcmax25Cc, Jmax25Cc, k25Ci, k25Cc
  real(dp) :: Rd
  real(dp) :: Kc_ci, Ko_ci, gammastar_ci, Km_ci
  real(dp) :: Kc_cc, Ko_cc, gammastar_cc, Km_cc

  ! gm LUT
  real(dp), dimension(:,:,:,:), allocatable :: LUT_VcmaxJmax ! Lookup table with Cc-based Vcmax and Jmax
  real(dp), dimension(:),       allocatable :: LUT_gm        ! gm values in gm LUT
  real(dp), dimension(:),       allocatable :: LUT_Vcmax     ! Vcmax_ci values in gm LUT
  real(dp), dimension(:),       allocatable :: LUT_Rd        ! Rd values in gm LUT

CONTAINS

  SUBROUTINE adjust_JV_gm(veg, p)

    IMPLICIT NONE

    type(veg_parameter_type), intent(inout) :: veg ! vegetation parameters
    integer,                  intent(in)    :: p   ! vegetation type

    ! local variables
    logical  :: Cc_based_OK, sw ! sw = stability switch
    integer  :: i, k, z
    integer  :: kmax = 20  ! maximum nr of iterations (inner loop)
    integer  :: zmax = 8   ! maximum nr of iterations (outer loop)
    integer  :: lAn
    real(dp) :: vstart, v
    real(dp) :: Vcmax25Cct1  ! Vcmax25Cc of previous iteration
    real(dp) :: Vcmax_diff
    real(dp) :: maxdiff = 0.002e-6_dp
    real(dp), dimension(nrci)  :: An1, Ci1
    real(dp), dimension(:), allocatable :: An, Ci, Cc, An_Cc

    ! MINPACK params
    integer, parameter      :: N = 2 ! Number of variables
    real(dp), dimension(N)  :: X
    real(dp), allocatable   :: fvec(:)
    integer                 :: info
    real(dp)                :: tol = 0.00001_dp

    ! assign local ptrs to constants defined in cable_data_module
    call point2constants(C)

    Ci1          = (/(real(i,dp), i=1, nrci, 1)/) / 2.0_dp * 1.0e-6_dp ! 1-1500 umol mol-1
    Rd           = real(veg%cfrd(p) * veg%vcmax(p) * light_inhibition(1200.0), dp)
    gmmax25      = real(veg%gm(p), dp)
    Vcmax25Ci    = real(veg%vcmax(p), dp)
    Jmax25Ci     = real(veg%ejmax(p), dp)
    Kc_ci        = real(C%conkc0, dp)
    Ko_ci        = real(C%conko0, dp)
    gammastar_ci = real(C%gam0, dp)
    Kc_cc        = real(C%conkc0cc, dp)
    Ko_cc        = real(C%conko0cc, dp)
    gammastar_cc = real(C%gam0cc, dp)

    Km_ci = Kc_ci * (1.0_dp + 0.21_dp / Ko_ci)
    Km_cc = Kc_cc * (1.0_dp + 0.21_dp / Ko_cc)

    if (veg%frac4(p) .lt. 0.001) then ! not C4

      !! 1) Calculate An-Ci curve
      call photosyn25(Ci1, nrci, Vcmax25Ci, Jmax25Ci, Rd, Km_ci, gammastar_ci, An1)

      !! 2) Exclude negative parts of the An-Ci curve
      lAn = count(An1 > 0.0_dp)

      allocate(An(lAn))
      allocate(An_Cc(lAn))
      allocate(Ci(lAn))
      allocate(Cc(lAn))
      allocate(fvec(lAn))

      An = pack(An1, An1 > 0.0_dp)
      Ci = pack(Ci1, An1 > 0.0_dp)

      Cc_based_OK = .false.
      z = 0
      ! 3) calculate Cc based on gm and An
      do while (.not. Cc_based_OK .and. z < zmax) ! if it iterates more than once, check gm and Vcmax, Jmax

         z = z + 1
         k = 0
         X(:) = [Vcmax25Ci,Jmax25Ci]
         sw = .false.
         vstart = 1.0_dp
         Vcmax25Cct1 = Vcmax25Ci
         Vcmax_diff = 1.0e-6_dp
         An_Cc = An

         do while (Vcmax_diff > maxdiff .AND. k < kmax)
            k  = k + 1
            Cc = Ci - An_Cc / gmmax25

            call lmdif1(photosyn25_f, lAn, N, X, fvec, tol, info, An_Cc, Cc, Rd, Km_cc, gammastar_cc)
            Vcmax25Cc = X(1)
            Jmax25Cc  = X(2)

            Vcmax_diff  = abs(Vcmax25Cc - Vcmax25Cct1)
            Vcmax25Cct1 = Vcmax25Cc

            call photosyn25(Cc, lAn, Vcmax25Cc, Jmax25Cc, Rd, Km_cc, gammastar_cc, An_Cc)

            ! safety switch ensuring stability
            if (minval(An_Cc) < 0.0_dp .and. (.not. sw)) then
               sw = .true.
               v  = vstart
            endif

            if (sw) then
               v = max(v - (vstart/(0.8_dp*kmax)),0.0_dp)
               An_Cc = v * An + (1.0_dp-v) * An_Cc
            endif
         end do

         !! Avoid unrealistic Vcmax and Jmax values
         if ((Vcmax25Cc < 0.9_dp*Vcmax25Ci) .or. (Vcmax25Cc > 2.5_dp*Vcmax25Ci) &
             .or. (Jmax25Cc < 0.9_dp*Jmax25Ci) .or. (jmax25cc > 1.5_dp*jmax25ci)) then
            gmmax25 = 1.2_dp * gmmax25 ! If no solution, try again with higher gmmax25
         else
            Cc_based_OK = .true.
            veg%vcmaxcc(p) = real(Vcmax25Cc)
            veg%ejmaxcc(p) = real(Jmax25Cc)
         endif

      end do

      deallocate(An)
      deallocate(An_Cc)
      deallocate(Ci)
      deallocate(Cc)
      deallocate(fvec)

    ELSE ! C4 (Vcmax and Jmax do not change with gm in C4 plants)

      veg%vcmaxcc(p) = real(Vcmax25Ci)
      veg%ejmaxcc(p) = real(Jmax25Ci)

    endif ! c4 flag

  END SUBROUTINE adjust_JV_gm


  ! Function to use within LMDIF1
  subroutine photosyn25_f(M, N, X, fvec, iflag, Anx, Cix, Rd, Km, gammastar)

    integer,                intent(in)    :: M, N, iflag
    real(dp), dimension(N), intent(inout) :: X
    real(dp), dimension(M), intent(out)   :: fvec
    real(dp), dimension(M), intent(in)    :: Anx, Cix
    real(dp),               intent(in)    :: Rd, Km, gammastar
    ! local
    real(dp), dimension(M) :: Ac, Aj

    Ac = (X(1) * (Cix - gammastar) / (Cix + Km))
    Aj = (X(2) * (Cix - gammastar) / 4.0 / (Cix + 2.0 * gammastar))

    ! avoid discontinuity (e.g. Duursma 2015, PLOS ONE)
    fvec = Anx  - ( (Ac + Aj - SQRT((Ac + Aj)**2 - 4.0*0.99999999_dp*Ac*Aj)) / &
                    (2.0*0.99999999_dp) - Rd )

  END SUBROUTINE photosyn25_f


  ! Function to calculate An-Ci curve under standard conditions
  subroutine photosyn25(Ciz, nrci, Vcmax25, Jmax25, Rd, Km, gammastar, Anz)

    integer,                   intent(in)  :: nrci
    real(dp), dimension(nrci), intent(in)  :: Ciz
    real(dp),                  intent(in)  :: Vcmax25, Jmax25, Rd, Km, gammastar
    real(dp), dimension(nrci), intent(out) :: Anz
    ! local
    real(dp), dimension(nrci)  :: Wc, We

    ! Rubisco-limited
    Wc =  Vcmax25 * (Ciz - gammastar) / (Ciz + Km)

    ! RuBP regeneration-limited
    We =  Jmax25 * (Ciz - gammastar) / 4.0 / (Ciz + 2.0 * gammastar)

    ! Net photosynthesis
    Anz = Min(Wc,We) - Rd

  END SUBROUTINE photosyn25


  SUBROUTINE read_gm_LUT(gm_LUT_file, LUT_VcmaxJmax, LUT_gm, LUT_vcmax, LUT_Rd)
    ! Read lookup table needed for parameter conversion of photosynthetic parameters
    ! from Ci- to Cc-based values (the latter considering gm explicitly)
    use netcdf, only: nf90_open, nf90_nowrite, nf90_noerr, &
         nf90_inq_dimid, nf90_inquire_dimension, &
         nf90_inq_varid, nf90_get_var, &
         nf90_close

    implicit none

    character(len=*),                          intent(in)  :: gm_LUT_file
    real(dp), dimension(:,:,:,:), allocatable, intent(out) :: LUT_VcmaxJmax
    real(dp), dimension(:),       allocatable, intent(out) :: LUT_gm
    real(dp), dimension(:),       allocatable, intent(out) :: LUT_vcmax
    real(dp), dimension(:),       allocatable, intent(out) :: LUT_Rd

    ! local
    integer  :: ncid_gmlut                      ! netcdf ID
    integer  :: ok                              ! netcdf error status
    integer  :: gm_dimid, vcmax_dimid, Rd_dimid ! dimension IDs
    integer  :: gm_len, vcmax_len, Rd_len       ! dimensions of LUT
    integer  :: vcmax_id, jmax_id
    real(dp), dimension(:,:,:), allocatable :: tmp ! for reading

    ok = nf90_open(trim(gm_LUT_file), nf90_nowrite, ncid_gmlut)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error opening gm lookup table.')
    ok = nf90_inq_dimid(ncid_gmlut, 'gm', gm_dimid)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error inquiring dimension gm from LUT.')
    ok = nf90_inq_dimid(ncid_gmlut, 'Vcmax_Ci', vcmax_dimid)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error inquiring dimension Vcmax_Ci from LUT.')
    ok = nf90_inq_dimid(ncid_gmlut, 'Rd', Rd_dimid)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error inquiring dimension Rd from LUT.')

    ok = nf90_inquire_dimension(ncid_gmlut, gm_dimid, len=gm_len)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error inquiring length of dimension gm from LUT.')
    ok = nf90_inquire_dimension(ncid_gmlut, vcmax_dimid, len=vcmax_len)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error inquiring length of dimension Vcmax_Ci from LUT.')
    ok = nf90_inquire_dimension(ncid_gmlut, Rd_dimid, len=Rd_len)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error inquiring length of dimension Rd from LUT.')

    write(*,*) 'gm LUT dimensions:'
    write(*,*) 'gm_len:', gm_len
    write(*,*) 'vcmax_len:', vcmax_len
    write(*,*) 'Rd_len:', Rd_len

    ! allocate variables in veg structure
    allocate(tmp(Rd_len, vcmax_len, gm_len))
    allocate(LUT_VcmaxJmax(2, Rd_len, vcmax_len, gm_len))
    allocate(LUT_gm(gm_len))
    allocate(LUT_vcmax(vcmax_len))
    allocate(LUT_Rd(Rd_len))

    ok = nf90_inq_varid(ncid_gmlut, 'Vcmax', vcmax_id)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error inquiring variable Vcmax from LUT.')
    ok = nf90_inq_varid(ncid_gmlut, 'Jmax', jmax_id)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error inquiring variable Jmax from LUT.')

    ok = nf90_get_var(ncid_gmlut, vcmax_id, tmp)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error getting variable Vcmax from LUT.')
    LUT_VcmaxJmax(1, :, :, :) = tmp
    ok = nf90_get_var(ncid_gmlut, jmax_id, tmp)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error getting variable Jmax from LUT.')
    LUT_VcmaxJmax(2, :, :, :) = tmp
    ok = nf90_get_var(ncid_gmlut, gm_dimid, LUT_gm)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error getting dimension gm from LUT.')
    ok = nf90_get_var(ncid_gmlut, vcmax_dimid, LUT_vcmax)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error getting dimension Vcmax_Ci from LUT.')
    ok = nf90_get_var(ncid_gmlut, Rd_dimid, LUT_Rd)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error getting dimension Rd from LUT.')

    ! convert values to those used in CABLE
    LUT_VcmaxJmax = LUT_VcmaxJmax * 1.0e-06
    LUT_vcmax     = LUT_vcmax     * 1.0e-06
    LUT_Rd        = LUT_Rd        * 1.0e-06

    ok = nf90_close(ncid_gmlut)
    if (ok /= NF90_NOERR) call nc_abort(ok, 'Error closing gm lookup table.')

    deallocate(tmp)

  END SUBROUTINE read_gm_LUT


  SUBROUTINE find_Vcmax_Jmax_LUT(veg, p, LUT_VcmaxJmax, LUT_gm, LUT_vcmax, LUT_Rd)
    ! JK: note that LUT is single precision at the moment
    implicit none

    type(veg_parameter_type),     intent(inout) :: veg           ! vegetation parameters
    integer,                      intent(in)    :: p             ! veg type (tile)
    real(dp), dimension(:,:,:,:), intent(in)    :: LUT_VcmaxJmax ! Lookup table with Cc-based Vcmax and Jmax
    real(dp), dimension(:),       intent(in)    :: LUT_gm        ! gm values in gm LUT
    real(dp), dimension(:),       intent(in)    :: LUT_vcmax     ! Vcmax_ci values in gm LUT
    real(dp), dimension(:),       intent(in)    :: LUT_Rd        ! Rd values in gm LUT

    ! local
    logical :: val_ok        ! check for NAs
    integer :: maxit, i      ! maximum nr of iterations in while loop, loop counter
    integer :: igm, ivc, ird ! indices for LUT

    if (veg%frac4(p) < 0.001) then ! not C4
       ! determine current Ci-based values
       Rd        = veg%cfrd(p) * veg%vcmax(p) * light_inhibition(1200.0)
       gmmax25   = veg%gm(p)
       Vcmax25Ci = veg%vcmax(p)  ! LUT assumes a given Jmax/Vcmax ratio! see details in nc LUT

       ! determine right indices of LUT
       igm = minloc(abs(gmmax25 - LUT_gm), 1)
       ivc = minloc(abs(Vcmax25Ci - LUT_vcmax), 1)
       ird = minloc(abs(Rd - LUT_Rd), 1)

       i = 0
       maxit = 50
       val_ok = .false.

       do while (.not. val_ok .and. i .le. maxit)
          i = i + 1
          veg%vcmaxcc(p) = LUT_VcmaxJmax(1, ird, ivc, igm)
          veg%ejmaxcc(p) = LUT_VcmaxJmax(2, ird, ivc, igm)

          ! check for implausible parameter combinations that result in NAs
          if ((veg%vcmaxcc(p) > 0.0) .and. (veg%ejmaxcc(p) > 0.0)) then
             val_ok = .true.
          else
             write(64,*) "gm-Vcmax relationship does not work! Vcmax_ci:", Vcmax25Ci
             write(64,*) "iteration:", i
             igm = igm + 1
          endif
       end do
    else ! C4 (Vcmax and Jmax do not change with gm in C4 plants)
       veg%vcmaxcc(p) = veg%vcmax(p)
       veg%ejmaxcc(p) = veg%ejmax(p)
    endif

  END SUBROUTINE find_Vcmax_Jmax_LUT


  ! conversion of k Parameter in Collatz et al. 1992 from implicit gm model
  ! to explicit gm model.
  Subroutine adjust_k_Collatz(veg, p)

    implicit none

    type(veg_parameter_type), intent(inout) :: veg   ! vegetation parameters
    integer,                  intent(in)    :: p     ! vegetation type

    ! local
    integer  :: i,k
    integer  :: kmax = 1000
    integer  :: lAn
    real(dp) :: diff, diffx
    real(dp), dimension(nrcic4) :: An_Ci1, Ci1, Aj_Ci, Ae_Ci
    real(dp), dimension(:), allocatable :: An_Ci, An_Cc, Ci, Cc
    real(dp) :: kinc = 0.001_dp  ! increment of k

    if (veg%frac4(p) > 0.001) then ! C4
       Ci1       = (/(real(i,dp), i=1, nrcic4, 1)/) / 4.0_dp * 1.0e-6_dp
       Rd        = real(veg%cfrd(p) * veg%vcmax(p) * light_inhibition(1200.0), dp)
       gmmax25   = real(veg%gm(p), dp)
       Vcmax25Ci = real(veg%vcmax(p), dp)
       k25Ci     = real(veg%c4kci(p), dp)

       ! 1) calculate An-ci curves (no light limitation)
       Aj_Ci = Vcmax25Ci - Rd
       Ae_Ci = k25Ci * Ci1 - Rd

       An_Ci1 = min(Aj_Ci, Ae_Ci)

       ! 2) exclude negative An values and those not limited by Ci
       lAn = count((An_Ci1 > 0.0_dp) .and. (An_Ci1 == Ae_Ci))

       allocate(An_Ci(lAn))
       allocate(An_Cc(lAn))
       allocate(Ci(lAn))
       allocate(Cc(lAn))

       An_Ci = pack(An_Ci1, (An_Ci1 > 0.0_dp) .and. (An_Ci1 == Ae_Ci))
       Ci    = pack(Ci1, (An_Ci1 > 0.0_dp) .and. (An_Ci1 == Ae_Ci))

       ! 3) calculate Cc
       Cc = Ci - An_Ci / gmmax25

       ! 4) fit k to Cc-based model (a poor man's optimisation...)
       k     = 0
       diffx = 1.0e6_dp
       diff  = 0.0_dp
       k25Cc = k25Ci
       do while ((diff < diffx) .and. (k < kmax))
          if (k > 0) then
            diffx = diff
          endif
          An_Cc = k25Cc * Cc - Rd
          diff = sqrt(sum((An_Cc - An_Ci)**2) / lAn)
          k25Cc = k25Cc + kinc
          k = k + 1
       end do

       veg%c4kcc(p) = real(k25Cc - 2.0_dp*kinc) ! single precision
       ! subtract 2x the increment to get the right value
    else   ! C3 (not used)
      veg%c4kcc(p) = 0.0
    endif

  end subroutine adjust_k_Collatz

END MODULE cable_adjust_JV_gm_module
