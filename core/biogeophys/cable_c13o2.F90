!> \file cable_c13o2

!> \brief 13CO2 calculations in CABLE

!> \details Routines for calculating, shuffling, reading and writing 13C within Cable.

!> \author Matthias Cuntz, Juergen Knauer
!> \date Apr 2019

MODULE cable_c13o2

  implicit none

  private

  ! ------------------------------------------------------------------
  ! Public

  ! Photosynthesis
  public :: c13o2_init_flux     ! Initialise all 13C flux variables

  ! Casa
  public :: c13o2_init_pools    ! Initialise all 13C pools to 0
  public :: c13o2_save_casapool ! Save Casa pools before Casa update
  public :: c13o2_update_pools  ! Update 13C Casa pools

  ! LUC
  public :: c13o2_init_luc      ! Initialise all 13C LUC pools to 0
  public :: c13o2_save_luc      ! Save Casa pools, LUC pools, and land-use fractions
  public :: c13o2_update_luc    ! Update 13C Casa and LUC pools

  ! Output
  public :: c13o2_create_output ! Create output netcdf file
  public :: c13o2_write_output  ! Write to output netcdf file
  public :: c13o2_close_output  ! Close output netcdf file

  ! Restart
  public :: c13o2_read_restart_pools  ! Read 13CO2 restart_in_pools file
  public :: c13o2_write_restart_pools ! Write 13CO2 restart_out_pools file
  public :: c13o2_read_restart_luc    ! Read 13CO2 restart_in_luc file
  public :: c13o2_write_restart_luc   ! Write 13CO2 restart_out_luc file

  ! General
  public :: c13o2_print_delta_pools ! Print delta values of all Casa pools on screen
  public :: c13o2_print_delta_luc   ! Print delta values of all LUC pools on screen

  ! ------------------------------------------------------------------
  ! Private

  ! Casa
  private :: c13o2_fluxmatrix_pools ! Flux matrix between Casa pools
  private :: c13o2_sinks_pools      ! Sinks of isotope pool model
  private :: c13o2_sources_pools    ! Sources of isotope pool model
  private :: c13o2_sources_pools_nofrac ! Same but RAn=1, i.e. simply 12C sources
  private :: c13o2_save_c13o2pools  ! Save 13C Casa pools before update
  private :: c13o2_c13o2pools_back  ! Write 13C Casa pools to c13o2pools after update

  ! General
  private :: c13o2_err_handler      ! Internal error handler

  ! ------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------

  ! Public Routines

  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------

  ! Photosynthesis

  ! ------------------------------------------------------------------

  ! Initialise all 13CO2 flux variables
  subroutine c13o2_init_flux(met, c13o2flux)

    use cable_def_types_mod, only: dp => r_2
    use cable_def_types_mod, only: met_type
    use cable_c13o2_def,     only: c13o2_flux
    use mo_isotope,          only: vpdbc13

    implicit none

    type(met_type),   intent(inout) :: met
    type(c13o2_flux), intent(inout) :: c13o2flux

    c13o2flux%ca  = met%ca * vpdbc13
    c13o2flux%RAn = 1.0_dp ! vpdbc13 / vpdbc13 ! Divide by 13C so that about same numerical precision as 12C
    ! c13o2flux%Vstarch  = 0.0_dp
    ! c13o2flux%Rsucrose = vpdbc13
    ! c13o2flux%Rphoto   = vpdbc13
    ! c13o2flux%Rstarch  = vpdbc13

  end subroutine c13o2_init_flux

  ! ------------------------------------------------------------------

  ! Casa

  ! ------------------------------------------------------------------

  ! Initialise all 13CO2 pools to 0.
  subroutine c13o2_init_pools(casapool, casaflux, c13o2pools)

    use casavariable,        only: casa_pool, casa_flux
    use cable_c13o2_def,     only: c13o2_pool
    ! use mo_isotope,          only: vpdbc13

    implicit none

    type(casa_pool),  intent(in)    :: casapool
    type(casa_flux),  intent(in)    :: casaflux
    type(c13o2_pool), intent(inout) :: c13o2pools

    c13o2pools%cplant   = casapool%cplant   ! * vpdbc13 / vpdbc13 ! Divide by 13C
    c13o2pools%clitter  = casapool%clitter  ! * vpdbc13 / vpdbc13 ! so that about same numerical precision as 12C
    c13o2pools%csoil    = casapool%csoil    ! * vpdbc13 / vpdbc13
    c13o2pools%clabile  = casapool%clabile  ! * vpdbc13 / vpdbc13
    c13o2pools%charvest = casaflux%Charvest ! * vpdbc13 / vpdbc13

  end subroutine c13o2_init_pools

  ! ------------------------------------------------------------------

  ! Save old C concentrations of Casa pools
  subroutine c13o2_save_casapool(casapool, casasave)

    use cable_def_types_mod, only: dp => r_2
    use casavariable,        only: casa_pool

    implicit none

    type(casa_pool),          intent(in)  :: casapool
    real(dp), dimension(:,:), intent(out) :: casasave

    integer :: nplant, nlitter, nsoil, nstart, nend

    nplant  = size(casapool%cplant,2)
    nlitter = size(casapool%clitter,2)
    nsoil   = size(casapool%csoil,2)

    casasave = 0.0_dp

    ! plant
    nstart = 1
    nend   = nplant
    casasave(:,nstart:nend) = casapool%cplant
    ! litter
    nstart = nstart + nplant
    nend   = nend + nlitter
    casasave(:,nstart:nend) = casapool%clitter
    ! soil
    nstart = nstart + nlitter
    nend   = nend + nsoil
    casasave(:,nstart:nend) = casapool%csoil
    ! labile
    nstart = nstart + nsoil
    nend   = nstart
    casasave(:,nstart) = casapool%clabile

  end subroutine c13o2_save_casapool

  ! ------------------------------------------------------------------

  ! Calculate the generic isotope pool model, updating 13C concentrations in Casa pools
#ifdef C13DEBUG
  subroutine c13o2_update_pools(casasave, casaflux, c13o2flux, c13o2pools, casapool)
#else
  subroutine c13o2_update_pools(casasave, casaflux, c13o2flux, c13o2pools)
#endif

    use cable_def_types_mod,   only: dp => r_2
    use casadimension,         only: deltpool
    use casavariable,          only: casa_flux
    use casaparm,              only: leaf
    use cable_c13o2_def,       only: c13o2_flux, c13o2_pool
    use mo_isotope_pool_model, only: isotope_pool_model
#ifdef C13DEBUG
    use casavariable,          only: casa_pool
#endif

    implicit none

    real(dp), dimension(:,:), intent(in)    :: casasave
    type(casa_flux),          intent(in)    :: casaflux
    type(c13o2_flux),         intent(inout) :: c13o2flux
    type(c13o2_pool),         intent(inout) :: c13o2pools
#ifdef C13DEBUG
    type(casa_pool), optional,intent(in)    :: casapool
#endif

    real(dp), dimension(size(casasave,1),size(casasave,2)) :: c13o2save, casasources, casasinks
    real(dp), dimension(size(casasave,1),size(casasave,2),size(casasave,2)) :: fluxmatrix
#ifdef C13DEBUG
    real(dp), dimension(size(casasave,1),size(casasave,2)) :: casatmp, casasourcestmp
    real(dp) :: inew
    integer :: i, j, mp
#endif

    ! update 13C of cumulative harvest
    where (casasave(:,leaf) > 0.0_dp)
       c13o2pools%charvest = c13o2pools%charvest + &
            casaflux%FluxFromPtoHarvest * c13o2pools%cplant(:,leaf) / casasave(:,leaf)
    elsewhere
       c13o2pools%charvest = c13o2pools%charvest + casaflux%FluxFromPtoHarvest ! * vpdbc13 / vpdbc13
    endwhere

    ! save the old 13C concentration in casa pools
    call c13o2_save_c13o2pools(c13o2pools, c13o2save)

    ! Prepare isotope pool model by getting pool exchange fluxes, sources and sinks
    ! exchange fluxes between pools
    call c13o2_fluxmatrix_pools(casaflux, fluxmatrix)
    ! sources such as photosynthesis
    call c13o2_sources_pools(c13o2flux, casaflux, casasources)
    ! #ifdef C13DEBUG
    ! call c13o2_sources_pools_nofrac(c13o2flux, casaflux, casasourcestmp)
    ! #endif
    ! sinks such as respiration
    call c13o2_sinks_pools(casaflux, casasinks)

    ! #ifdef C13DEBUG
    ! ! Check C fluxes
    ! mp = size(casasave,1)
    ! do i=1, mp
    !    do j=1, 3 ! pools
    !       ! inew = casasave(i,j) + casaflux%fracCalloc(i,j) * casaflux%Cnpp(i) &
    !       !      - casaflux%FluxFromPtoCO2(i,j) - sum(casaflux%FluxFromPtoL(i,j,:))
    !       ! print*, 'Plant ', i, j, casapool%cplant(i,j)-inew
    !       ! inew = casasave(i,3+j) &
    !       !      - casaflux%FluxFromLtoCO2(i,j) + sum(casaflux%FluxFromPtoL(i,:,j)) - sum(casaflux%FluxFromLtoS(i,j,:))
    !       ! print*, 'Litter ', i, j, casapool%clitter(i,j)-inew
    !       inew = casasave(i,6+j) &
    !            - casaflux%FluxFromStoCO2(i,j) + sum(casaflux%FluxFromLtoS(i,:,j)) - sum(casaflux%FluxFromStoS(i,j,:)) &
    !            + sum(casaflux%FluxFromStoS(i,:,j))
    !       print*, '    Soil ', i, j, mydiff(casapool%csoil(i,j),inew)
    !       ! if (abs(casapool%csoil(i,j)-inew) > 1e-12_dp) then
    !       !    print*, 's01 ', casasave(i,6+j), casaflux%FluxFromStoCO2(i,j), sum(casaflux%FluxFromLtoS(i,:,j))
    !       !    print*, 's02 ', sum(casaflux%FluxFromStoS(i,j,:)), sum(casaflux%FluxFromStoS(i,:,j))
    !       ! endif
    !    end do
    ! end do
    ! #endif

    ! Calc the isotope pool model
    ! print*, '    DSoil11 ', mydiff(casasave(:,7:9), c13o2save(:,7:9))
    ! print*, '    Sour01 ', casasources(:,7:9)
    ! print*, '    Sink01 ', mydiff(casasinks(:,7:9), casaflux%FluxFromStoCO2)
    ! print*, '    Flux01 ',  sum(casaflux%FluxFromLtoS,2) + sum(casaflux%FluxFromStoS,2) - sum(fluxmatrix(:,4:9,7:9),2)
    ! print*, '    Flux02 ',  sum(casaflux%FluxFromStoS,3) - sum(fluxmatrix(:,7:9,4:9), dim=3)
    ! print*, '    iSour01 ', mydiff(sum(casaflux%FluxFromLtoS,2) + sum(casaflux%FluxFromStoS,2), sum(fluxmatrix(:,:,7:9),2))
    ! print*, '    iSink01 ', mydiff(sum(casaflux%FluxFromStoS,3), sum(fluxmatrix(:,7:9,:), dim=3))
    ! print*, 'Ci01 ', casapool%cplant, casapool%clitter, casapool%csoil, casapool%clabile
    call isotope_pool_model(deltpool, c13o2save, casasave, fluxmatrix, S=casasources, T=casasinks, trans=.true.)
    ! print*, 'Ci03 ', c13o2save
    ! #ifdef C13DEBUG
    ! print*, '    Diff1 plant ', mydiff(casapool%cplant,c13o2save(:,1:3))
    ! print*, '    Diff1 litter ', mydiff(casapool%clitter,c13o2save(:,4:6))
    ! print*, '    Diff1 soil ', mydiff(casapool%csoil,c13o2save(:,7:9))
    ! print*, '    Diff1 labile ', casapool%clabile-c13o2save(:,10)
    ! casatmp = casasave
    ! call isotope_pool_model(deltpool, casatmp, casasave, fluxmatrix, S=casasourcestmp, T=casasinks, trans=.true.)
    ! print*, '    DSoil12 ', mydiff(casatmp(:,7:9),c13o2save(:,7:9))
    ! print*, '    Diff2 plant ', mydiff(casapool%cplant,casatmp(:,1:3))
    ! print*, '    Diff2 litter ', mydiff(casapool%clitter,casatmp(:,4:6))
    ! print*, '    Diff2 soil ', mydiff(casapool%csoil,casatmp(:,7:9))
    ! print*, '    Diff2 labile ', casapool%clabile-casatmp(:,10)
    ! print*, '    Sources01 ', casasources
    ! print*, '    Sources01 ', casasourcestmp
    ! if (any(abs(c13o2save(:,1:3) - casatmp(:,1:3)) > c13o2save(:,1:3)*epsilon(1.0_dp)*10._dp)) then
    !    print*, '    DPlant01 ', casapool%cplant
    !    print*, '    DPlant02 ', c13o2save(:,1:3)
    !    print*, '    DPlant03 ', casatmp(:,1:3)
    !    print*, '    DPlant04 ', mydiff(casapool%cplant,c13o2save(:,1:3))
    ! endif
    ! if (any(abs(c13o2save(:,4:6) - casatmp(:,4:6)) > c13o2save(:,4:6)*epsilon(1.0_dp)*10._dp)) then
    !    print*, '    DLitter01 ', casapool%clitter
    !    print*, '    DLitter02 ', c13o2save(:,4:6)
    !    print*, '    DLitter03 ', casatmp(:,4:6)
    !    print*, '    DLitter04 ', casapool%clitter-c13o2save(:,4:6)
    ! endif
    ! if (any(abs(c13o2save(:,7:9) - casatmp(:,7:9)) > c13o2save(:,7:9)*epsilon(1.0_dp)*10._dp)) then
    !    print*, '    DSoil01 ', casapool%csoil
    !    print*, '    DSoil02 ', c13o2save(:,7:9)
    !    print*, '    DSoil03 ', casatmp(:,7:9)
    !    print*, '    DSoil04 ', mydiff(casapool%csoil,c13o2save(:,7:9))
    ! endif
    ! if (any(abs(c13o2save(:,10) - casatmp(:,10)) > c13o2save(:,10)*epsilon(1.0_dp)*10._dp)) then
    !    print*, '    DLabile01 ', casapool%clabile
    !    print*, '    DLabile02 ', c13o2save(:,10)
    !    print*, '    DLabile03 ', casatmp(:,10)
    !    print*, '    DLabile04 ', casapool%clabile-c13o2save(:,10)
    ! endif
    ! #endif

    ! put new solution into initial c13o2_pool type
    call c13o2_c13o2pools_back(c13o2pools, c13o2save)

  end subroutine c13o2_update_pools

  ! ------------------------------------------------------------------

  ! LUC

  ! ------------------------------------------------------------------

  ! Initialise all 13CO2 LUC pools to 0.
  subroutine c13o2_init_luc(c13o2luc, c13o2pools, veg, np)

    use cable_def_types_mod, only: veg_parameter_type, i4 => i_d, dp => r_2
    use cable_common_module, only: cable_user
    use casavariable,        only: casa_pool
    use casaparm,            only: leaf, wood, froot
    use cable_c13o2_def,     only: c13o2_luc, c13o2_pool, c13o2_alloc_luc, c13o2_zero_luc
    ! use mo_isotope,          only: vpdbc13

    implicit none

    type(c13o2_luc),          intent(inout) :: c13o2luc
    type(c13o2_pool),         intent(inout) :: c13o2pools
    type(veg_parameter_type), intent(in)    :: veg
    integer(i4),              intent(in)    :: np

    call c13o2_alloc_luc(c13o2luc, np)
    call c13o2_zero_luc(c13o2luc)

    if (cable_user%POPLUC_RunType == 'init') then
       ! zero biomass in secondary forest tiles (both CASA and POP variables)
       where (veg%iLU(:) == 2)
          c13o2pools%cplant(:,leaf)  = 0.01_dp ! * vpdbc13 / vpdbc13 ! Divide by 13C
          c13o2pools%cplant(:,wood)  = 0.01_dp ! * vpdbc13 / vpdbc13 ! so that about same numerical precision as 12C
          c13o2pools%cplant(:,froot) = 0.01_dp ! * vpdbc13 / vpdbc13
       end where
    else if (cable_user%POPLUC_RunType == 'restart') then
       call c13o2_read_restart_luc(cable_user%c13o2_restart_in_luc, c13o2luc)
    endif

  end subroutine c13o2_init_luc

  ! ------------------------------------------------------------------

  ! Save old C concentrations of Casa pools
  subroutine c13o2_save_luc(casapool, popluc, casasave, lucsave)

    use cable_def_types_mod, only: dp => r_2
    use casavariable,        only: casa_pool
    use popluc_types,        only: popluc_type

    implicit none

    type(casa_pool),          intent(in)  :: casapool
    type(popluc_type),        intent(in)  :: popluc
    real(dp), dimension(:,:), intent(out) :: casasave
    real(dp), dimension(:,:), intent(out) :: lucsave

    integer :: nharvest, nclearance, nstart, nend

    nharvest   = size(popluc%HarvProd,2)
    nclearance = size(popluc%ClearProd,2)

    call c13o2_save_casapool(casapool, casasave)

    lucsave = 0.0_dp

    ! harvest
    nstart = 1
    nend   = nharvest
    lucsave(:,nstart:nend) = popluc%HarvProd
    ! clearance
    nstart = nstart + nharvest
    nend   = nend + nclearance
    lucsave(:,nstart:nend) = popluc%ClearProd
    ! agric
    nstart = nstart + nclearance
    nend   = nstart
    lucsave(:,nstart) = popluc%AgProd

  end subroutine c13o2_save_luc

  ! ------------------------------------------------------------------

  ! Calculate the generic isotope pool model, updating 13C concentrations in Casa pools
#ifdef C13DEBUG
    subroutine c13o2_update_luc(casasave, lucsave, popluc, prim_only, c13o2pools, c13o2luc, casapool)
#else
    subroutine c13o2_update_luc(casasave, lucsave, popluc, prim_only, c13o2pools, c13o2luc)
#endif

    use cable_def_types_mod,   only: dp => r_2
    use cable_io_vars_module,  only: landpt, patch
    use casavariable,          only: casa_flux
    use popluc_types,          only: popluc_type
    use popluc_constants,      only: nLU
    use popluc_module,         only: kHarvProd, kClearProd, kAgProd
    use cable_c13o2_def,       only: c13o2_pool, c13o2_luc
    use mo_isotope_pool_model, only: isotope_pool_model
    use mo_isotope_luc_model,  only: isotope_luc_model
#ifdef C13DEBUG
    use casavariable,          only: casa_pool
    use mo_isotope,            only: delta1000
#endif

    implicit none

    real(dp), dimension(:,:), intent(in)    :: casasave
    real(dp), dimension(:,:), intent(in)    :: lucsave
    type(popluc_type),        intent(in)    :: popluc
    logical,  dimension(:),   intent(in)    :: prim_only
    type(c13o2_pool),         intent(inout) :: c13o2pools
    type(c13o2_luc),          intent(inout) :: c13o2luc
#ifdef C13DEBUG
    type(casa_pool), optional,intent(in)    :: casapool
#endif

    real(dp), dimension(size(casasave,1),size(casasave,2))   :: c13o2savepools, rsavepools ! 13C conc, isotope ratio
    real(dp), dimension(nLU)                                 :: A, Anew  ! patch fractions
    real(dp), dimension(nLU,nLU)                             :: dA, dAp  ! area change matrix
    real(dp), dimension(nLU)                                 :: c13sluc, sluc, tluc ! 13C source, C sink
    real(dp), dimension(c13o2luc%nharvest)                   :: c13sharv ! 13C source
    real(dp), dimension(c13o2luc%nharvest,c13o2luc%nharvest) :: fharv    ! fake harvest flux matrix
    real(dp), dimension(c13o2luc%nclearance)                     :: c13sclear ! 13C source
    real(dp), dimension(c13o2luc%nclearance,c13o2luc%nclearance) :: fclear    ! fake clearance flux matrix
    real(dp), dimension(1)                                       :: c13sag    ! 13C source
    real(dp), dimension(1,1)                                     :: fag       ! fake clearance flux matrix

    integer :: nplant, nlitter, nsoil, nharvest, nclearance
    integer :: g, j, l, c, cs, ce
#ifdef C13DEBUG
    real(dp), dimension(size(casasave,1),size(casasave,2)) :: casatmp
    real(dp), dimension(nLU)                               :: csluc
    integer :: iwtile, iwpool
#endif

    nplant     = c13o2pools%nplant
    nlitter    = c13o2pools%nlitter
    nsoil      = c13o2pools%nsoil
    nharvest   = c13o2luc%nharvest
    nclearance = c13o2luc%nclearance

    fharv  = 0.0_dp
    fclear = 0.0_dp
    fag    = 0.0_dp

    ! save old 13C concentrations in Casa and LUC
    call c13o2_save_c13o2pools(c13o2pools, c13o2savepools)
    ! isotope ratio in old pools
    rsavepools = 1.0_dp
    where(casasave > 0.0_dp) rsavepools = c13o2savepools / casasave

#ifdef C13DEBUG
    casatmp = casasave
#endif
    do g=1, popluc%np ! loop over popluc gridcells == Cable gridcells
       if (.not. prim_only(g)) then ! only then land-use change is possible
          j = landpt(g)%cstart ! start index of Cable tiles for grid cell g
          l = landpt(g)%cend   ! end index of Cable tiles for grid cell g
          if ((l-j+1) /= nLU) then
             write(*,*) 'I do still not understand the tiling: ', g, j, l, nLU
             stop 9
          endif
#ifdef C13DEBUG
          iwtile = l-2
          iwpool = 1
#endif
          if ((popluc%ptos(g) + popluc%ptog(g) + popluc%stog(g) + popluc%gtos(g)) > 0.0_dp) then
             A       = patch(j:l)%frac
             dA      = 0.0_dp
             dA(1,2) = popluc%ptos(g)
             dA(1,3) = popluc%ptog(g)
             dA(2,3) = popluc%stog(g)
             dA(3,2) = popluc%gtos(g)
             ! plant
             dAp  = 0.0
             Anew = A - sum(dA, dim=2) + sum(dA, dim=1)
             do c=1, nplant
#ifdef C13DEBUG
                if (c==iwpool) then
                   print*, 'LA01 ', A(iwtile), A(iwtile) - sum(dA(iwtile,:)) + sum(dA(:,iwtile)),  Anew(iwtile)
                   print*, 'LA02 ', sum(dA(:,iwtile))
                endif
#endif
                cs = c ! cs is index c in casasave
                tluc    = casasave(j:l,cs) * sum(dA, dim=2)
                sluc    = 0.0_dp
                c13sluc = 1.0_dp
                if (c==2) then
                   sluc(j+1)    = popluc%dcSHarvClear(g) * Anew(j+1)
                   c13sluc(j+1) = rsavepools(j+1,cs)
                endif
                call isotope_luc_model(c13o2pools%cplant(j:l,c), A, dAp, C=casasave(j:l,cs), S=sluc, Rs=c13sluc, T=tluc, At=Anew)
                ! update in popluc only done if Anew>1e-5
                where (Anew <= 1.e-5_dp) c13o2pools%cplant(j:l,c) = c13o2savepools(j:l,cs)
#ifdef C13DEBUG
                if (c==iwpool) then
                   print*, 'LP00 ', rsavepools(iwtile,cs)
                   print*, 'LP01 ', casasave(iwtile,cs)
                   print*, 'LP02 ', 0.0_dp
                   print*, 'LP03 ', tluc(iwtile) - sluc(iwtile), tluc(iwtile), sluc(iwtile)
                   print*, 'LP04 ', sum(dA(:,iwtile) * casasave(j:l,cs), dim=1) - sluc(iwtile) + tluc(iwtile)
                   print*, 'LP05 ', casasave(iwtile,cs)*sum(dA(:,iwtile))
                   ! Cnew = iC * (A - sum(dA, dim=2)) + sum(dA * spread(iC, dim=2, ncopies=nn), dim=1) + iS - iT
                   ! where (Anew > 0._dp) Cnew = Cnew / Anew
                   print*, 'LP06 ', ( casasave(iwtile,cs) * A(iwtile) + sluc(iwtile) - tluc(iwtile) ) / Anew(iwtile)
                endif
#endif
             end do
#ifdef C13DEBUG
             print*, 'LU01 ', delta1000(c13o2pools%cplant(j:l,:), casapool%cplant(j:l,:), 1.0_dp, -999._dp, tiny(1.0_dp))
#endif
             ! litter
             do c=1, nlitter
                cs = nplant + c
                c13sluc = 0.0_dp
#ifdef C13DEBUG
                csluc = 0.0_dp
#endif
                if (c==2) then
                   c13sluc(:) = sum(dA*spread(c13o2savepools(j:l,1), 2, nLU), dim=1) + &
                        sum(dA*spread(c13o2savepools(j:l,3), 2, nLU), dim=1)
#ifdef C13DEBUG
                   csluc(:) = sum(dA*spread(casasave(j:l,1), 2, nLU), dim=1) + &
                        sum(dA*spread(casasave(j:l,3), 2, nLU), dim=1)
#endif
                else if (c==3) then
                   c13sluc(2) = popluc%FluxPHarvResidtoLitter(g) * rsavepools(j,2) + &
                        popluc%FluxSHarvResidtoLitter(g) * rsavepools(j+1,2)
                   c13sluc(3) = popluc%FluxPClearResidtoLitter(g) * rsavepools(j,2) + &
                        popluc%FluxSClearResidtoLitter(g) * rsavepools(j+1,2)
#ifdef C13DEBUG
                   csluc(2) = popluc%FluxPHarvResidtoLitter(g) + &
                        popluc%FluxSHarvResidtoLitter(g)
                   csluc(3) = popluc%FluxPClearResidtoLitter(g) + &
                        popluc%FluxSClearResidtoLitter(g)
#endif
                endif
                call isotope_luc_model(c13o2pools%clitter(j:l,c), A, dA, C=casasave(j:l,cs), S=c13sluc)
                where (Anew <= 1.e-5_dp) c13o2pools%clitter(j:l,c) = c13o2savepools(j:l,cs)
#ifdef C13DEBUG
                if (c==iwpool) then
                   print*, 'LL00 ', rsavepools(j,2), rsavepools(j+1,2)
                   print*, 'LL01 ', casasave(iwtile,cs)
                   print*, 'LL02 ', sum(dA(:,iwtile) * casasave(j:l,cs))
                   print*, 'LL03 ', csluc(iwtile)
                   print*, 'LL04 ', sum(dA(:,iwtile) * casasave(j:l,cs), dim=1) + csluc(iwtile)
                   print*, 'LL05 ', casasave(iwtile,cs)*sum(dA(:,iwtile))
                   ! Cnew = iC * (A - sum(dA, dim=2)) + sum(dA * spread(iC, dim=2, ncopies=nn), dim=1) + iS - iT
                   ! where (Anew > 0._dp) Cnew = Cnew / Anew
                   print*, 'LL06 ', ( casasave(iwtile,cs) * (A(iwtile) - sum(dA(iwtile,:))) + &
                        sum(dA(:,iwtile) * casasave(j:l,cs)) + csluc(iwtile) ) / &
                        (A(iwtile) - sum(dA(iwtile,:)) + sum(dA(:,iwtile)))
                endif
                call isotope_luc_model(casatmp(j:l,cs), A, dA, C=casasave(j:l,cs), S=csluc)
#endif
             end do
#ifdef C13DEBUG
             print*, 'LU02 ', delta1000(c13o2pools%clitter(j:l,:), casapool%clitter(j:l,:), 1.0_dp, -999._dp, tiny(1.0_dp))
#endif
             ! soil
             do c=1, nsoil
                cs = nplant + nlitter + c
                call isotope_luc_model(c13o2pools%csoil(j:l,c), A, dA, C=casasave(j:l,cs))
                where (Anew <= 1.e-5_dp) c13o2pools%csoil(j:l,c) = c13o2savepools(j:l,cs)
#ifdef C13DEBUG
                if (c==iwpool) then
                   print*, 'LS01 ', casasave(iwtile,cs)
                   print*, 'LS02 ', sum(dA(:,iwtile) * casasave(j:l,cs))
                   print*, 'LS03 ', sum(dA(:,iwtile) * casasave(j:l,cs), dim=1)
                   print*, 'LS04 ', casasave(iwtile,cs)*sum(dA(:,iwtile))
                   ! Cnew = iC * (A - sum(dA, dim=2)) + sum(dA * spread(iC, dim=2, ncopies=nn), dim=1) + iS - iT
                   ! where (Anew > 0._dp) Cnew = Cnew / Anew
                   print*, 'LS05 ', ( casasave(iwtile,cs) * (A(iwtile) - sum(dA(iwtile,:))) + &
                        sum(dA(:,iwtile) * casasave(j:l,cs)) ) / &
                        (A(iwtile) - sum(dA(iwtile,:)) + sum(dA(:,iwtile)))
                endif
#endif
             end do
#ifdef C13DEBUG
             print*, 'LU03 ', delta1000(c13o2pools%csoil(j:l,:), casapool%csoil(j:l,:), 1.0_dp, -999._dp, tiny(1.0_dp))
#endif
             ! labile
             cs = nplant + nlitter + nsoil + 1
             call isotope_luc_model(c13o2pools%clabile(j:l), A, dA, C=casasave(j:l,cs))
             where (Anew <= 1.e-5_dp) c13o2pools%clabile(j:l) = c13o2savepools(j:l,cs)
#ifdef C13DEBUG
             print*, 'LU04 ', delta1000(c13o2pools%clabile(j:l), casapool%clabile(j:l), 1.0_dp, -999._dp, tiny(1.0_dp))
#endif
             ! harvest
             cs = 1
             ce = nharvest
             c13sharv = popluc%fracHarvProd(g,:) * sum(popluc%FHarvest(g,:) * rsavepools(j:l,2))
             call isotope_pool_model(1.0_dp, c13o2luc%charvest(g,:), lucsave(g,cs:ce), fharv, S=c13sharv, beta=kHarvProd)
#ifdef C13DEBUG
             print*, 'LU05 ', delta1000(c13o2luc%charvest(g,:), popluc%HarvProd(g,:), 1.0_dp, -999._dp, tiny(1.0_dp))
#endif
             ! clearance
             cs = nharvest + 1
             ce = nharvest + nclearance
             c13sclear = popluc%fracClearProd(g,:) * sum(popluc%FClearance(g,:) * rsavepools(j:l,2))
             call isotope_pool_model(1.0_dp, c13o2luc%cclearance(g,:), lucsave(g,cs:ce), fclear, S=c13sclear, beta=kClearProd)
#ifdef C13DEBUG
             print*, 'LU06 ', delta1000(c13o2luc%cclearance(g,:), popluc%ClearProd(g,:), 1.0_dp, -999._dp, tiny(1.0_dp))
#endif
             ! agric
             cs = nharvest + nclearance + 1
             ce = cs
             c13sag = c13o2pools%charvest(l) * patch(l)%frac
             call isotope_pool_model(1.0_dp, c13o2luc%cagric(g:g), lucsave(g,cs:ce), fag, S=c13sag, beta=(/kAgProd/))
#ifdef C13DEBUG
             print*, 'LU07 ', delta1000(c13o2luc%cagric(g:g), popluc%AgProd(g), 1.0_dp, -999._dp, tiny(1.0_dp))
#endif
          end if
       end if
       c13o2pools%charvest(l) = 0.0_dp
    end do

  end subroutine c13o2_update_luc

  ! ------------------------------------------------------------------

  ! Output

  ! ------------------------------------------------------------------

  ! Create 13C Casa output file
  subroutine c13o2_create_output(casamet, c13o2pools, file_id, vars, var_ids)

    use cable_common_module,  only: cable_user, filename
    use cable_io_vars_module, only: timeunits, calendar
    use casavariable,         only: casa_met, casafile
    use cable_c13o2_def,      only: c13o2_pool
    use netcdf,               only: nf90_create, nf90_clobber, nf90_noerr, &
         nf90_put_att, nf90_global, nf90_def_dim, nf90_unlimited, &
         nf90_def_var, nf90_double, nf90_int, nf90_enddef, nf90_put_var

    implicit none

    type(casa_met),                  intent(in)  :: casamet
    type(c13o2_pool),                intent(in)  :: c13o2pools
    integer,                         intent(out) :: file_id
    integer, parameter :: nvars = 7
    character(len=20), dimension(nvars), intent(out) :: vars
    integer,           dimension(nvars), intent(out) :: var_ids

    ! local variables
    integer :: i, status
    character(len=200) :: fname, dum
    integer :: olen

    ! dimensions
    integer, parameter :: ndims = 5
    character(len=20), dimension(ndims) :: dims    ! names
    integer,           dimension(ndims) :: dim_ids ! ids
    integer,           dimension(ndims) :: ldims   ! lengths

    ! variables
    character(len=40), dimension(nvars) :: lvars ! long names
    character(len=20), dimension(nvars) :: uvars ! units
    integer, dimension(nvars) :: dvars           ! number of dimensions
    integer, dimension(nvars) :: tvars           ! type
    integer, dimension(3)     :: idids           ! tmp for dim ids

    ! dimension names
    dims(1) = 'ntile'
    dims(2) = 'nplant'
    dims(3) = 'nlitter'
    dims(4) = 'nsoil'
    dims(5) = 'time'
    ! dimension lengths
    ldims(1) = c13o2pools%ntile
    ldims(2) = c13o2pools%nplant
    ldims(3) = c13o2pools%nlitter
    ldims(4) = c13o2pools%nsoil
    ldims(5) = nf90_unlimited

    ! variable names
    vars(1) = 'time'
    vars(2) = 'latitude'
    vars(3) = 'longitude'
    vars(4) = 'cplant'
    vars(5) = 'clitter'
    vars(6) = 'csoil'
    vars(7) = 'clabile'
    ! vars(8) = 'charvest'
    ! variable long_name
    lvars(1) = 'Time'
    lvars(2) = 'Latitude'
    lvars(3) = 'Longitude'
    lvars(4) = '13C content of plant pools'
    lvars(5) = '13C content of litter pools'
    lvars(6) = '13C content of soil pools'
    lvars(7) = '13C content of excess carbon pool'
    ! lvars(8) = '13C content of accumulated harvest'
    ! variable units
    uvars(1) = trim(timeunits)
    uvars(2) = 'degrees_north'
    uvars(3) = 'degrees_east'
    uvars(4) = 'kg(13C)/m^2'
    uvars(5) = 'kg(13C)/m^2'
    uvars(6) = 'kg(13C)/m^2'
    uvars(7) = 'kg(13C)/m^2'
    ! uvars(8) = 'kg(13C)/m^2'
    ! number of dimensions
    dvars(1) = 1 ! time
    dvars(2) = 1 ! ntile
    dvars(3) = 1 ! ntile
    dvars(4) = 3 ! ntile, nplant, time
    dvars(5) = 3 ! ntile, nlitter, time
    dvars(6) = 3 ! ntile, nsoil, time
    dvars(7) = 2 ! ntile, time
    ! dvars(8) = 2 ! ntile, time
    ! variable type
    tvars(1) = nf90_int
    tvars(2) = nf90_double
    tvars(3) = nf90_double
    tvars(4) = nf90_double
    tvars(5) = nf90_double
    tvars(6) = nf90_double
    tvars(7) = nf90_double
    ! tvars(8) = nf90_double

    ! output file name
    if (len_trim(cable_user%c13o2_outfile) > 0) then
       fname = trim(cable_user%c13o2_outfile)
    else
       if (len_trim(casafile%out) > 0) then
          olen = len_trim(casafile%out)
          fname = casafile%out(1:olen-3)//'_c13o2.nc'
       else
          if (len_trim(cable_user%mettype) > 0) then
             write(dum, fmt="(i4,'_',i4)") cable_user%yearstart, cable_user%yearend
             if ((cable_user%yearstart < 1000) .and. (cable_user%yearend < 1000)) then
                write(dum, fmt="(i3,'_',i3)") cable_user%yearstart, cable_user%yearend
             else if (cable_user%yearstart < 1000) then
                write(dum, fmt="(i3,'_',i4)") cable_user%yearstart, cable_user%yearend
             endif
             fname = trim(filename%path)//'/'//trim(cable_user%RunIden)//'_'//trim(dum)//'_casa_out_c13o2.nc'
          else
             fname = trim(filename%path)//'/'//trim(cable_user%RunIden)//'_casa_out_c13o2.nc'
          endif
       endif
    endif

    ! create output file
    write(*,*) 'Defining 13CO2 output file: ', trim(fname)
    status = nf90_create(trim(fname), nf90_clobber, file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not open c13o2 output file: '//trim(fname))
    ! status = nf90_redef(file_id)
    ! if (status /= nf90_noerr) &
    !      call c13o2_err_handler('Could not redef c13o2 output file: '//trim(fname))

    ! global attributes
    status = nf90_put_att(file_id, nf90_global, "StartYear", cable_user%yearstart)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not set global attribute StartYear in c13o2 output file: '//trim(fname))
    status = nf90_put_att(file_id, nf90_global, "EndYear", cable_user%yearend)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not set global attribute EndYear in c13o2 output file: '//trim(fname))
    status = nf90_put_att(file_id, nf90_global, "RunIden", cable_user%RunIden)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not set global attribute RunIden in c13o2 output file: '//trim(fname))

    ! define dimensions
    do i=1, ndims
       status = nf90_def_dim(file_id, trim(dims(i)), ldims(i), dim_ids(i))
       if (status /= nf90_noerr) then
          write(*,*) 'Dim: ', trim(dims(i)), ', Size: ', ldims(i)
          call c13o2_err_handler('Could not define dimension in c13o2 output file: '//trim(fname))
       endif
    end do

    ! define variables
    do i=1, nvars
       if (trim(vars(i)) == 'time') then
          idids(1) = dim_ids(5)
       else if (trim(vars(i)) == 'latitude') then
          idids(1) = dim_ids(1)
       else if (trim(vars(i)) == 'longitude') then
          idids(1) = dim_ids(1)
       else if (trim(vars(i)) == 'cplant') then
          idids(1) = dim_ids(1)
          idids(2) = dim_ids(2)
          idids(3) = dim_ids(5)
       else if (trim(vars(i)) == 'clitter') then
          idids(1) = dim_ids(1)
          idids(2) = dim_ids(3)
          idids(3) = dim_ids(5)
       else if (trim(vars(i)) == 'csoil') then
          idids(1) = dim_ids(1)
          idids(2) = dim_ids(4)
          idids(3) = dim_ids(5)
       else if (trim(vars(i)) == 'clabile') then
          idids(1) = dim_ids(1)
          idids(2) = dim_ids(5)
       ! else if (trim(vars(i)) == 'charvest') then
       !    idids(1) = dim_ids(1)
       !    idids(2) = dim_ids(5)
       else
          call c13o2_err_handler('Dimension not known for variable '//trim(vars(i))// &
               ' in c13o2 output file: '//trim(fname))
       endif
       status = nf90_def_var(file_id, trim(vars(i)), tvars(i), idids(1:dvars(i)), var_ids(i))
       if (status /= nf90_noerr) then
          write(*,*) 'Var: ', trim(vars(i)), ', dim_ids: ', idids(1:dvars(i))
          call c13o2_err_handler('Could not define variable in c13o2 output file: '//trim(fname))
       endif
       status = nf90_put_att(file_id, var_ids(i), 'long_name', lvars(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define long_name for '//trim(vars(i))// &
            ' in c13o2 output file: '//trim(fname))
       status = nf90_put_att(file_id, var_ids(i), 'units', uvars(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define units for '//trim(vars(i))//' in c13o2 output file: '//trim(fname))
       if (trim(vars(i)) == 'time') then
          status = nf90_put_att(file_id, var_ids(i), 'calendar', calendar)
          if (status /= nf90_noerr) &
               call c13o2_err_handler('Could not define calendar for variable time in c13o2 output file: '//trim(fname))
       endif
    end do ! nvars

    ! end definition phase
    status = nf90_enddef(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not end definition phase of c13o2 output file: '//trim(fname))

    ! put static variables lat and lon
    do i=1, nvars
       if (trim(vars(i)) == 'latitude') then
          status = nf90_put_var(file_id, var_ids(i), casamet%lat)
          if (status /= nf90_noerr) &
               call c13o2_err_handler('Could not put latitudes to c13o2 output file: '//trim(fname))
       else if (trim(vars(i)) == 'longitude') then
          status = nf90_put_var(file_id, var_ids(i), casamet%lon)
          if (status /= nf90_noerr) &
               call c13o2_err_handler('Could not put longitudes to c13o2 output file: '//trim(fname))
       endif
    end do

  end subroutine c13o2_create_output

  ! ------------------------------------------------------------------

  ! Write into 13C Casa output file
  subroutine c13o2_write_output(file_id, vars, var_ids, timestep, c13o2pools)

    use cable_c13o2_def, only: c13o2_pool
    use netcdf,          only: nf90_put_var, nf90_noerr

    implicit none

    integer,                         intent(in) :: file_id
    character(len=20), dimension(:), intent(in) :: vars
    integer,           dimension(:), intent(in) :: var_ids
    integer,                         intent(in) :: timestep
    type(c13o2_pool),                intent(in) :: c13o2pools

    ! local variables
    integer :: i, status
    integer :: nvars, nland, nplant, nlitter, nsoil
    integer, parameter :: sp = kind(1.0)

    nvars   = size(vars,1)
    nland   = c13o2pools%ntile
    nplant  = c13o2pools%nplant
    nlitter = c13o2pools%nlitter
    nsoil   = c13o2pools%nsoil
    ! define variables
    do i=1, nvars
       if (trim(vars(i)) == 'time') then
          status = nf90_put_var(file_id, var_ids(i), timestep, &
               start=(/timestep/))
       else if (trim(vars(i)) == 'cplant') then
          status = nf90_put_var(file_id, var_ids(i), real(c13o2pools%cplant,sp), &
               start=(/1,1,timestep/), count=(/nland,nplant,1/))
       else if (trim(vars(i)) == 'clitter') then
          status = nf90_put_var(file_id, var_ids(i), real(c13o2pools%clitter,sp), &
               start=(/1,1,timestep/), count=(/nland,nlitter,1/))
       else if (trim(vars(i)) == 'csoil') then
          status = nf90_put_var(file_id, var_ids(i), real(c13o2pools%csoil,sp), &
               start=(/1,1,timestep/), count=(/nland,nsoil,1/))
       else if (trim(vars(i)) == 'clabile') then
          status = nf90_put_var(file_id, var_ids(i), real(c13o2pools%clabile,sp), &
               start=(/1,timestep/), count=(/nland,1/))
       ! else if (trim(vars(i)) == 'charvest') then
       !    status = nf90_put_var(file_id, var_ids(i), real(c13o2pools%charvest,sp), &
       !         start=(/1,timestep/), count=(/nland,1/))
       endif
       if (status /= nf90_noerr) then
          write(*,*) 'Var: ', trim(vars(i)), ', var_id: ', var_ids(i)
          call c13o2_err_handler('Could not put variable in c13o2 output file')
       endif
    end do ! nvars

  end subroutine c13o2_write_output

  ! ------------------------------------------------------------------

  ! Close 13C Casa output file
  subroutine c13o2_close_output(file_id)

    use netcdf, only: nf90_close, nf90_noerr

    implicit none

    integer, intent(in) :: file_id

    integer :: status

    status = nf90_close(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not close 13CO2 output file.')

  end subroutine c13o2_close_output

  ! ------------------------------------------------------------------

  ! Restart

  ! ------------------------------------------------------------------

#ifndef UM_BUILD
  ! Read 13C Casa pools from restart file
  subroutine c13o2_read_restart_pools(c13o2_restart_in_pools, c13o2pools)

    use cable_c13o2_def,     only: c13o2_pool
    use netcdf,              only: nf90_open, nf90_nowrite, nf90_noerr, &
         nf90_inq_varid, nf90_get_var, nf90_close

    implicit none

    character(len=*), intent(in)    :: c13o2_restart_in_pools
    type(c13o2_pool), intent(inout) :: c13o2pools

    logical :: iexist
    integer :: status, file_id, var_id

    inquire(file=trim(c13o2_restart_in_pools), exist=iexist)
    if (.not. iexist) &
         call c13o2_err_handler('13CO2 restart_in_pools file does not exist: '//trim(c13o2_restart_in_pools))

    write(*,*) 'Read 13CO2 restart_in_pools file: ', trim(c13o2_restart_in_pools)

    status = nf90_open(trim(c13o2_restart_in_pools), nf90_nowrite, file_id)
    if (status /= nf90_noerr) call c13o2_err_handler('Error 13CO2 restart_in_pools file: '//trim(c13o2_restart_in_pools))

    !MC13 ToDo - dimension check
    ! cplant
    status = nf90_inq_varid(file_id, 'cplant', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable cplant not found in restart_in_pools file: '//trim(c13o2_restart_in_pools))
    status = nf90_get_var(file_id, var_id, c13o2pools%cplant, &
         start=(/1,1/), count=(/c13o2pools%ntile, c13o2pools%nplant/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable cplant from restart_in_pools file: '//trim(c13o2_restart_in_pools))
    ! clitter
    status = nf90_inq_varid(file_id, 'clitter', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable clitter not found in restart_in_pools file: '//trim(c13o2_restart_in_pools))
    status = nf90_get_var(file_id, var_id, c13o2pools%clitter, &
         start=(/1,1/), count=(/c13o2pools%ntile, c13o2pools%nlitter/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable clitter from restart_in_pools file: '//trim(c13o2_restart_in_pools))
    ! csoil
    status = nf90_inq_varid(file_id, 'csoil', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable csoil not found in restart_in_pools file: '//trim(c13o2_restart_in_pools))
    status = nf90_get_var(file_id, var_id, c13o2pools%csoil, &
         start=(/1,1/), count=(/c13o2pools%ntile, c13o2pools%nsoil/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable csoil from restart_in_pools file: '//trim(c13o2_restart_in_pools))
    ! clabile
    status = nf90_inq_varid(file_id, 'clabile', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable clabile not found in restart_in_pools file: '//trim(c13o2_restart_in_pools))
    status = nf90_get_var(file_id, var_id, c13o2pools%clabile, start=(/1/), count=(/c13o2pools%ntile/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable clabile from restart_in_pools file: '//trim(c13o2_restart_in_pools))
    ! charvest
    status = nf90_inq_varid(file_id, 'charvest', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable charvest not found in restart_in_pools file: '//trim(c13o2_restart_in_pools))
    status = nf90_get_var(file_id, var_id, c13o2pools%charvest, start=(/1/), count=(/c13o2pools%ntile/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable charvest from restart_in_pools file: '//trim(c13o2_restart_in_pools))

    status = nf90_close(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not close the c13o2 restart_in_pools file: '//trim(c13o2_restart_in_pools))

  end subroutine c13o2_read_restart_pools
#endif

  ! ------------------------------------------------------------------

#ifndef UM_BUILD
  ! Write 13C Casa pools into restart file
  subroutine c13o2_write_restart_pools(c13o2pools)

    use cable_common_module, only: cable_user, filename, CurYear
    use cable_c13o2_def,     only: c13o2_pool
    use netcdf,              only: nf90_create, nf90_clobber, nf90_noerr, &
         nf90_put_att, nf90_global, nf90_def_dim, &
         nf90_def_var, nf90_double, nf90_enddef, nf90_put_var, nf90_close

    implicit none

    type(c13o2_pool), intent(in) :: c13o2pools

    ! local variables
    integer :: file_id, land_id, plant_id, litter_id, soil_id
    integer :: i, status
    character(len=200) :: fname, cyear

    ! 1 dim arrays (npt,)
    integer, parameter :: nlandvars = 2
    character(len=20), dimension(nlandvars)  :: landvars
    ! 2 dim arrays (npt,nplant)
    integer, parameter :: nplantvars = 1
    character(len=20), dimension(nplantvars)  :: plantvars
    ! 2 dim arrays (npt,nlitter)
    integer, parameter :: nlittervars = 1
    character(len=20), dimension(nlittervars) :: littervars
    ! 2 dim arrays (npt,nsoil)
    integer, parameter :: nsoilvars = 1
    character(len=20), dimension(nsoilvars)   :: soilvars
    ! variable ids
    integer, dimension(nlandvars)   :: landvars_id
    integer, dimension(nplantvars)  :: plantvars_id
    integer, dimension(nlittervars) :: littervars_id
    integer, dimension(nsoilvars)   :: soilvars_id

    ! Variables
    landvars(1)   = 'clabile'
    landvars(2)   = 'charvest'
    plantvars(1)  = 'cplant'
    littervars(1) = 'clitter'
    soilvars(1)   = 'csoil'

    ! restart file name
    if (len_trim(cable_user%c13o2_restart_out_pools) > 0) then
       fname = trim(cable_user%c13o2_restart_out_pools)
    else
       fname = trim(filename%path)//'/'//trim(cable_user%RunIden)//'_c13o2_pools_rst.nc'
    endif

    ! create restart file
    write(*,*) 'Writing 13CO2 Casa restart file: ', trim(fname)
    status = nf90_create(trim(fname), nf90_clobber, file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not open c13o2 restart_out_pools file: '//trim(fname))
    ! status = nf90_redef(file_id)
    ! if (status /= nf90_noerr) &
    !      call c13o2_err_handler('Could not redef c13o2 restart_out_pools file: '//trim(fname))

    ! global attributes
    write(cyear, fmt='(I4)') CurYear + 1
    status = nf90_put_att(file_id, nf90_global, "Valid restart date", "01.01."//cyear)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not set global date attribute in c13o2 restart_out_pools file: '//trim(fname))

    ! define dimensions
    status = nf90_def_dim(file_id, 'nland', c13o2pools%ntile, land_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define land dimension in c13o2 restart_out_pools file: '//trim(fname))
    status = nf90_def_dim(file_id, 'nplant', c13o2pools%nplant, plant_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define plant dimension in c13o2 restart_out_pools file: '//trim(fname))
    status = nf90_def_dim(file_id, 'nlitter', c13o2pools%nlitter, litter_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define litter dimension in c13o2 restart_out_pools file: '//trim(fname))
    status = nf90_def_dim(file_id, 'nsoil', c13o2pools%nsoil, soil_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define soil dimension in c13o2 restart_out_pools file: '//trim(fname))

    ! define variables
    do i=1, nlandvars
       status = nf90_def_var(file_id, trim(landvars(i)), nf90_double, &
            (/land_id/), landvars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(landvars(i))// &
            ' in c13o2 restart_out_pools file: '//trim(fname))
    end do
    do i=1, nplantvars
       status = nf90_def_var(file_id, trim(plantvars(i)), nf90_double, &
            (/land_id,plant_id/), plantvars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(plantvars(i))// &
            ' in c13o2 restart_out_pools file: '//trim(fname))
    end do
    do i=1, nlittervars
       status = nf90_def_var(file_id, trim(littervars(i)), nf90_double, &
            (/land_id,litter_id/), littervars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(littervars(i))// &
            ' in c13o2 restart_out_pools file: '//trim(fname))
    end do
    do i=1, nsoilvars
       status = nf90_def_var(file_id, trim(soilvars(i)), nf90_double, &
            (/land_id,soil_id/), soilvars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(soilvars(i))// &
            ' in c13o2 restart_out_pools file: '//trim(fname))
    end do

    status = nf90_enddef(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not end definition phase of c13o2 restart_out_pools file: '//trim(fname))

    ! put variables
    status = nf90_put_var(file_id, landvars_id(1), c13o2pools%clabile)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(landvars(1))//' to c13o2 restart_out_pools file: '//trim(fname))
    status = nf90_put_var(file_id, landvars_id(2), c13o2pools%charvest)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(landvars(2))//' to c13o2 restart_out_pools file: '//trim(fname))
    status = nf90_put_var(file_id, plantvars_id(1), c13o2pools%cplant)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(plantvars(1))//' to c13o2 restart_out_pools file: '//trim(fname))
    status = nf90_put_var(file_id, littervars_id(1), c13o2pools%clitter)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(littervars(1))//' to c13o2 restart_out_pools file: '//trim(fname))
    status = nf90_put_var(file_id, soilvars_id(1), c13o2pools%csoil)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(soilvars(1))//' to c13o2 restart_out_pools file: '//trim(fname))

    ! close restart file
    status = nf90_close(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not close c13o2 restart_out_pools file: '//trim(fname))

  end subroutine c13o2_write_restart_pools
#endif

  ! ------------------------------------------------------------------

  ! Read 13C LUC pools from restart file
  subroutine c13o2_read_restart_luc(c13o2_restart_in_luc, c13o2luc)

    use cable_c13o2_def,     only: c13o2_luc
    use netcdf,              only: nf90_open, nf90_nowrite, nf90_noerr, &
         nf90_inq_varid, nf90_get_var, nf90_close

    implicit none

    character(len=*), intent(in)    :: c13o2_restart_in_luc
    type(c13o2_luc),  intent(inout) :: c13o2luc

    logical :: iexist
    integer :: status, file_id, var_id

    inquire(file=trim(c13o2_restart_in_luc), exist=iexist)
    if (.not. iexist) &
         call c13o2_err_handler('13CO2 restart_in_luc file does not exist: '//trim(c13o2_restart_in_luc))

    write(*,*) 'Read 13CO2 restart_in_luc file: ', trim(c13o2_restart_in_luc)

    status = nf90_open(trim(c13o2_restart_in_luc), nf90_nowrite, file_id)
    if (status /= nf90_noerr) call c13o2_err_handler('Error 13CO2 restart_in_luc file: '//trim(c13o2_restart_in_luc))

    !MC13 ToDo - dimension check
    ! charvest
    status = nf90_inq_varid(file_id, 'charvest', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable charvest not found in restart_in_luc file: '//trim(c13o2_restart_in_luc))
    status = nf90_get_var(file_id, var_id, c13o2luc%charvest, &
         start=(/1,1/), count=(/c13o2luc%nland, c13o2luc%nharvest/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable charvest from restart_in_luc file: '//trim(c13o2_restart_in_luc))
    ! cclearance
    status = nf90_inq_varid(file_id, 'cclearance', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable cclearance not found in restart_in_luc file: '//trim(c13o2_restart_in_luc))
    status = nf90_get_var(file_id, var_id, c13o2luc%cclearance, &
         start=(/1,1/), count=(/c13o2luc%nland, c13o2luc%nclearance/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable cclearance from restart_in_luc file: '//trim(c13o2_restart_in_luc))
    ! cagric
    status = nf90_inq_varid(file_id, 'cagric', var_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Variable cagric not found in restart_in_luc file: '//trim(c13o2_restart_in_luc))
    status = nf90_get_var(file_id, var_id, c13o2luc%cagric, &
         start=(/1/), count=(/c13o2luc%nland/))
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Error reading variable cagric from restart_in_luc file: '//trim(c13o2_restart_in_luc))

    status = nf90_close(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not close the c13o2 restart_in_luc file: '//trim(c13o2_restart_in_luc))

  end subroutine c13o2_read_restart_luc

  ! ------------------------------------------------------------------

  ! Write 13C Casa pools into restart file
  subroutine c13o2_write_restart_luc(c13o2luc)

    use cable_common_module, only: cable_user, filename, CurYear
    use cable_c13o2_def,     only: c13o2_luc
    use netcdf,              only: nf90_create, nf90_clobber, nf90_noerr, &
         nf90_put_att, nf90_global, nf90_def_dim, &
         nf90_def_var, nf90_double, nf90_enddef, nf90_put_var, nf90_close

    implicit none

    type(c13o2_luc), intent(in) :: c13o2luc

    ! local variables
    integer :: file_id, land_id, harvest_id, clearance_id
    integer :: i, status
    character(len=200) :: fname, cyear

    ! 1 dim arrays (npt,)
    integer, parameter :: nlandvars = 1
    character(len=20), dimension(nlandvars)      :: landvars
    ! 2 dim arrays (npt,nharvest)
    integer, parameter :: nharvestvars = 1
    character(len=20), dimension(nharvestvars)   :: harvestvars
    ! 2 dim arrays (npt,nclearance)
    integer, parameter :: nclearancevars = 1
    character(len=20), dimension(nclearancevars) :: clearancevars
    ! variable ids
    integer, dimension(nlandvars)      :: landvars_id
    integer, dimension(nharvestvars)   :: harvestvars_id
    integer, dimension(nclearancevars) :: clearancevars_id

    ! Variables
    landvars(1)      = 'cagric'
    harvestvars(1)   = 'charvest'
    clearancevars(1) = 'cclearance'

    ! restart file name
    if (len_trim(cable_user%c13o2_restart_out_luc) > 0) then
       fname = trim(cable_user%c13o2_restart_out_luc)
    else
       fname = trim(filename%path)//'/'//trim(cable_user%RunIden)//'_c13o2_luc_rst.nc'
    endif

    ! create restart file
    write(*,*) 'Writing 13CO2 LUC restart file: ', trim(fname)
    status = nf90_create(trim(fname), nf90_clobber, file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not open c13o2 restart_out_luc file: '//trim(fname))
    ! status = nf90_redef(file_id)
    ! if (status /= nf90_noerr) &
    !      call c13o2_err_handler('Could not redef c13o2 restart_out_luc file: '//trim(fname))

    ! global attributes
    write(cyear, fmt='(I4)') CurYear + 1
    status = nf90_put_att(file_id, nf90_global, "Valid restart date", "01.01."//cyear)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not set global date attribute in c13o2 restart_out_luc file: '//trim(fname))

    ! define dimensions
    status = nf90_def_dim(file_id, 'nland', c13o2luc%nland, land_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define land dimension in c13o2 restart_out_luc file: '//trim(fname))
    status = nf90_def_dim(file_id, 'nharvest', c13o2luc%nharvest, harvest_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define harvest dimension in c13o2 restart_out_luc file: '//trim(fname))
    status = nf90_def_dim(file_id, 'nclearance', c13o2luc%nclearance, clearance_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not define clearance dimension in c13o2 restart_out_luc file: '//trim(fname))

    ! define variables
    do i=1, nlandvars
       status = nf90_def_var(file_id, trim(landvars(i)), nf90_double, &
            (/land_id/), landvars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(landvars(i))// &
            ' in c13o2 restart_out_luc file: '//trim(fname))
    end do
    do i=1, nharvestvars
       status = nf90_def_var(file_id, trim(harvestvars(i)), nf90_double, &
            (/land_id,harvest_id/), harvestvars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(harvestvars(i))// &
            ' in c13o2 restart_out_luc file: '//trim(fname))
    end do
    do i=1, nclearancevars
       status = nf90_def_var(file_id, trim(clearancevars(i)), nf90_double, &
            (/land_id,clearance_id/), clearancevars_id(i))
       if (status /= nf90_noerr) &
            call c13o2_err_handler('Could not define variable '//trim(clearancevars(i))// &
            ' in c13o2 restart_out_luc file: '//trim(fname))
    end do

    status = nf90_enddef(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not end definition phase of c13o2 restart_out_luc file: '//trim(fname))

    ! put variables
    status = nf90_put_var(file_id, landvars_id(1), c13o2luc%cagric)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(landvars(1))//' to c13o2 restart_out_luc file: '//trim(fname))
    status = nf90_put_var(file_id, harvestvars_id(1), c13o2luc%charvest)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(harvestvars(1))//' to c13o2 restart_out_luc file: '//trim(fname))
    status = nf90_put_var(file_id, clearancevars_id(1), c13o2luc%cclearance)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could put variable '//trim(clearancevars(1))//' to c13o2 restart_out_luc file: '//trim(fname))

    ! close restart file
    status = nf90_close(file_id)
    if (status /= nf90_noerr) &
         call c13o2_err_handler('Could not close c13o2 restart_out_luc file: '//trim(fname))

  end subroutine c13o2_write_restart_luc

  ! ------------------------------------------------------------------

  ! Generell

  ! ------------------------------------------------------------------

  ! Print 13C delta values of Casa pools on screen
  subroutine c13o2_print_delta_pools(casapool, casaflux, c13o2pools)

    use cable_def_types_mod, only: dp => r_2
    use casavariable,        only: casa_pool, casa_flux
    use cable_c13o2_def,     only: c13o2_pool
    use mo_isotope,          only: delta1000!, vpdbc13

    implicit none

    type(casa_pool),  intent(in) :: casapool
    type(casa_flux),  intent(in) :: casaflux
    type(c13o2_pool), intent(in) :: c13o2pools

    integer :: nplant, nlitter, nsoil, mp
    character(len=30) :: form1

    mp      = size(casapool%cplant,1)
    nplant  = size(casapool%cplant,2)
    nlitter = size(casapool%clitter,2)
    nsoil   = size(casapool%csoil,2)

    write(*,*) '    delta-13C of Casa pools'
    !write(form1,'(A,I3,A)') '(a,', nplant*mp, 'f20.14)'
    !write(form1,'(A,I3,A)') '(a,', nplant*mp, 'es22.14)'
    write(form1,'(A,I3,A)') '(a,', nplant*mp, 'g15.6e3)'
    ! write(*,*) '        plant12:    ', casapool%cplant
    ! write(*,*) '        plant13:    ', c13o2pools%cplant
    ! write(*,form1) '        plant:      ', delta1000(c13o2pools%cplant,   casapool%cplant,   vpdbc13, -999._dp, tiny(1.0_dp))
    write(*,form1) '        plant:      ', delta1000(c13o2pools%cplant,   casapool%cplant,   1.0_dp, -999._dp, tiny(1.0_dp))
    write(form1,'(A,I3,A)') '(a,', nlitter*mp, 'g15.6e3)'
    ! write(*,*) '        litter12:   ', casapool%clitter
    ! write(*,*) '        litter13:   ', c13o2pools%clitter
    write(*,form1) '        litter:     ', delta1000(c13o2pools%clitter,  casapool%clitter,  1.0_dp, -999._dp, tiny(1.0_dp))
    write(form1,'(A,I3,A)') '(a,', nsoil*mp, 'g15.6e3)'
    ! write(*,*) '        soil12:     ', casapool%csoil
    ! write(*,*) '        soil13:     ', c13o2pools%csoil
    write(*,form1) '        soil:       ', delta1000(c13o2pools%csoil,    casapool%csoil,    1.0_dp, -999._dp, tiny(1.0_dp))
    write(form1,'(A,I3,A)') '(a,', 1*mp, 'g15.6e3)'
    ! write(*,*) '        labile12:   ', casapool%clabile
    ! write(*,*) '        labile13:   ', c13o2pools%clabile
    write(*,form1) '        labile:     ', delta1000(c13o2pools%clabile,  casapool%clabile,  1.0_dp, -999._dp, tiny(1.0_dp))
    write(form1,'(A,I3,A)') '(a,', 1*mp, 'g15.6e3)'
    ! write(*,*) '        Charvest12: ', casaflux%Charvest
    ! write(*,*) '        Charvest13: ', c13o2pools%charvest
    write(*,form1) '        Charvest:   ', delta1000(c13o2pools%charvest, casaflux%Charvest, 1.0_dp, -999._dp, tiny(1.0_dp))

  end subroutine c13o2_print_delta_pools

  ! ------------------------------------------------------------------

  ! Print 13C delta values of LUC pools on screen
  subroutine c13o2_print_delta_luc(popluc, c13o2luc)

    use cable_def_types_mod, only: dp => r_2
    use popluc_types,        only: popluc_type
    use cable_c13o2_def,     only: c13o2_luc
    use mo_isotope,          only: delta1000!, vpdbc13

    implicit none

    type(popluc_type), intent(in) :: popluc
    type(c13o2_luc),   intent(in) :: c13o2luc

    integer :: nharvest, nclearance, mp
    character(len=30) :: form1

    mp         = size(popluc%HarvProd,1)
    nharvest   = size(popluc%HarvProd,2)
    nclearance = size(popluc%ClearProd,2)

    write(*,*) '    delta-13C of LUC pools'
    write(form1,'(A,I3,A)') '(a,', nharvest*mp, 'g15.6e3)'
    ! write(*,form1) '        harvest:   ', delta1000(c13o2luc%charvest,   popluc%HarvProd,  vpdbc13, -999._dp, tiny(1.0_dp))
    write(*,form1) '        harvest:   ', delta1000(c13o2luc%charvest,   popluc%HarvProd,  1.0_dp, -999._dp, tiny(1.0_dp))
    write(form1,'(A,I3,A)') '(a,', nclearance*mp, 'g15.6e3)'
    write(*,form1) '        clearance: ', delta1000(c13o2luc%cclearance, popluc%ClearProd, 1.0_dp, -999._dp, tiny(1.0_dp))
    write(form1,'(A,I3,A)') '(a,', 1*mp, 'g15.6e3)'
    write(*,form1) '        agric:     ', delta1000(c13o2luc%cagric,     popluc%AgProd,    1.0_dp, -999._dp, tiny(1.0_dp))

  end subroutine c13o2_print_delta_luc

  ! ------------------------------------------------------------------

  ! Private Routines

  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------

  ! Casa

  ! ------------------------------------------------------------------

  ! Exchange matrix of all fluxes between Casa pools
  subroutine c13o2_fluxmatrix_pools(casaflux, fluxmatrix)

    use cable_def_types_mod, only: dp => r_2
    use casavariable,        only: casa_flux

    implicit none

    type(casa_flux),            intent(in)  :: casaflux
    real(dp), dimension(:,:,:), intent(out) :: fluxmatrix

    integer :: nplant, nlitter, nsoil, nstartr, nendr, nstartc, nendc

    nplant  = size(casaflux%FluxFromPtoL,2)
    nlitter = size(casaflux%FluxFromPtoL,3)
    nsoil   = size(casaflux%FluxFromLtoS,3)

    fluxmatrix = 0.0_dp

    ! plant to litter
    nstartr = 1
    nendr   = nplant
    nstartc = nplant + 1
    nendc   = nplant + nlitter
    fluxmatrix(:,nstartr:nendr,nstartc:nendc) = casaflux%FluxFromPtoL
    ! litter to soil
    nstartr = nstartr + nplant
    nendr   = nendr + nlitter
    nstartc = nstartc + nlitter
    nendc   = nendc + nsoil
    fluxmatrix(:,nstartr:nendr,nstartc:nendc) = casaflux%FluxFromLtoS
    ! soil to soil
    nstartr = nstartr + nlitter
    nendr   = nendr + nsoil
    nstartc = nstartc
    nendc   = nendc
    fluxmatrix(:,nstartr:nendr,nstartc:nendc) = casaflux%FluxFromStoS
    ! labile = 0.

  end subroutine c13o2_fluxmatrix_pools

  ! ------------------------------------------------------------------

  ! Sinks of Casa pools others than exchanges to other Casa pools, such as respiration
  subroutine c13o2_sinks_pools(casaflux, casasinks)

    use cable_def_types_mod, only: dp => r_2
    use casavariable,        only: casa_flux
    use casaparm,            only: leaf

    implicit none

    type(casa_flux),          intent(in)  :: casaflux
    real(dp), dimension(:,:), intent(out) :: casasinks

    integer :: nplant, nlitter, nsoil, nstart, nend

    nplant  = size(casaflux%FluxFromPtoCO2,2)
    nlitter = size(casaflux%FluxFromLtoCO2,2)
    nsoil   = size(casaflux%FluxFromStoCO2,2)

    casasinks = 0.0_dp

    ! plant
    nstart = 1
    nend   = nplant
    casasinks(:,nstart:nend) = casaflux%FluxFromPtoCO2
    casasinks(:,leaf)        = casasinks(:,leaf) + casaflux%FluxFromPtoHarvest
    ! litter
    nstart = nstart + nplant
    nend   = nend + nlitter
    casasinks(:,nstart:nend) = casaflux%FluxFromLtoCO2
    ! soil
    nstart = nstart + nlitter
    nend   = nend + nsoil
    casasinks(:,nstart:nend) = casaflux%FluxFromStoCO2
    ! labile
    nstart = nstart + nsoil
    nend   = nstart
    casasinks(:,nstart) = casaflux%clabloss

  end subroutine c13o2_sinks_pools

  ! ------------------------------------------------------------------

  ! Sources of Casa pools others than exchanges to other Casa pools, such as from photosynthesis
  subroutine c13o2_sources_pools(c13o2flux, casaflux, casasources)

    use cable_def_types_mod, only: dp => r_2
    use cable_c13o2_def,     only: c13o2_flux
    use casavariable,        only: casa_flux

    implicit none

    type(c13o2_flux),         intent(inout) :: c13o2flux
    type(casa_flux),          intent(in)    :: casaflux
    real(dp), dimension(:,:), intent(out)   :: casasources

    integer :: nplant, nlitter, nsoil, nstart, nend

    nplant  = size(casaflux%fracCalloc,2)
    nlitter = size(casaflux%FluxFromLtoCO2,2)
    nsoil   = size(casaflux%FluxFromStoCO2,2)

    ! Stay with old isotope ratio if net assimilation < 0. because no NSC pool
    where (abs(c13o2flux%cAn12) > 0.0_dp) c13o2flux%RAn = c13o2flux%cAn / c13o2flux%cAn12

    casasources = 0.0_dp

    ! Use leaf discrimination for net assimilation and GPP
    ! plant
    casasources(:,1:nplant) = casaflux%fracCalloc * spread(casaflux%Cnpp * c13o2flux%RAn, 2, nplant)
    ! litter = 0.
    ! soil = 0.
    ! labile
    nstart = nplant + nlitter + nsoil + 1
    nend   = nstart
    casasources(:,nstart) = casaflux%fracClabile * casaflux%Cgpp * c13o2flux%RAn

  end subroutine c13o2_sources_pools

  ! Same but with c13o2flux%RAn=1 - to check isotope code
  subroutine c13o2_sources_pools_nofrac(c13o2flux, casaflux, casasources)

    use cable_def_types_mod, only: dp => r_2
    use cable_c13o2_def,     only: c13o2_flux
    use casavariable,        only: casa_flux

    implicit none

    type(c13o2_flux),         intent(in)  :: c13o2flux
    type(casa_flux),          intent(in)  :: casaflux
    real(dp), dimension(:,:), intent(out) :: casasources

    integer :: nplant, nlitter, nsoil, nstart, nend

    nplant  = size(casaflux%fracCalloc,2)
    nlitter = size(casaflux%FluxFromLtoCO2,2)
    nsoil   = size(casaflux%FluxFromStoCO2,2)

    casasources = 0.0_dp

    ! Use leaf discrimination for net assimilation and GPP
    ! plant - spread(x, dim, ncopies)
    casasources(:,1:nplant) = casaflux%fracCalloc * spread(casaflux%Cnpp * 1.0_dp, 2, nplant)
    ! litter = 0.
    ! soil = 0.
    ! labile
    nstart = nplant + nlitter + nsoil + 1
    nend   = nstart
    casasources(:,nstart) = casaflux%fracClabile * casaflux%Cgpp * 1.0_dp

  end subroutine c13o2_sources_pools_nofrac

  ! ------------------------------------------------------------------

  ! Save old 13C concentration of Casa pools
  subroutine c13o2_save_c13o2pools(c13o2pools, c13o2save)

    use cable_def_types_mod, only: dp => r_2
    use cable_c13o2_def,     only: c13o2_pool

    implicit none

    type(c13o2_pool),         intent(in)  :: c13o2pools
    real(dp), dimension(:,:), intent(out) :: c13o2save

    integer :: nplant, nlitter, nsoil, nstart, nend

    nplant  = c13o2pools%nplant
    nlitter = c13o2pools%nlitter
    nsoil   = c13o2pools%nsoil

    ! plant
    nstart = 1
    nend   = nplant
    c13o2save(:,nstart:nend) = c13o2pools%cplant
    ! litter
    nstart = nstart + nplant
    nend   = nend + nlitter
    c13o2save(:,nstart:nend) = c13o2pools%clitter
    ! soil
    nstart = nstart + nlitter
    nend   = nend + nsoil
    c13o2save(:,nstart:nend) = c13o2pools%csoil
    ! labile
    nstart = nstart + nsoil
    nend   = nstart
    c13o2save(:,nstart) = c13o2pools%clabile

  end subroutine c13o2_save_c13o2pools

  ! ------------------------------------------------------------------

  ! Write back new 13C concentration of Casa pools c13o2_pool type
  subroutine c13o2_c13o2pools_back(c13o2pools, c13o2save)

    use cable_def_types_mod, only: dp => r_2
    use cable_c13o2_def,     only: c13o2_pool

    implicit none

    type(c13o2_pool),         intent(inout) :: c13o2pools
    real(dp), dimension(:,:), intent(in)    :: c13o2save

    integer :: nplant, nlitter, nsoil, nstart, nend

    nplant  = c13o2pools%nplant
    nlitter = c13o2pools%nlitter
    nsoil   = c13o2pools%nsoil

    ! plant
    nstart = 1
    nend   = nplant
    c13o2pools%cplant = c13o2save(:,nstart:nend)
    ! litter
    nstart = nstart + nplant
    nend   = nend + nlitter
    c13o2pools%clitter = c13o2save(:,nstart:nend)
    ! soil
    nstart = nstart + nlitter
    nend   = nend + nsoil
    c13o2pools%csoil = c13o2save(:,nstart:nend)
    ! labile
    nstart = nstart + nsoil
    nend   = nstart
    c13o2pools%clabile = c13o2save(:,nstart)

  end subroutine c13o2_c13o2pools_back

  ! ------------------------------------------------------------------

  ! General

  ! ------------------------------------------------------------------

  ! Write out error message and stop program
  subroutine c13o2_err_handler(message)

    implicit none

    character(len=*), intent(in) :: message

    write(*,*) trim(message)
    stop 9

  end subroutine c13o2_err_handler

  ! give diff in order of epsilon
  elemental pure function mydiff(v1, v2)

    use mo_kind, only: dp

    implicit none

    real(dp), intent(in) :: v1, v2
    real(dp)             :: mydiff

    ! integer  :: n
    real(dp) :: nn

    if (v1 > 0.0_dp) then
       nn = log10(v1)
       mydiff = (v1-v2) * 10._dp**(-nn) / epsilon(1.0_dp)
       ! n  = int(nn)
       ! mydiff = (v1-v2) * 10._dp**(-n) / epsilon(1.0_dp)
    else
       mydiff = (v1-v2) / epsilon(1.0_dp)
    endif

    return

  end function mydiff


END MODULE cable_c13o2
