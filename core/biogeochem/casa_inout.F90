!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Input and output code for CASA-CNP when run offline
!          ACCESS version may use some of this code but split into different files?
!
! Contact: Yingping.Wang@csiro.au and Bernard.Pak@csiro.au
!
! History: Developed for offline code.  Expect to re-write for MPI and ACCESS
!          versions
!
! Matthias Cuntz, 24/2/2020: Make all new netCDF files 64bit_offset. Do also not reopen
! the define mode because they are already in define mode once they are created. Read
! and write N/P variables from and to dump files only if icycle>1 or icycle>2. Also pass
! the variables in MPI code only if icycle>1,2. Rewrote write_casa_output_nc so that N/P
! variables are only written in case of icycle>1,2. Initialised more variables in Cable and Casa,
! which might otherwise be undefined during write_output. Also some more consistencies of variable
! kinds during calculations, mostly casa. Reinstated ESM-SnowMIP output variables that were present
! in an earlier branch trunk4252_wales. Only pass bal to close_output_file, not all
! other unnecessary variables. Use a generic conversion routine toreal4() for output to catch
! floating point underflow in conversions (real(variable,4)).
! ==============================================================================
! casa_inout.f90
!
! the following routines are used when "casacnp" is coupled to "cable"
!   zero_casa_dump
!   casa_readbiome
!   casa_readphen
!   casa_readpoint   (removed, now done in parameter_module)
!   casa_init
!   casa_poolout
!   casa_cnpflux  (zeros casabal quantites on doy 1 and updates casabal at end of biogeochem)
!   biogeochem
module casa_inout

  implicit none

contains

  subroutine zero_casa_dump(casamet, phen)
    ! Zero variables that will be read from dump files.
    ! This is done to have zero in MPI and non MPI code
    ! if .not. cable_user%casa_dump_read.

    use cable_def_types_mod, only: r_2
    use casavariable,        only: casa_met
    use phenvariable,        only: phen_variable
    use cable_common_module, only: cable_user

    implicit none

    type(casa_met),      intent(inout) :: casamet
    type(phen_variable), intent(inout) :: phen

    casamet%Tairkspin      = 0.0_r_2
    casamet%Tsoilspin_1    = 0.0_r_2
    casamet%Tsoilspin_2    = 0.0_r_2
    casamet%Tsoilspin_3    = 0.0_r_2
    casamet%Tsoilspin_4    = 0.0_r_2
    casamet%Tsoilspin_5    = 0.0_r_2
    casamet%Tsoilspin_6    = 0.0_r_2
    casamet%moistspin_1    = 0.0_r_2
    casamet%moistspin_2    = 0.0_r_2
    casamet%moistspin_3    = 0.0_r_2
    casamet%moistspin_4    = 0.0_r_2
    casamet%moistspin_5    = 0.0_r_2
    casamet%moistspin_6    = 0.0_r_2
    casamet%cgppspin       = 0.0_r_2
    casamet%crmplantspin_1 = 0.0_r_2
    casamet%crmplantspin_2 = 0.0_r_2
    casamet%crmplantspin_3 = 0.0_r_2
    phen%phasespin      = 0
    phen%doyphasespin_1 = 0
    phen%doyphasespin_2 = 0
    phen%doyphasespin_3 = 0
    phen%doyphasespin_4 = 0
    casamet%mtempspin = 0.0_r_2
    casamet%frecspin  = 0.0_r_2
    if (cable_user%c13o2) then
       casamet%cAn12spin = 0.0_r_2
       casamet%cAn13spin = 0.0_r_2
    end if
    if (cable_user%call_blaze) then
       casamet%dprecip_spin      = 0.0_r_2
       casamet%aprecip_av20_spin = 0.0_r_2
       casamet%du10_max_spin     = 0.0_r_2
       casamet%drhum_spin        = 0.0_r_2
       casamet%dtemp_max_spin    = 0.0_r_2
       casamet%dtemp_max_spin    = 0.0_r_2
       casamet%KBDI_spin         = 0.0_r_2
       casamet%D_MacArthur_spin  = 0.0_r_2
       casamet%FFDI_spin         = 0.0_r_2
       casamet%DSLR_spin         = 0
       casamet%last_precip_spin  = 0.0_r_2
    end if

  end subroutine zero_casa_dump

  ! #define UM_BUILD YES
  SUBROUTINE casa_readbiome(veg,soil,casabiome,casapool,casaflux,casamet,phen)
    ! mst actually not used in this routine (BP sep2010)
    !SUBROUTINE casa_readbiome(mvt,mst,veg,soil, &
    !                          casabiome,casapool,casaflux,casamet,phen)
    USE cable_def_types_mod
    USE casadimension
    USE casaparm
    USE casavariable
    USE phenvariable
    !! vh_js !!
    USE cable_common_module, only: cable_user

    IMPLICIT NONE
    !  INTEGER,               INTENT(IN)    :: mvt,mst
    TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    TYPE (casa_biome),          INTENT(INOUT) :: casabiome
    TYPE (casa_pool),           INTENT(INOUT) :: casapool
    TYPE (casa_flux),           INTENT(INOUT) :: casaflux
    TYPE (casa_met),            INTENT(INOUT) :: casamet
    TYPE (phen_variable),       INTENT(INOUT) :: phen

    ! local variables
    REAL(r_2), DIMENSION(mvtype)       :: leafage,frootage,woodage
    REAL(r_2), DIMENSION(mvtype)       :: cwdage,metage,strage
    REAL(r_2), DIMENSION(mvtype)       :: micage,slowage,passage,clabileage,slax
    REAL(r_2), DIMENSION(mvtype,mplant):: ratioCNplant
    REAL(r_2), DIMENSION(mvtype,msoil) :: ratioCNsoil,ratioCNsoilmin,ratioCNsoilmax
    REAL(r_2), DIMENSION(ms)           :: depthsoila,depthsoilb
    REAL(r_2), DIMENSION(mvtype)       :: xfNminloss, xfNminleach, xnfixrate
    REAL(r_2), DIMENSION(mvtype)       :: cleaf,cwood,cfroot,      &
         cmet,cstr,ccwd,          &
         cmic,cslow,cpass
    REAL(r_2), DIMENSION(mvtype)       :: nleaf,nwood,nfroot,      &
         nmet,nstr,ncwd,          &
         nmic,nslow,npass,xnsoilmin
    REAL(r_2), DIMENSION(mvtype)       :: xpleaf, xpwood, xpfroot, &
         xpmet, xpstr, xpcwd,     &
         xpmic,xpslow,xppass,xplab,xpsorb,xpocc
    REAL(r_2), DIMENSION(mso)       :: xkmlabp,xpsorbmax,xfPleach
    REAL(r_2), DIMENSION(mso,msoil) :: ratioNPsoil
    REAL(r_2), DIMENSION(mvtype)       :: xfherbivore,xxkleafcoldmax, xxkleafdrymax
    REAL(r_2), DIMENSION(mvtype,ms)    :: fracroot
    REAL(r_2) ::  xratioNPleafmin,xratioNPleafmax,         &
         xratioNPwoodmin,xratioNPwoodmax,         &
         xratioNPfrootmin,xratioNPfrootmax
    INTEGER :: i,iv1,nv,ns,npt,iso
    INTEGER :: nv0,nv1,nv2,nv3,nv4,nv6,nv7,nv8,nv9,nv10,nv11,nv12,nv13
    REAL(r_2), DIMENSION(mvtype)       :: xxnpmax,xq10soil,xxkoptlitter,xxkoptsoil,xprodptase, &
         xcostnpup,xmaxfinelitter,xmaxcwd,xnintercept,xnslope


    REAL(r_2), DIMENSION(mvtype)       :: la_to_sa, vcmax_scalar, disturbance_interval
    REAL(r_2), DIMENSION(mvtype)       :: DAMM_EnzPool, DAMM_KMO2,DAMM_KMcp, DAMM_Ea, DAMM_alpha
    REAL(r_2), DIMENSION(mso)          :: xxkplab,xxkpsorb,xxkpocc


    OPEN(101,file=casafile%cnpbiome)
    DO i=1,3
       READ(101,*)
    ENDDO

    DO nv=1,mvtype
       READ(101,*) nv0,casabiome%ivt2(nv)
    ENDDO

    READ(101,*)
    READ(101,*)
    DO nv=1,mvtype
       READ(101,*) nv1,casabiome%kroot(nv),casabiome%rootdepth(nv),      &
            casabiome%kuptake(nv),casabiome%krootlen(nv),         &
            casabiome%kminN(nv), casabiome%kuplabP(nv),           &
            xfherbivore(nv),leafage(nv),woodage(nv),frootage(nv), &
            metage(nv),strage(nv),cwdage(nv),  &
            micage(nv),slowage(nv),passage(nv),clabileage(nv),slax(nv)
       ! write(59,*) nv, leafage(nv),woodage(nv),frootage(nv), &
       !             metage(nv),strage(nv),cwdage(nv),  &
       !             micage(nv),slowage(nv),passage(nv)
    ENDDO

    READ(101,*)
    READ(101,*)
    DO nv=1,mvtype
       READ(101,*) nv2, &
            casabiome%fracnpptoP(nv,leaf),casabiome%fracnpptoP(nv,wood), &
            casabiome%fracnpptoP(nv,froot),casabiome%rmplant(nv,leaf),   &
            casabiome%rmplant(nv,wood),casabiome%rmplant(nv,froot)
    ENDDO

    READ(101,*)
    READ(101,*)
    DO nv=1,mvtype
       READ(101,*) nv2, ratioCNplant(nv,leaf),ratioCNplant(nv,wood),   &
            ratioCNplant(nv,froot),                                         &
            casabiome%ftransNPtoL(nv,leaf), casabiome%ftransNPtoL(nv,wood), &
            casabiome%ftransNPtoL(nv,froot),                                &
            casabiome%fracligninplant(nv,leaf),                             &
            casabiome%fracligninplant(nv,wood),                             &
            casabiome%fracligninplant(nv,froot),                            &
            ratioCNsoil(nv,mic),ratioCNsoil(nv,slow),ratioCNsoil(nv,pass),  &
            ratioCNsoilmin(nv,mic),ratioCNsoilmin(nv,slow),ratioCNsoilmin(nv,pass),  &
            ratioCNsoilmax(nv,mic),ratioCNsoilmax(nv,slow),ratioCNsoilmax(nv,pass),  &
                                ! xfherbivore(nv),casabiome%ratiofrootleaf(nv),                  &
            casabiome%glaimax(nv),casabiome%glaimin(nv)
    ENDDO

    READ(101,*)
    READ(101,*)
    DO nv=1,mvtype
       READ(101,*) nv3, cleaf(nv),cwood(nv),cfroot(nv),cmet(nv),   &
            cstr(nv),ccwd(nv), cmic(nv), cslow(nv),cpass(nv)
    ENDDO

    READ(101,*)
    READ(101,*)
    DO nv=1,mvtype
       READ(101,*) nv4, &
            phen%TKshed(nv),xxkleafcoldmax(nv),casabiome%xkleafcoldexp(nv), &
            xxkleafdrymax(nv),casabiome%xkleafdryexp(nv)
    ENDDO
    !  READ(101,*)
    !  READ(101,*)
    !  DO nv=1,mvtype
    !    READ(101,*) nv5, &
    !         xxkleafcoldmax(nv),casabiome%xkleafcoldexp(nv),   &
    !         xxkleafdrymax(nv),casabiome%xkleafdryexp(nv)
    !  ENDDO

    READ(101,*)
    READ(101,*)
    DO nv=1,mvtype
       READ(101,*) nv6, &
            casabiome%ratioNCplantmin(nv,leaf),casabiome%ratioNCplantmax(nv,leaf), &
            casabiome%ratioNCplantmin(nv,wood),casabiome%ratioNCplantmax(nv,wood), &
            casabiome%ratioNCplantmin(nv,froot),casabiome%ratioNCplantmax(nv,froot), &
            xfNminloss(nv), xfNminleach(nv),xnfixrate(nv)
    ENDDO

    READ(101,*)
    READ(101,*)
    DO nv=1,mvtype
       READ(101,*) nv7,nleaf(nv),nwood(nv),nfroot(nv), &
            nmet(nv),nstr(nv), ncwd(nv), &
            nmic(nv),nslow(nv),npass(nv),xnsoilmin(nv)
    ENDDO

    READ(101,*)
    READ(101,*)
    DO nv=1,mvtype
       READ(101,*) nv8,xratioNPleafmin,xratioNPleafmax,      &
            xratioNPwoodmin,xratioNPwoodmax,                      &
            xratioNPfrootmin,xratioNPfrootmax,                    &
            casabiome%ftransPPtoL(nv,leaf), casabiome%ftransPPtoL(nv,wood), &
            casabiome%ftransPPtoL(nv,froot)
       casabiome%ratioPcplantmin(nv,leaf)  = 1.0/(xratioNPleafmax*ratioCNplant(nv,leaf))
       casabiome%ratioPcplantmax(nv,leaf)  = 1.0/(xratioNPleafmin*ratioCNplant(nv,leaf))
       casabiome%ratioPcplantmin(nv,wood)  = 1.0/(xratioNPwoodmax*ratioCNplant(nv,wood))
       casabiome%ratioPcplantmax(nv,wood)  = 1.0/(xratioNPwoodmin*ratioCNplant(nv,wood))
       casabiome%ratioPcplantmin(nv,froot) = 1.0/(xratioNPfrootmax*ratioCNplant(nv,froot))
       casabiome%ratioPcplantmax(nv,froot) = 1.0/(xratioNPfrootmin*ratioCNplant(nv,froot))


       casabiome%ratioNPplantmin(nv,leaf)  = xratioNPleafmin
       casabiome%ratioNPplantmax(nv,leaf)  = xratioNPleafmax
       casabiome%ratioNPplantmin(nv,wood)  = xratioNPwoodmin
       casabiome%ratioNPplantmax(nv,wood)  = xratioNPwoodmax
       casabiome%ratioNPplantmin(nv,froot) = xratioNPfrootmin
       casabiome%ratioNPplantmax(nv,froot) = xratioNPfrootmax

    ENDDO

    READ(101,*)
    READ(101,*)
    DO iso=1,mso
       READ(101,*) nv9,xkmlabp(iso),xpsorbmax(iso),xfPleach(iso), &
            ratioNPsoil(iso,mic),ratioNPsoil(iso,slow),ratioNPsoil(iso,pass), &
            xxkplab(iso),xxkpsorb(iso),xxkpocc(iso)
    ENDDO
    READ(101,*)
    READ(101,*)
    DO nv=1,mvtype
       READ(101,*) nv10, &
            xpleaf(nv),xpwood(nv),xpfroot(nv),xpmet(nv),xpstr(nv),xpcwd(nv), &
            xpmic(nv),xpslow(nv),xppass(nv),xplab(nv),xpsorb(nv),xpocc(nv)
    ENDDO

    READ(101,*)
    READ(101,*)
    DO nv=1,mvtype
       READ(101,*) nv11, &
            xxnpmax(nv),xq10soil(nv),xxkoptlitter(nv),xxkoptsoil(nv),xprodptase(nv), &
            xcostnpup(nv),xmaxfinelitter(nv),xmaxcwd(nv),xnintercept(nv),xnslope(nv)
    ENDDO

    READ(101,*)
    READ(101,*)
    DO nv=1,mvtype
       READ(101,*) nv12, &
            la_to_sa(nv),disturbance_interval(nv),vcmax_scalar(nv)
    ENDDO

    READ(101,*)
    READ(101,*)
    DO nv=1,mvtype
       READ(101,*) nv13, &
            DAMM_EnzPool(nv), DAMM_KMO2(nv),DAMM_KMcp(nv), DAMM_Ea(nv), DAMM_alpha(nv)
    ENDDO

    CLOSE(101)

    fracroot   = 0.0_r_2
    depthsoila = 0.0_r_2
    depthsoilb = 0.0_r_2
    DO ns=1, ms
       depthsoilb(ns) = depthsoilb(ns) + soil%zse(ns)
       IF (ns==1) THEN
          depthsoila(ns) = 0.0_r_2
       ELSE
          depthsoila(ns) = depthsoilb(ns-1)
       ENDIF
    ENDDO

    DO nv=1,mvtype
       casabiome%sla(nv)             = slax(nv)
       casabiome%fraclabile(nv,leaf) = deltcasa*0.6_r_2    !1/day
       casabiome%fraclabile(nv,froot)= deltcasa*0.4_r_2    !1/day
       casabiome%fraclabile(nv,wood) = deltcasa*0.0_r_2
       casabiome%plantrate(nv,leaf)  = deltcasa/(leafage(nv)*(1.0_r_2-xfherbivore(nv)))
       casabiome%plantrate(nv,froot) = deltcasa/frootage(nv)
       casabiome%plantrate(nv,wood)  = deltcasa/woodage(nv)
       casabiome%litterrate(nv,metb) = deltcasa/metage(nv)
       casabiome%litterrate(nv,str)  = deltcasa/strage(nv)
       casabiome%litterrate(nv,cwd)  = deltcasa/cwdage(nv)
       casabiome%soilrate(nv,mic)    = deltcasa/micage(nv)
       casabiome%soilrate(nv,slow)   = deltcasa/slowage(nv)
       casabiome%soilrate(nv,pass)   = deltcasa/passage(nv)
       casabiome%xkleafcoldmax(nv)   = deltcasa * xxkleafcoldmax(nv)
       casabiome%xkleafdrymax(nv)    = deltcasa * xxkleafdrymax(nv)
       ! casabiome%kuplabp(nv)         = xkuplabp(nv)
       casabiome%rmplant(nv,:)       = casabiome%rmplant(nv,:)*deltcasa
       casabiome%kclabrate(nv)       = deltcasa/clabileage(nv)

       casabiome%xnpmax(nv)          = xxnpmax(nv)
       casabiome%q10soil(nv)         = xq10soil(nv)
       casabiome%xkoptlitter(nv)     = xxkoptlitter(nv)
       casabiome%xkoptsoil(nv)       = xxkoptsoil(nv)
       casabiome%prodptase(nv)       = xprodptase(nv)/365.0_r_2   ! convert from yearly to daily
       casabiome%costnpup(nv)        = xcostnpup(nv)
       casabiome%maxfinelitter(nv)   = xmaxfinelitter(nv)
       casabiome%maxcwd(nv)          = xmaxcwd(nv)
       casabiome%nintercept(nv)      = xnintercept(nv)
       casabiome%nslope(nv)          = xnslope(nv)

       casabiome%la_to_sa(nv)        = la_to_sa(nv)
       casabiome%vcmax_scalar(nv)        = vcmax_scalar(nv)
       casabiome%disturbance_interval(nv)        = disturbance_interval(nv)
       casabiome%DAMM_EnzPool(nv) = DAMM_EnzPool(nv)
       casabiome%DAMM_KMO2(nv) = DAMM_KMO2(nv)
       casabiome%DAMM_KMcp(nv) = DAMM_KMcp(nv)
       casabiome%DAMM_Ea(nv) = DAMM_Ea(nv)
       casabiome%DAMM_alpha(nv) = DAMM_alpha(nv)
    ENDDO

    DO ns=1,mso
       casabiome%xkplab(ns)          =  xxkplab(ns)
       casabiome%xkpsorb(ns)         =  xxkpsorb(ns)
       casabiome%xkpocc(ns)          =  xxkpocc(ns)
    ENDDO

    DO npt = 1, mp
       iv1=veg%iveg(npt)
       iso=casamet%isorder(npt)
       ! The following to be commented out when coupled to CABLE
       !    veg%froot(npt,:) = fracroot(iv1,:)
       casamet%iveg2(npt) =casabiome%ivt2(iv1)
       casamet%lnonwood(npt) = 1
       casapool%cplant(npt,wood)  = 0.0_r_2
       casapool%clitter(npt,cwd)  = 0.0_r_2
       casapool%nplant(npt,wood)  = 0.0_r_2
       casapool%nlitter(npt,cwd)  = 0.0_r_2
       casapool%pplant(npt,wood)  = 0.0_r_2
       casapool%plitter(npt,cwd)  = 0.0_r_2
       IF (casamet%iveg2(npt)==forest.or.casamet%iveg2(npt)==shrub) THEN
          casamet%lnonwood(npt) = 0
          casapool%cplant(npt,wood)  = Cwood(iv1)
          casapool%clitter(npt,cwd)  = ccwd(iv1)
          casapool%nplant(npt,wood)  = nwood(iv1)
          casapool%nlitter(npt,cwd)  = ncwd(iv1)
          casapool%pplant(npt,wood)  = xpwood(iv1)
          casapool%plitter(npt,cwd)  = xpcwd(iv1)
          !! vh_js !!
          IF (cable_user%CALL_POP) THEN  ! initialise very small wood pool, so POP can start from zero.
             casapool%cplant(npt,wood) = 0.01_r_2
             casapool%nplant(npt,wood)= casabiome%ratioNCplantmin(iv1,wood)* casapool%cplant(npt,wood)
             casapool%pplant(npt,wood)= casabiome%ratioPCplantmin(iv1,wood)* casapool%cplant(npt,wood)
          ENDIF
          !! vh_js
       ENDIF
       casapool%cplant(npt,leaf)     = cleaf(iv1)
       casapool%cplant(npt,froot)    = cfroot(iv1)
       casapool%clabile(npt)         = 0.0_r_2
       casapool%clitter(npt,metb)    = cmet(iv1)
       casapool%clitter(npt,str)     = cstr(iv1)
       casapool%csoil(npt,mic)       = cmic(iv1)
       casapool%csoil(npt,slow)      = cslow(iv1)
       casapool%csoil(npt,pass)      = cpass(iv1)
       IF (icycle==1) THEN
          casapool%ratioNCplant(npt,:)  = 1.0_r_2/ratioCNplant(iv1,:)
       ENDIF

       ! initializing glai in case not reading pool file (eg. during spin)
       casamet%glai(npt) = MAX(casabiome%glaimin(iv1), &
            casabiome%sla(iv1) * casapool%cplant(npt,leaf))

       casaflux%fNminloss(npt)   = xfNminloss(iv1)
       ! comment out by ypw 12/07/2009
       casaflux%fNminleach(npt)  = 10.0_r_2*xfNminleach(iv1) * deltcasa
       ! casaflux%fNminleach(npt)  = xfNminleach(iv1)
       casapool%nplant(npt,leaf) = nleaf(iv1)
       casapool%nplant(npt,froot)= nfroot(iv1)
       casapool%nlitter(npt,metb) = nmet(iv1)
       ! casapool%nlitter(npt,str) = nstr(iv1)
       casapool%nlitter(npt,str) = cstr(iv1)*ratioNCstrfix
       casapool%nsoil(npt,mic)   = nmic(iv1)
       casapool%nsoil(npt,slow)  = nslow(iv1)
       casapool%nsoil(npt,pass)  = npass(iv1)
       casapool%nsoilmin(npt)    = xnsoilmin(iv1)
       casapool%pplant(npt,leaf) = xpleaf(iv1)
       casapool%pplant(npt,froot)= xpfroot(iv1)
       casapool%plitter(npt,metb) = xpmet(iv1)
       ! casapool%plitter(npt,str) = xpstr(iv1)
       casapool%plitter(npt,str) = casapool%nlitter(npt,str)/ratioNPstrfix
       casapool%psoil(npt,mic)   = xpmic(iv1)
       casapool%psoil(npt,slow)  = xpslow(iv1)
       casapool%psoil(npt,pass)  = xppass(iv1)
       casapool%psoillab(npt)    = xplab(iv1)
       casapool%psoilsorb(npt)   = xpsorb(iv1)
       casapool%psoilocc(npt)    = xpocc(iv1)
       casapool%Psoillab(npt)    = 0.0_r_2
       casaflux%kmlabp(npt)      = xkmlabp(iso)
       casaflux%psorbmax(npt)    = xpsorbmax(iso)
       casaflux%fpleach(npt)     = xfPleach(iso) /(365.0_r_2)    ! convert from 1/year to 1/day
       ! we used the spatially explicit estimate N fixation by Wang and Houlton (GRL)
       ! casaflux%Nminfix(npt)     = xnfixrate(iv1)/365.0_r_2
       casapool%ratioNCplant(npt,:)  = 1.0_r_2/ratioCNplant(iv1,:)
       casapool%ratioNPplant(npt,:)  = casabiome%ratioNPplantmin(iv1,:)
       casapool%ratioNClitter(npt,:) = casapool%nlitter(npt,:)/(casapool%clitter(npt,:)+1.0e-10_r_2)
       casapool%ratioNPlitter(npt,:) = casapool%nlitter(npt,:)/(casapool%plitter(npt,:)+1.0e-10_r_2)
       casapool%ratioNCsoil(npt,:)   = 1.0_r_2/ratioCNsoil(iv1,:)
       casapool%ratioNPsoil(npt,:)   = ratioNPsoil(iso,:)
       casapool%ratioNCsoilmin(npt,:)   = 1.0_r_2/ratioCNsoilmax(iv1,:)
       casapool%ratioNCsoilmax(npt,:)   = 1.0_r_2/ratioCNsoilmin(iv1,:)
       casapool%ratioNCsoilnew(npt,:)   = casapool%ratioNCsoilmax(npt,:)
    ENDDO

    if(icycle<2) then
       casapool%Nplant(:,:)  = casapool%Cplant(:,:) * casapool%ratioNCplant(:,:)
       casapool%Nsoil(:,:)   = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
    endif
    if(icycle<3) then
       casapool%Psoil(:,:)   = casapool%Nsoil(:,:)/ casapool%ratioNPsoil(:,:)
       casapool%Psoilsorb(:) = casaflux%psorbmax(:) * casapool%psoillab(:) &
            /(casaflux%kmlabp(:)+casapool%psoillab(:))
    endif

  END SUBROUTINE casa_readbiome


  SUBROUTINE casa_readphen(veg, casamet, phen)
    !SUBROUTINE casa_readphen(mvt,veg,casamet,phen)
    ! read in the tabulated modis-derived leaf phenology data
    ! for latitude bands of 79.75 to -55.25
    USE cable_def_types_mod
    USE casadimension
    USE casaparm
    USE casavariable
    USE phenvariable
    IMPLICIT NONE
    !  INTEGER,              INTENT(IN)    :: mvt
    TYPE(veg_parameter_type), INTENT(IN)    :: veg  ! vegetation parameters
    TYPE(casa_met),           INTENT(IN)    :: casamet
    TYPE(phen_variable),      INTENT(INOUT) :: phen

    ! local variables
    INTEGER, PARAMETER             :: npheno=8 ! was 10(IGBP). changed by Q.Zhang @01/12/2011
    INTEGER :: np,nx,ilat
    INTEGER, DIMENSION(271,mvtype) :: greenup, fall,  phendoy1
    INTEGER, DIMENSION(npheno)     :: greenupx,fallx,xphendoy1
    INTEGER, DIMENSION(npheno)     :: ivtx
    REAL(r_2), DIMENSION(271)      :: xlat

    ! initilize for evergreen PFTs
    greenup(:,:) = -50
    fall(:,:)    = 367
    phendoy1(:,:)= 2

    OPEN(101, file=casafile%phen)
    READ(101, *)
    READ(101, *) (ivtx(nx), nx=1, npheno) ! fixed at 10, as only 10 of 17 IGBP PFT
    ! have seasonal leaf phenology
    DO ilat=271,1,-1
       READ(101,*) xlat(ilat), (greenupx(nx), nx=1, npheno), &
            (fallx(nx), nx=1, npheno), (xphendoy1(nx), nx=1, npheno)
       DO nx=1, npheno
          greenup(ilat,ivtx(nx)) = greenupx(nx)
          fall(ilat,ivtx(nx))    = fallx(nx)
          phendoy1(ilat,ivtx(nx))= xphendoy1(nx)
       ENDDO
    ENDDO

    CLOSE(101)

    DO np=1,mp
       ilat = nint((casamet%lat(np)+55.25)/0.5)+1
       ilat = MIN(271,MAX(1,ilat))
       phen%phase(np) = phendoy1(ilat,veg%iveg(np))
       phen%doyphase(np,1) = greenup(ilat,veg%iveg(np)) ! DOY for greenup
       phen%doyphase(np,2) = phen%doyphase(np,1) +14    ! DOY for steady LAI
       phen%doyphase(np,3) = fall(ilat,veg%iveg(np))    ! DOY for leaf senescence
       phen%doyphase(np,4) = phen%doyphase(np,3) +14    ! DOY for minimal LAI season
       IF (phen%doyphase(np,2) > 365) phen%doyphase(np,2)=phen%doyphase(np,2)-365
       IF (phen%doyphase(np,4) > 365) phen%doyphase(np,4)=phen%doyphase(np,4)-365
    ENDDO

  END SUBROUTINE casa_readphen


  SUBROUTINE casa_init(casabiome, casamet, casaflux, casapool, casabal, phen)
    ! mst not used (BP sep2010)
    !! for first time reading file *_1220.csv  (BP may2010)
    !SUBROUTINE casa_init(mst,casapool,casabal,veg)
    !!SUBROUTINE casa_init(mst,casapool,casabal)
    !! end addition (BP may2010)
    !  initialize some values in phenology parameters and leaf growth phase
    USE casadimension
    USE casaparm
    USE casavariable
    USE phenvariable
    ! for first time reading file *_1220.csv  (BP may2010)
    USE cable_def_types_mod
    USE cable_io_vars_module, ONLY: patch
    USE cable_common_module, only: cable_user

    ! end addition (BP may2010)
    IMPLICIT NONE

    TYPE(casa_biome),    INTENT(INOUT) :: casabiome
    TYPE(casa_met),      INTENT(INOUT) :: casamet
    TYPE(casa_flux),     INTENT(INOUT) :: casaflux
    TYPE(casa_pool),     INTENT(INOUT) :: casapool
    TYPE(casa_balance),  INTENT(INOUT) :: casabal
    ! for first time reading file *_1220.csv  (BP may2010)
    TYPE(phen_variable), INTENT(INOUT) :: phen
    ! end addition (BP may2010)

    ! local variables
    INTEGER :: npt

    if (.not. cable_user%casa_fromzero) then
       write(*,*) 'initial pool from restart file'
    endif
    write(*,*) 'icycle, initcasa, mp ', icycle, initcasa, mp
    !phen%phase = 2

    !MC - remove zeroing because done in load_parameters before
    ! !CLN initialise all !!!!! THIS NEEDS FIXING because of e.g. ICE-WATER
    ! casaflux%Cgpp         = 0.0_r_2
    ! casaflux%Cnpp         = 0.0_r_2
    ! casaflux%Crp          = 0.0_r_2
    ! casaflux%Crgplant     = 0.0_r_2
    ! ! casaflux%Nminfix      = 0.0_r_2
    ! casaflux%Nminuptake   = 0.0_r_2
    ! casaflux%Plabuptake   = 0.0_r_2
    ! casaflux%Clabloss     = 0.0_r_2
    ! casaflux%fracClabile  = 0.0_r_2
    ! casaflux%stemnpp      = 0.0_r_2
    ! casaflux%frac_sapwood = 0.0_r_2
    ! casaflux%sapwood_area = 0.0_r_2
    ! casaflux%FluxCtohwp   = 0.0_r_2
    ! casaflux%FluxCtoClear = 0.0_r_2
    ! casaflux%fracCalloc   = 0.0_r_2
    ! casaflux%fracNalloc   = 0.0_r_2
    ! casaflux%fracPalloc   = 0.0_r_2
    ! casaflux%Crmplant     = 0.0_r_2
    ! casaflux%kplant       = 0.0_r_2
    ! casaflux%kplant_fire  = 0.0_r_2
    ! casaflux%kplant_tot   = 0.0_r_2

    ! casaflux%fromPtoL      = 0.0_r_2
    ! casaflux%fromPtoL_fire = 0.0_r_2

    ! casaflux%Cnep         = 0.0_r_2
    ! casaflux%Crsoil       = 0.0_r_2
    ! casapool%dClabiledt   = 0.0_r_2
    ! ! casaflux%Nmindep      =  casaflux%Nmindep /2.0_r_2
    ! ! casaflux%Nmindep      = 0.0_r_2
    ! casaflux%Nminloss     = 0.0_r_2
    ! casaflux%Nminleach    = 0.0_r_2
    ! casaflux%Nupland      = 0.0_r_2
    ! casaflux%Nlittermin   = 0.0_r_2
    ! casaflux%Nsmin        = 0.0_r_2
    ! casaflux%Nsimm        = 0.0_r_2
    ! casaflux%Nsnet        = 0.0_r_2
    ! ! casaflux%fNminloss    = 0.0_r_2
    ! ! casaflux%fNminleach   = 0.0_r_2
    ! ! casaflux%Pdep         = 0.0_r_2
    ! ! casaflux%Pwea         = 0.0_r_2
    ! casaflux%Pleach       = 0.0_r_2
    ! casaflux%Ploss        = 0.0_r_2
    ! casaflux%Pupland      = 0.0_r_2
    ! casaflux%Plittermin   = 0.0_r_2
    ! casaflux%Psmin        = 0.0_r_2
    ! casaflux%Psimm        = 0.0_r_2
    ! casaflux%Psnet        = 0.0_r_2
    ! ! casaflux%fPleach      = 0.0_r_2 !vh ! this should be a parameter, not a flux variable
    ! casaflux%kplab        = 0.0_r_2
    ! casaflux%kpsorb       = 0.0_r_2
    ! casaflux%kpocc        = 0.0_r_2
    ! ! casaflux%kmlabp       = 0.0_r_2  !vh ! this should be a paramter, not a flux variable
    ! ! casaflux%Psorbmax     = 0.0_r_2 !vh ! this should be a paramter, not a flux variable

    ! casaflux%klitter       = 0.0_r_2
    ! casaflux%klitter_fire  = 0.0_r_2
    ! casaflux%ksoil         = 0.0_r_2
    ! casaflux%fromLtoS      = 0.0_r_2
    ! casaflux%fromStoS      = 0.0_r_2
    ! casaflux%fromLtoCO2    = 0.0_r_2
    ! casaflux%fromStoCO2    = 0.0_r_2
    ! casaflux%FluxCtolitter = 0.0_r_2
    ! casaflux%FluxNtolitter = 0.0_r_2
    ! casaflux%FluxPtolitter = 0.0_r_2
    ! casaflux%FluxCtosoil   = 0.0_r_2
    ! casaflux%FluxNtosoil   = 0.0_r_2
    ! casaflux%FluxPtosoil   = 0.0_r_2
    ! casaflux%FluxCtoCO2    = 0.0_r_2
    ! casaflux%FluxCtoCO2_plant_fire  = 0.0_r_2
    ! casaflux%FluxCtoCO2_litter_fire = 0.0_r_2

    ! casaflux%FluxFromPtoL       = 0.0_r_2
    ! casaflux%FluxFromLtoS       = 0.0_r_2
    ! casaflux%FluxFromStoS       = 0.0_r_2
    ! casaflux%FluxFromPtoCO2     = 0.0_r_2
    ! casaflux%FluxFromLtoCO2     = 0.0_r_2
    ! casaflux%FluxFromStoCO2     = 0.0_r_2
    ! casaflux%FluxFromPtoHarvest = 0.0_r_2

    ! casaflux%Cplant_turnover                     = 0.0_r_2
    ! casaflux%Cplant_turnover_disturbance         = 0.0_r_2
    ! casaflux%Cplant_turnover_crowding            = 0.0_r_2
    ! casaflux%Cplant_turnover_resource_limitation = 0.0_r_2
    ! casaflux%fHarvest = 0.0_r_2
    ! casaflux%Charvest = 0.0_r_2
    ! casaflux%Nharvest = 0.0_r_2
    ! casaflux%fcrop    = 0.0_r_2

    phen%doyphase(:,1) = -50
    phen%doyphase(:,2) = phen%doyphase(:,1) + 14
    phen%doyphase(:,3) = 367
    phen%doyphase(:,4) = phen%doyphase(:,3) + 14
    phen%phase(:) = 2
    phen%phen(:)  = 1
    phen%aphen(:) = 0

    ! casapool%dCplantdt     = 0.0_r_2
    ! casapool%dNplantdt     = 0.0_r_2
    ! casapool%dPplantdt     = 0.0_r_2
    ! casapool%dNsoilmindt   = 0.0_r_2
    ! casapool%dPsoillabdt   = 0.0_r_2
    ! casapool%dPsoilsorbdt  = 0.0_r_2
    ! casapool%dPsoiloccdt   = 0.0_r_2
    ! casapool%dClitterdt    = 0.0_r_2
    ! casapool%dNlitterdt    = 0.0_r_2
    ! casapool%dPlitterdt    = 0.0_r_2
    ! casapool%dCsoildt      = 0.0_r_2
    ! casapool%dNsoildt      = 0.0_r_2
    ! casapool%dPsoildt      = 0.0_r_2
    ! casapool%ratioNClitter = 0.0_r_2
    ! casapool%ratioNCsoil   = 0.0_r_2
    ! casapool%ratioPClitter = 0.0_r_2
    ! casapool%ratioPCsoil   = 0.0_r_2
    ! !CLN add more if necessary

    casapool%Nsoilmin(:) = 2.5_r_2
    !MC - Check if next commented 2 if-blocks are necessary
    ! if (icycle > 1) then
    !    casapool%nplant     = MAX(1.e-6_r_2, casapool%nplant)
    !    casapool%nlitter    = MAX(1.e-6_r_2, casapool%nlitter)
    !    casapool%nsoil      = MAX(1.e-6_r_2, casapool%nsoil)
    !    casapool%nsoilmin   = MAX(1.e-6_r_2, casapool%nsoilmin)
    !    casabal%cplantlast  = casapool%cplant
    !    casabal%clitterlast = casapool%clitter
    !    casabal%csoillast   = casapool%csoil
    !    casabal%clabilelast = casapool%clabile
    ! endif
    ! if (icycle > 2) then
    !    casapool%pplant       = MAX(1.e-7_r_2, casapool%pplant)
    !    casapool%plitter      = MAX(1.e-7_r_2, casapool%plitter)
    !    casapool%psoil        = MAX(1.e-7_r_2, casapool%psoil)
    !    casapool%Psoillab     = MAX(1.e-7_r_2, casapool%psoillab)
    !    casapool%psoilsorb    = MAX(1.e-7_r_2, casapool%psoilsorb)
    !    casapool%psoilocc     = MAX(1.e-7_r_2, casapool%psoilocc)
    !    casabal%pplantlast    = casapool%pplant
    !    casabal%plitterlast   = casapool%plitter
    !    casabal%psoillast     = casapool%psoil
    !    casabal%psoillablast  = casapool%psoillab
    !    casabal%psoilsorblast = casapool%psoilsorb
    !    casabal%psoilocclast  = casapool%psoilocc
    ! endif
    if ((initcasa == 1) .and. cable_user%casa_fromzero) then
       write(*,*) 'casa_init: not using restart file!'
       write(*,*) 'Using input from readbiome.!!!'
       write(*,*) 'initialising frac_sapwood=1 and sapwood_area = 0)'
       casaflux%frac_sapwood(:) = 1.0_r_2
       casaflux%sapwood_area(:) = 0.0_r_2
#ifndef UM_BUILD
    else if ((initcasa == 1) .and. (.not. cable_user%casa_fromzero)) then
       call read_casa_restart_nc(casabiome, casamet, casapool, casaflux, &
            casabal, phen)
       call zero_casa_dump(casamet, phen)
#endif
    else
       do npt=1, mp
          casamet%lon(npt) = real(patch(npt)%longitude, r_2)
          casamet%lat(npt) = real(patch(npt)%latitude, r_2)
       enddo
    endif

    !MC - remove zeroing because done in load_parameters before
    ! where(casamet%lnonwood==1) casapool%cplant(:,WOOD) = 0.0_r_2
    ! where(casamet%lnonwood==1) casapool%nplant(:,WOOD) = 0.0_r_2
    ! where(casamet%lnonwood==1) casapool%pplant(:,WOOD) = 0.0_r_2

    ! ! reset labile C pool, comment out by Q.Zhang 10/09/2011: casapool%clabile = 0.0
    ! ! check pool sizes
    ! casapool%Ctot_0     = 0.0_r_2
    ! casapool%Ctot       = 0.0_r_2
    ! casapool%cplant     = MAX(0.0_r_2, casapool%cplant)
    ! casapool%clitter    = MAX(0.0_r_2, casapool%clitter)
    ! casapool%csoil      = MAX(0.0_r_2, casapool%csoil)
    ! casabal%cplantlast  = casapool%cplant
    ! casabal%clitterlast = casapool%clitter
    ! casabal%csoillast   = casapool%csoil
    ! casabal%clabilelast = casapool%clabile
    ! casabal%sumcbal     = 0.0_r_2
    ! casabal%FCgppyear   = 0.0_r_2
    ! casabal%FCrpyear    = 0.0_r_2
    ! casabal%FCnppyear   = 0.0_r_2
    ! casabal%FCrsyear    = 0.0_r_2
    ! casabal%FCneeyear   = 0.0_r_2
    ! ! vh !
    ! WHERE (casamet%lnonwood==1) casapool%cplant(:,WOOD) = 0.0_r_2
    ! IF (icycle == 1) THEN
    !    casapool%Nplant(:,:) = casapool%cplant(:,:) * casapool%ratioNCplant(:,:)
    !    casapool%Nsoil(:,:)  = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
    !    casapool%Psoil(:,:)  = casapool%Nsoil(:,:) / casapool%ratioNPsoil(:,:)
    !    casapool%Nsoilmin(:) = 2.5_r_2
    ! ENDIF

    ! IF (icycle >= 1) THEN
    !    casapool%nplant     = MAX(1.0e-6_r_2, casapool%nplant)
    !    casapool%nlitter    = MAX(1.0e-6_r_2, casapool%nlitter)
    !    casapool%nsoil      = MAX(1.0e-6_r_2, casapool%nsoil)
    !    casapool%nsoilmin   = MAX(1.0e-6_r_2, casapool%nsoilmin)
    !    casabal%nplantlast  = casapool%nplant
    !    casabal%nlitterlast = casapool%nlitter
    !    casabal%nsoillast   = casapool%nsoil
    !    casabal%nsoilminlast= casapool%nsoilmin
    !    casabal%sumnbal     = 0.0_r_2
    !    casabal%FNdepyear   = 0.0_r_2
    !    casabal%FNfixyear   = 0.0_r_2
    !    casabal%FNsnetyear  = 0.0_r_2
    !    casabal%FNupyear    = 0.0_r_2
    !    casabal%FNleachyear = 0.0_r_2
    !    casabal%FNlossyear  = 0.0_r_2
    !    ! vh !
    !    WHERE(casamet%lnonwood==1) casapool%nplant(:,WOOD) = 0.0_r_2
    ! ENDIF

    ! IF (icycle >=1 ) THEN
    !    casapool%pplant       = MAX(1.0e-7_r_2, casapool%pplant)
    !    casapool%plitter      = MAX(1.0e-7_r_2, casapool%plitter)
    !    casapool%psoil        = MAX(1.0e-7_r_2, casapool%psoil)
    !    casapool%Psoillab     = MAX(1.0e-7_r_2, casapool%psoillab)  ! was 2.0, changed according to  YP
    !    casapool%psoilsorb    = MAX(1.0e-7_r_2, casapool%psoilsorb) ! was 10.0, -
    !    casapool%psoilocc     = MAX(1.0e-7_r_2, casapool%psoilocc)  ! was 50.0, -
    !    casabal%pplantlast    = casapool%pplant
    !    casabal%plitterlast   = casapool%plitter
    !    casabal%psoillast     = casapool%psoil
    !    casabal%psoillablast  = casapool%psoillab
    !    casabal%psoilsorblast = casapool%psoilsorb
    !    casabal%psoilocclast  = casapool%psoilocc
    !    casabal%sumpbal       = 0.0_r_2
    !    casabal%FPweayear     = 0.0_r_2
    !    casabal%FPdustyear    = 0.0_r_2
    !    casabal%FPsnetyear    = 0.0_r_2
    !    casabal%FPupyear      = 0.0_r_2
    !    casabal%FPleachyear   = 0.0_r_2
    !    casabal%FPlossyear    = 0.0_r_2
    !    ! vh !
    !    WHERE(casamet%lnonwood==1) casapool%pplant(:,WOOD) = 0.0_r_2
    ! ENDIF

  END SUBROUTINE casa_init


  SUBROUTINE casa_poolout(ktau, veg, soil, casabiome, casapool, casaflux, casamet, &
       casabal, phen)

    USE cable_def_types_mod
    USE casadimension
    USE casaparm
    USE casavariable
    USE phenvariable
    USE cable_common_module, only: cable_user

    IMPLICIT NONE

    INTEGER,                   INTENT(IN)    :: ktau
    TYPE(veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    TYPE(casa_biome),          INTENT(INOUT) :: casabiome
    TYPE(casa_pool),           INTENT(INOUT) :: casapool
    TYPE(casa_flux),           INTENT(INOUT) :: casaflux
    TYPE(casa_met),            INTENT(INOUT) :: casamet
    TYPE(casa_balance),        INTENT(INOUT) :: casabal
    TYPE(phen_variable),       INTENT(INOUT) :: phen

    ! local variables
    REAL(r_2), DIMENSION(mso) :: Psorder, xPsoil50 ! , Pweasoil
    REAL(r_2), DIMENSION(mso) :: fracPlab, fracPocc ! , fracPsorb, fracPorg
    REAL(r_2), DIMENSION(mp)  :: totPsoil
    INTEGER :: npt, nout, nso

    ! Soiltype     soilnumber soil P(g P/m2)
    ! Alfisol     1       61.3
    ! Andisol     2       103.9
    ! Aridisol    3       92.8
    ! Entisol     4       136.9
    ! Gellisol    5       98.2
    ! Histosol    6       107.6
    ! Inceptisol  7       84.1
    ! Mollisol    8       110.1
    ! Oxisol      9       35.4
    ! Spodosol    10      41.0
    ! Ultisol     11      51.5
    ! Vertisol    12      190.6
    DATA Psorder/61.3_r_2, 103.9_r_2, 92.8_r_2, 136.9_r_2, 98.2_r_2, 107.6_r_2, 84.1_r_2, &
         110.1_r_2, 35.4_r_2, 41.0_r_2, 51.5_r_2, 190.6_r_2/
    ! DATA Pweasoil/0.05_r_2, 0.04_r_2, 0.03_r_2, 0.02_r_2, 0.01_r_2, 0.009_r_2, 0.008_r_2, &
    !      0.007_r_2, 0.006_r_2, 0.005_r_2, 0.004_r_2, 0.003_r_2/
    DATA fracPlab/0.08_r_2, 0.08_r_2, 0.10_r_2, 0.02_r_2, 0.08_r_2, 0.08_r_2, 0.08_r_2, &
         0.06_r_2, 0.02_r_2, 0.05_r_2, 0.09_r_2, 0.05_r_2/
    ! DATA fracPsorb/0.32_r_2, 0.37_r_2, 0.57_r_2, 0.67_r_2, 0.37_r_2, 0.37_r_2, 0.37_r_2, &
    !      0.32_r_2, 0.24_r_2, 0.22_r_2, 0.21_r_2, 0.38_r_2/
    DATA fracPocc/0.36_r_2, 0.38_r_2, 0.25_r_2, 0.26_r_2, 0.38_r_2, 0.38_r_2, 0.38_r_2, &
         0.44_r_2, 0.38_r_2, 0.38_r_2, 0.37_r_2, 0.45_r_2/
    ! DATA fracPorg/0.25_r_2, 0.17_r_2, 0.08_r_2, 0.05_r_2, 0.17_r_2, 0.17_r_2, 0.17_r_2, &
    !      0.18_r_2, 0.36_r_2, 0.35_r_2, 0.34_r_2, 0.12_r_2/
    DATA xPsoil50/7.6_r_2, 4.1_r_2, 4.2_r_2, 3.4_r_2, 4.1_r_2, 4.1_r_2, 4.8_r_2, 4.1_r_2, &
         6.9_r_2, 6.9_r_2, 6.9_r_2, 1.7_r_2/
    !
    ! estimated based on Yang, Post and Jain (2013)
    !   Soiltype     soilnumber soil P(g P/m2  top 50 cm)
    !   Alfisol     1       400
    !   Andisol     2       426
    !   Aridisol    3       352
    !   Entisol     4       490
    !   Gellisol    5       403
    !   Histosol    6       441
    !   Inceptisol  7       501
    !   Mollisol    8       358
    !   Oxisol      9       96
    !   Spodosol    10      364
    !   Ultisol     11      272
    !   Vertisol    12      430
    !  DATA psorder/400.0,426.0,352.0,490.0,403.0,441.0,501.0,358.0,96.0,364.0,272.0,430.0/
    !  DATA pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/
    !  DATA fracpLab/0.07,0.04,0.08,0.10,0.08,0.10,0.12,0.05,0.05,0.06,0.06,0.05/
    !  DATA fracPsorb/0.30,0.44,0.69,0.53,0.37,0.14,0.24,0.32,0.15,0.21,0.17,0.35/
    !  DATA fracPocc/0.38,0.22,0.18,0.22,0.38,0.42,0.23,0.44,0.60,0.30,0.51,0.48/
    !  DATA fracPorg/0.25,0.30,0.05,0.15,0.17,0.34,0.41,0.19,0.20,0.43,0.26,0.12/
    !  DATA xpsoil50/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/

    write(*,*) 'Within casa_poolout, mp = ', mp
    nout=103
    OPEN(nout,file=casafile%cnpepool)
    write(*,*) 'Opened file ', casafile%cnpepool

    casabal%sumcbal = MIN(9999.0_r_2, MAX(-9999.0_r_2, casabal%sumcbal))
    casabal%sumnbal = MIN(9999.0_r_2, MAX(-9999.0_r_2, casabal%sumnbal))
    casabal%sumpbal = MIN(9999.0_r_2, MAX(-9999.0_r_2, casabal%sumpbal))

    DO npt =1, mp
       nso = casamet%isorder(npt)
       totpsoil(npt) = psorder(nso) *xpsoil50(nso)
       if (casamet%iveg2(npt)>0) then
          IF (icycle<2) THEN
             casapool%Nplant(npt,:) = casapool%ratioNCplant(npt,:)  &
                  * casapool%cplant(npt,:)
             casapool%Nlitter(npt,:)= casapool%ratioNClitter(npt,:) &
                  * casapool%clitter(npt,:)
             casapool%Nsoil(npt,:)  = casapool%ratioNCsoil(npt,:)   &
                  * casapool%Csoil(npt,:)
             casapool%nsoilmin(npt) = 2.0_r_2
             casabal%sumnbal(npt)   = 0.0_r_2
             if(casamet%iveg2(npt)==grass) then
                casapool%nplant(npt,wood) = 0.0_r_2
                casapool%nlitter(npt,cwd) = 0.0_r_2
             endif
          ENDIF

          IF (icycle<3) THEN
             casabal%sumpbal(npt)    = 0.0_r_2
             casapool%pplant(npt,:)  = casapool%Nplant(npt,:)/casapool%ratioNPplant(npt,:)
             casapool%plitter(npt,:) = casapool%Nlitter(npt,:)/(casapool%ratioNPlitter(npt,:) + 1.e-10_r_2)
             casapool%psoil(npt,:)   = casapool%Nsoil(npt,:)/casapool%ratioNPsoil(npt,:)
             casapool%psoillab(npt)  = totpsoil(npt) *fracpLab(nso)
             casapool%psoilsorb(npt) = casaflux%psorbmax(npt) * casapool%psoillab(npt) &
                  / (casaflux%kmlabp(npt)+casapool%psoillab(npt))
             casapool%psoilocc(npt) = totpsoil(npt) *fracPocc(nso)
             if(casamet%iveg2(npt)==grass) then
                casapool%pplant(npt,wood) = 0.0_r_2
                casapool%plitter(npt,cwd) = 0.0_r_2
             endif
          ENDIF
       else
          casapool%cplant(npt,:)  = 0.0_r_2
          casapool%clitter(npt,:) = 0.0_r_2
          casapool%csoil(npt,:)   = 0.0_r_2
          casapool%clabile(npt)   = 0.0_r_2
          casapool%nplant(npt,:)  = 0.0_r_2
          casapool%nlitter(npt,:) = 0.0_r_2
          casapool%nsoil(npt,:)   = 0.0_r_2
          casapool%nsoilmin(npt)  = 0.0_r_2
          casapool%pplant(npt,:)  = 0.0_r_2
          casapool%plitter(npt,:) = 0.0_r_2
          casapool%psoil(npt,:)   = 0.0_r_2
          casapool%psoillab(npt)  = 0.0_r_2
          casapool%psoilsorb(npt) = 0.0_r_2
          casapool%psoilocc(npt)  = 0.0_r_2
          casabal%sumcbal(npt)    = 0.0_r_2
          casabal%sumnbal(npt)    = 0.0_r_2
          casabal%sumpbal(npt)    = 0.0_r_2
       endif

       !! vh_js  !!
       IF (cable_user%CALL_POP) THEN
          WRITE(nout,92) ktau,npt,veg%iveg(npt),soil%isoilm(npt) ,     &
               casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
               casamet%areacell(npt)*(1.0e-9_r_2),casamet%glai(npt),       &
               casabiome%sla(veg%iveg(npt)), phen%phase(npt), &
               phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
               casapool%clabile(npt), &
               casapool%cplant(npt,:),casapool%clitter(npt,:),casapool%csoil(npt,:), &
               casaflux%frac_sapwood(npt), casaflux%sapwood_area(npt), &
               casapool%nplant(npt,:),casapool%nlitter(npt,:),casapool%nsoil(npt,:), &
               casapool%nsoilmin(npt),casapool%pplant(npt,:),          &
               casapool%plitter(npt,:), casapool%psoil(npt,:),         &
               casapool%psoillab(npt),casapool%psoilsorb(npt),casapool%psoilocc(npt), &
               casabal%sumcbal(npt),casabal%sumnbal(npt),casabal%sumpbal(npt)
       ELSE
          WRITE(nout,92) ktau,npt,veg%iveg(npt),soil%isoilm(npt),     &
               casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
               casamet%areacell(npt)*(1.0e-9_r_2),casamet%glai(npt),       &
               casabiome%sla(veg%iveg(npt)), phen%phase(npt), &
               phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
               casapool%clabile(npt), &
               casapool%cplant(npt,:),casapool%clitter(npt,:),casapool%csoil(npt,:), &
               casapool%nplant(npt,:),casapool%nlitter(npt,:),casapool%nsoil(npt,:), &
               casapool%nsoilmin(npt),casapool%pplant(npt,:),          &
               casapool%plitter(npt,:), casapool%psoil(npt,:),         &
               casapool%psoillab(npt),casapool%psoilsorb(npt),casapool%psoilocc(npt), &
               casabal%sumcbal(npt),casabal%sumnbal(npt),casabal%sumpbal(npt)
       ENDIF

    ENDDO

    CLOSE(nout)

92  format(5(i6,',',2x),5(f18.6,',',2x),2(i6,',',2x),100(f18.6,',',2x))

  END SUBROUTINE casa_poolout


  ! casa_fluxout output data for Julie Tang; comment out (BP apr2010)
  SUBROUTINE casa_fluxout(myear,veg,soil,casabal,casamet)
    !SUBROUTINE casa_fluxout(myear,clitterinput,csoilinput)

    USE cable_def_types_mod
    ! USE cableDeclare, ONLY: veg, soil
    USE casadimension
    USE casaparm
    USE casavariable
    USE phenvariable
    !  USE casaDeclare

    IMPLICIT NONE

    INTEGER,                    INTENT(IN)    :: myear
    TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    TYPE (casa_met),            INTENT(INOUT) :: casamet
    TYPE (casa_balance),        INTENT(INOUT) :: casabal
    !  REAL(r_2),    INTENT(IN) :: clitterinput(mp,3),csoilinput(mp,3)

    ! local variables
    INTEGER ::  npt,nout
    REAL(r_2) :: xyear, totGPP, totNPP

    totGPP = 0.0_r_2
    totNPP = 0.0_r_2
    nout   = 104
    xyear  = 1.0_r_2/real(myear,r_2)
    casabal%FCgppyear    = casabal%FCgppyear    * xyear
    casabal%FCnppyear    = casabal%FCnppyear    * xyear
    casabal%FCrmleafyear = casabal%FCrmleafyear * xyear
    casabal%FCrmwoodyear = casabal%FCrmwoodyear * xyear
    casabal%FCrmrootyear = casabal%FCrmrootyear * xyear
    casabal%FCrgrowyear  = casabal%FCrgrowyear  * xyear
    casabal%FCrsyear     = casabal%FCrsyear     * xyear
    casabal%FCneeyear    = casabal%FCneeyear    * xyear
    casabal%FNdepyear    = casabal%FNdepyear    * xyear
    casabal%FNfixyear    = casabal%FNfixyear    * xyear
    casabal%FNsnetyear   = casabal%FNsnetyear   * xyear
    casabal%FNupyear     = casabal%FNupyear     * xyear
    casabal%FNleachyear  = casabal%FNleachyear  * xyear
    casabal%FNlossyear   = casabal%FNlossyear   * xyear
    casabal%FPweayear    = casabal%FPweayear    * xyear
    casabal%FPdustyear   = casabal%FPdustyear   * xyear
    casabal%FPsnetyear   = casabal%FPsnetyear   * xyear
    casabal%FPupyear     = casabal%FPupyear     * xyear
    casabal%FPleachyear  = casabal%FPleachyear  * xyear
    casabal%FPlossyear   = casabal%FPlossyear   * xyear
    !  clitterinput = clitterinput * xyear
    !  csoilinput   = csoilinput   * xyear

    write(*,*) 'writing CNP fluxes out to file ', casafile%cnpflux
    OPEN(nout,file=casafile%cnpflux)
    DO npt =1,mp
       SELECT CASE(icycle)
       CASE(1)
          WRITE(nout,*) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
               casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
               casamet%areacell(npt)*(1.0e-9_r_2),casabal%Fcgppyear(npt),  &
               casabal%Fcnppyear(npt),  &
               casabal%Fcrmleafyear(npt),casabal%Fcrmwoodyear(npt),     &
               casabal%Fcrmrootyear(npt),casabal%Fcrgrowyear(npt),     &
               casabal%Fcrsyear(npt),casabal%Fcneeyear(npt) ! , &
          ! clitterinput(npt,:),csoilinput(npt,:)
       CASE(2)
          WRITE(nout,*) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
               casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
                                !MC - casamet%areacell is different between MPI and serial
               casamet%areacell(npt)*(1.0e-9_r_2),casabal%Fcgppyear(npt),  &
               casabal%FCnppyear(npt),                                 &
               casabal%Fcrmleafyear(npt),casabal%Fcrmwoodyear(npt),     &
               casabal%Fcrmrootyear(npt),casabal%Fcrgrowyear(npt),     &
               casabal%FCrsyear(npt), casabal%FCneeyear(npt),          &
               ! clitterinput(npt,:),csoilinput(npt,:), &
               casabal%FNdepyear(npt),casabal%FNfixyear(npt),casabal%FNsnetyear(npt), &
               casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt)
       CASE(3)
          WRITE(nout,*) myear,npt,veg%iveg(npt),soil%isoilm(npt), &
               casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt),  &
               casamet%areacell(npt)*(1.0e-9_r_2),casabal%Fcgppyear(npt), &
               casabal%FCnppyear(npt),                                  &
               casabal%Fcrmleafyear(npt),casabal%Fcrmwoodyear(npt),     &
               casabal%Fcrmrootyear(npt),casabal%Fcrgrowyear(npt),     &
               casabal%FCrsyear(npt),   casabal%FCneeyear(npt),         &
                                !        clitterinput(npt,:),csoilinput(npt,:), &
               casabal%FNdepyear(npt),casabal%FNfixyear(npt),  casabal%FNsnetyear(npt),&
               casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt),&
               casabal%FPweayear(npt),casabal%FPdustyear(npt), casabal%FPsnetyear(npt),&
               casabal%FPupyear(npt), casabal%FPleachyear(npt),casabal%FPlossyear(npt)
       END SELECT
       totGPP = totGPP+casabal%Fcgppyear(npt)* casamet%areacell(npt)
       totNPP = totNPP+casabal%Fcnppyear(npt)* casamet%areacell(npt)
    ENDDO

    write(*,*) 'totGPP global = ', totGPP * (1.0e-15_r_2)
    write(*,*) 'totNPP global = ', totNPP * (1.0e-15_r_2)

    CLOSE(nout)

  END SUBROUTINE casa_fluxout


  ! clitterinput and csoilinput are for Julie Tang; comment out (BP apr2010)
  !SUBROUTINE casa_cnpflux(clitterinput,csoilinput)
  SUBROUTINE casa_cnpflux(casaflux,casapool,casabal,zeroflux)

    USE cable_def_types_mod
    USE casadimension
    USE casaparm
    USE casavariable

    IMPLICIT NONE

    TYPE(casa_flux),    INTENT(INOUT) :: casaflux
    TYPE(casa_pool),    INTENT(IN)    :: casapool
    TYPE(casa_balance), INTENT(INOUT) :: casabal
    LOGICAL,            intent(in)    :: zeroflux

    IF (zeroflux) THEN
       casabal%FCgppyear    = 0.0_r_2
       casabal%FCrpyear     = 0.0_r_2
       casabal%FCrmleafyear = 0.0_r_2
       casabal%FCrmwoodyear = 0.0_r_2
       casabal%FCrmrootyear = 0.0_r_2
       casabal%FCrgrowyear  = 0.0_r_2
       casabal%FCnppyear    = 0.0_r_2
       casabal%FCrsyear     = 0.0_r_2
       casabal%FCneeyear    = 0.0_r_2
       casabal%dCdtyear     = 0.0_r_2

       casabal%FNdepyear   = 0.0_r_2
       casabal%FNfixyear   = 0.0_r_2
       casabal%FNsnetyear  = 0.0_r_2
       casabal%FNupyear    = 0.0_r_2
       casabal%FNleachyear = 0.0_r_2
       casabal%FNlossyear  = 0.0_r_2

       casabal%FPweayear   = 0.0_r_2
       casabal%FPdustyear  = 0.0_r_2
       casabal%FPsnetyear  = 0.0_r_2
       casabal%FPupyear    = 0.0_r_2
       casabal%FPleachyear = 0.0_r_2
       casabal%FPlossyear  = 0.0_r_2

       casaflux%FluxCtohwp   = 0.0_r_2
       casaflux%FluxNtohwp   = 0.0_r_2
       casaflux%FluxPtohwp   = 0.0_r_2
       casaflux%FluxCtoclear = 0.0_r_2
       casaflux%FluxNtoclear = 0.0_r_2
       casaflux%FluxPtoclear = 0.0_r_2
       casaflux%CtransferLUC = 0.0_r_2

       !MC - commented because otherwise serial and MPI versions different
       ! casaflux%Cplant_turnover_disturbance         = 0.0_r_2
       ! casaflux%Cplant_turnover_crowding            = 0.0_r_2
       ! casaflux%Cplant_turnover_resource_limitation = 0.0_r_2

    ELSE

       casabal%FCgppyear       = casabal%FCgppyear       + casaflux%Cgpp              * deltpool
       casabal%FCrpyear        = casabal%FCrpyear        + casaflux%Crp               * deltpool
       casabal%FCrmleafyear(:) = casabal%FCrmleafyear(:) + casaflux%Crmplant(:,leaf)  * deltpool
       casabal%FCrmwoodyear(:) = casabal%FCrmwoodyear(:) + casaflux%Crmplant(:,wood)  * deltpool
       casabal%FCrmrootyear(:) = casabal%FCrmrootyear(:) + casaflux%Crmplant(:,froot) * deltpool
       casabal%FCrgrowyear     = casabal%FCrgrowyear     + casaflux%Crgplant          * deltpool
       ! change made ypwang 17-nov-2013 to accoutn for change in labile carbon pool  size
       casabal%FCnppyear = casabal%FCnppyear + (casaflux%Cnpp+casapool%dClabiledt) * deltpool
       casabal%FCrsyear  = casabal%FCrsyear  + casaflux%Crsoil                     * deltpool
       !casabal%FCneeyear = casabal%FCneeyear &
       !     + (casaflux%Cnpp+casapool%dClabiledt-casaflux%Crsoil) * deltpool
       casabal%FCneeyear = casabal%FCneeyear + (casaflux%Cnpp-casaflux%Crsoil) * deltpool
       casabal%dCdtyear  = casabal%dCdtyear  + (casapool%Ctot-casapool%Ctot_0) * deltpool

       !  DO n=1,3
       !    clitterinput(:,n)= clitterinput(:,n) + casaflux%kplant(:,n) * casapool%cplant(:,n) * deltpool
       !    csoilinput(:,n) = csoilinput(:,n) + casaflux%fluxCtosoil(:,n) * deltpool
       !    !csoilinput(:,n) = csoilinput(:,n)+casaflux%fluxCtolitter(:,n)*deltpool
       !  ENDDO

       IF (icycle > 1) THEN
          casabal%FNdepyear   = casabal%FNdepyear   + casaflux%Nmindep    * deltpool
          casabal%FNfixyear   = casabal%FNfixyear   + casaflux%Nminfix    * deltpool
          casabal%FNsnetyear  = casabal%FNsnetyear  + casaflux%Nsnet      * deltpool
          casabal%FNupyear    = casabal%FNupyear    + casaflux%Nminuptake * deltpool
          casabal%FNleachyear = casabal%FNleachyear + casaflux%Nminleach  * deltpool
          casabal%FNlossyear  = casabal%FNlossyear  + casaflux%Nminloss   * deltpool
       ENDIF

       IF (icycle > 2) THEN
          casabal%FPweayear   = casabal%FPweayear   + casaflux%Pwea       * deltpool
          casabal%FPdustyear  = casabal%FPdustyear  + casaflux%Pdep       * deltpool
          casabal%FPsnetyear  = casabal%FPsnetyear  + casaflux%Psnet      * deltpool
          casabal%FPupyear    = casabal%FPupyear    + casaflux%Plabuptake * deltpool
          casabal%FPleachyear = casabal%FPleachyear + casaflux%Pleach     * deltpool
          casabal%FPlossyear  = casabal%FPlossyear  + casaflux%Ploss      * deltpool
       ENDIF

    ENDIF

  END SUBROUTINE casa_cnpflux


  ! changed by yp wang following Chris Lu 5/nov/2012
  SUBROUTINE biogeochem(idoY,LALLOC,veg,soil,casabiome,casapool,casaflux, &
       casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf,xkleafcold,xkleafdry, &
       cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd, &
       nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd, &
       pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

    USE cable_def_types_mod
    USE casadimension
    USE casa_cnp_module
    USE POP_TYPES,            ONLY: POP_TYPE

    IMPLICIT NONE

    INTEGER,                   INTENT(IN)    :: idoy
    INTEGER,                   INTENT(IN)    :: LALLOC
    TYPE(veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    TYPE(casa_biome),          INTENT(INOUT) :: casabiome
    TYPE(casa_pool),           INTENT(INOUT) :: casapool
    TYPE(casa_flux),           INTENT(INOUT) :: casaflux
    TYPE(casa_met),            INTENT(INOUT) :: casamet
    TYPE(casa_balance),        INTENT(INOUT) :: casabal
    TYPE(phen_variable),       INTENT(INOUT) :: phen
    TYPE(POP_TYPE),            INTENT(IN)    :: POP
    TYPE(climate_TYPE),        INTENT(IN)    :: climate
    real(r_2), dimension(mp),  intent(out)   :: xnplimit
    real(r_2), dimension(mp),  intent(out)   :: xkNlimiting
    real(r_2), dimension(mp),  intent(out)   :: xklitter, xksoil
    real(r_2), dimension(mp),  intent(out)   :: xkleafcold, xkleafdry, xkleaf
    ! added by ypwang following Chris Lu 5/nov/2012
    real(r_2), dimension(mp),  intent(out)   :: &
         cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd, &
         nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd, &
         pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd

    ! local variables
    REAL(r_2), DIMENSION(mp) :: xNPuptake
    INTEGER :: j
    REAL(r_2), ALLOCATABLE :: tmp(:)

    xKNlimiting = 1.0_r_2

    ! zero annual sums
    if (idoy==1) CALL casa_cnpflux(casaflux, casapool, casabal, .true.)

    IF (cable_user%PHENOLOGY_SWITCH.eq.'MODIS') THEN
       call phenology(idoy,veg,phen)
    ENDIF
    call avgsoil(veg,soil,casamet)
    call casa_rplant(veg,casabiome,casapool,casaflux,casamet,climate)

    IF (.NOT.cable_user%CALL_POP) THEN
       call casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen,LALLOC)
    ENDIF

    call casa_xrateplant(xkleafcold, xkleafdry, xkleaf, veg, casabiome, &
         casamet, phen)
    call casa_coeffplant(xkleafcold, xkleafdry, xkleaf, veg, casabiome, casapool, &
         casaflux, casamet)

    call casa_xnp(xnplimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)

    IF (cable_user%CALL_POP) THEN

       call casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen,LALLOC)
       WHERE (pop%pop_grid(:)%cmass_sum_old.gt.0.001_r_2 .and. pop%pop_grid(:)%cmass_sum.gt.0.001_r_2 )
          casaflux%frac_sapwood(POP%Iwood) = POP%pop_grid(:)%csapwood_sum / POP%pop_grid(:)%cmass_sum
          casaflux%sapwood_area(POP%Iwood) = max(POP%pop_grid(:)%sapwood_area/10000._r_2, 1.0e-6_r_2)
          veg%hc(POP%Iwood) = real(POP%pop_grid(:)%height_max)

          WHERE (pop%pop_grid(:)%LU ==2)
             casaflux%kplant(POP%Iwood,2) = 1.0_r_2 -  &
                  (1.0_r_2 -  max( min((POP%pop_grid(:)%stress_mortality + &
                  POP%pop_grid(:)%crowding_mortality ) &
                  /(POP%pop_grid(:)%cmass_sum+POP%pop_grid(:)%growth) + &
                  1.0_r_2/veg%disturbance_interval(POP%Iwood,1), 0.99_r_2), &
                  0.0_r_2))**(1.0_r_2/365.0_r_2)
          ELSEWHERE
             casaflux%kplant(POP%Iwood,2) =  1.0_r_2 -  &
                  (1.0_r_2 -  max( min((POP%pop_grid(:)%stress_mortality + &
                  POP%pop_grid(:)%crowding_mortality + &
                  POP%pop_grid(:)%cat_mortality  ) &
                  /(POP%pop_grid(:)%cmass_sum+POP%pop_grid(:)%growth), 0.99_r_2), &
                  0.0_r_2))**(1.0_r_2/365.0_r_2)
          ENDWHERE
          veg%hc(POP%Iwood) = real(POP%pop_grid(:)%height_max)
       ELSEWHERE
          casaflux%frac_sapwood(POP%Iwood) = 1.0_r_2
          casaflux%sapwood_area(POP%Iwood) = max(POP%pop_grid(:)%sapwood_area/10000._r_2, 1e-6_r_2)
          casaflux%kplant(POP%Iwood,2) = 0.0_r_2
          veg%hc(POP%Iwood) = real(POP%pop_grid(:)%height_max)
       ENDWHERE
       if (any(casaflux%kplant(:,leaf) < 0.0_r_2)) then
          do j=1, mp
             if (casaflux%kplant(j,leaf) < 0.0_r_2) then
                print*, 'ERR KK20 ', j, casaflux%kplant(j,leaf)
             endif
          enddo
       endif
       casaflux%kplant_tot(POP%Iwood,2) = casaflux%kplant(POP%Iwood,2) + &
            (1.0_r_2 -casaflux%kplant(POP%Iwood,2)) * casaflux%kplant_fire(POP%Iwood,2)
    ENDIF

    call casa_xratesoil(xklitter,xksoil,veg,soil,casamet,casabiome)
    call casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casaflux,casamet)

    IF (icycle>1) THEN
       call casa_xkN(xkNlimiting,casapool,casaflux,casamet,casabiome,veg)
       ! vh ! why is this here and not in casa_cnp.F90?

       DO j=1, mlitter
          casaflux%klitter(:,j) = casaflux%klitter(:,j)* xkNlimiting(:)
          ! MC - have to update klitter_tot as well, otherwise no C balance
          casaflux%klitter_tot(:,j) = casaflux%klitter(:,j) + &
               (1.0_r_2 -casaflux%klitter(:,j)) * casaflux%klitter_fire(:,j)
       ENDDO
       call casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
       IF (icycle >2) call casa_puptake(veg,xkNlimiting,casabiome, &
            casapool,casaflux,casamet)
    ENDIF

    ! changed by ypwang following Chris Lu on 5/nov/2012
    call casa_delplant(veg,casabiome,casapool,casaflux,casamet, &
         cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd, &
         nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd, &
         pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

    casaflux%Cplant_turnover_disturbance         = 0.0_r_2
    casaflux%Cplant_turnover_crowding            = 0.0_r_2
    casaflux%Cplant_turnover_resource_limitation = 0.0_r_2

    if (cable_user%CALL_POP) THEN
       if (.not.allocated(tmp)) allocate(tmp(size(POP%pop_grid)))
       tmp =  POP%pop_grid(:)%stress_mortality &
            + POP%pop_grid(:)%crowding_mortality &
            + POP%pop_grid(:)%cat_mortality &
            + POP%pop_grid(:)%fire_mortality
       ! where (tmp .gt. 1.e-12_r_2)
       where (tmp > real(epsilon(1.0),r_2))
          casaflux%Cplant_turnover_disturbance(POP%Iwood(:)) =  &
               casaflux%Cplant_turnover(POP%Iwood(:),2) * &
               (POP%pop_grid(:)%cat_mortality + POP%pop_grid(:)%fire_mortality) / tmp(:)
          casaflux%Cplant_turnover_crowding(POP%Iwood(:)) =  &
               casaflux%Cplant_turnover(POP%Iwood(:),2) * POP%pop_grid(:)%crowding_mortality / tmp(:)
          casaflux%Cplant_turnover_resource_limitation(POP%Iwood) = &
               casaflux%Cplant_turnover(POP%Iwood(:),2) * POP%pop_grid(:)%stress_mortality / tmp(:)
       endwhere
       deallocate(tmp)
    endif

    call casa_delsoil(veg, casapool, casaflux, casamet, casabiome)

    call casa_cnpcycle(veg, casabiome, casapool, casaflux, casamet)
    !! vh_js !!
    !CLN ndummy must be before pdummy!!!!
    IF (icycle<3) then
       IF (icycle<2) call casa_ndummy(casapool)
       call casa_pdummy(casapool)
    ENDIF

    call casa_cnpbal(casapool,casaflux,casabal,idoy)

    call casa_cnpflux(casaflux,casapool,casabal,.false.)

    ! for spinning up only
    IF (cable_user%limit_labile) THEN
       casapool%Nsoilmin = max(casapool%Nsoilmin, 5.0_r_2)
       casapool%Psoillab = max(casapool%Psoillab, 1.0_r_2)
    ENDIF

  END SUBROUTINE biogeochem


#ifndef UM_BUILD
  subroutine write_casa_restart_nc(casabiome, casamet, casapool, casaflux, &
       casabal, phen)

    use cable_common_module, only: cable_user, filename
    use casavariable, only: casa_biome, casa_met, casa_pool, casa_flux, &
         casa_balance, write_netcdf_casa_var, casafile
    use phenvariable, only: phen_variable, write_netcdf_phen_var

    implicit none

    type(casa_biome),    intent(in) :: casabiome
    type(casa_met),      intent(in) :: casamet
    type(casa_pool),     intent(in) :: casapool
    type(casa_flux),     intent(in) :: casaflux
    type(casa_balance),  intent(in) :: casabal
    type(phen_variable), intent(in) :: phen

    character(len=99) :: basename, fname

    if (len_trim(casafile%cnpepool) > 0) then
       basename = trim(casafile%cnpepool)
    else
       basename = trim(filename%path) // '/' // &
            trim(cable_user%runiden) // '_casa'
    end if

    ! casabiome
    fname = trim(basename) // '_biome_rst.nc'
    write(*,*) 'Writing casa_biome restart: ', trim(fname)
    call write_netcdf_casa_var(trim(fname), casabiome)

    ! casamet
    fname = trim(basename) // '_met_rst.nc'
    write(*,*) 'Writing casa_met restart: ', trim(fname)
    call write_netcdf_casa_var(trim(fname), casamet)

    ! casapool
    fname = trim(basename) // '_pool_rst.nc'
    write(*,*) 'Writing casa_pool restart: ', trim(fname)
    call write_netcdf_casa_var(trim(fname), casapool)

    ! casaflux
    fname = trim(basename) // '_flux_rst.nc'
    write(*,*) 'Writing casa_flux restart: ', trim(fname)
    call write_netcdf_casa_var(trim(fname), casaflux)

    ! casabal
    fname = trim(basename) // '_bal_rst.nc'
    write(*,*) 'Writing casa_bal restart: ', trim(fname)
    call write_netcdf_casa_var(trim(fname), casabal)

    ! phen
    fname = trim(basename) // '_phen_rst.nc'
    write(*,*) 'Writing casa_phen restart: ', trim(fname)
    call write_netcdf_phen_var(trim(fname), phen)

    ! write(*,*) 'Done writing casa restart files'

  end subroutine write_casa_restart_nc
#endif


#ifndef UM_BUILD
  subroutine read_casa_restart_nc(casabiome, casamet, casapool, casaflux, &
       casabal, phen)

    use cable_common_module, only: cable_user, filename
    use casavariable, only: casa_biome, casa_met, casa_pool, casa_flux, &
         casa_balance, read_netcdf_casa_var, casafile
    use phenvariable, only: phen_variable, read_netcdf_phen_var

    implicit none

    type(casa_biome),    intent(inout) :: casabiome
    type(casa_met),      intent(inout) :: casamet
    type(casa_pool),     intent(inout) :: casapool
    type(casa_flux),     intent(inout) :: casaflux
    type(casa_balance),  intent(inout) :: casabal
    type(phen_variable), intent(inout) :: phen

    character(len=99) :: basename, fname

    if (len_trim(casafile%cnpipool) > 0) then
       basename = trim(casafile%cnpipool)
    else
       basename = trim(filename%path) // '/' // &
            trim(cable_user%runiden) // '_casa'
    end if

    ! casabiome
    fname = trim(basename) // '_biome_rst.nc'
    write(*,*) 'Reading casa_biome restart: ', trim(fname)
    call read_netcdf_casa_var(trim(fname), casabiome)

    ! casamet
    fname = trim(basename) // '_met_rst.nc'
    write(*,*) 'Reading casa_met restart: ', trim(fname)
    call read_netcdf_casa_var(trim(fname), casamet)

    ! casapool
    fname = trim(basename) // '_pool_rst.nc'
    write(*,*) 'Reading casa_pool restart: ', trim(fname)
    call read_netcdf_casa_var(trim(fname), casapool)

    ! casaflux
    fname = trim(basename) // '_flux_rst.nc'
    write(*,*) 'Reading casa_flux restart: ', trim(fname)
    call read_netcdf_casa_var(trim(fname), casaflux)

    ! casabal
    fname = trim(basename) // '_bal_rst.nc'
    write(*,*) 'Reading casa_bal restart: ', trim(fname)
    call read_netcdf_casa_var(trim(fname), casabal)

    ! phen
    fname = trim(basename) // '_phen_rst.nc'
    write(*,*) 'Reading casa_phen restart: ', trim(fname)
    call read_netcdf_phen_var(trim(fname), phen)

  end subroutine read_casa_restart_nc
#endif


#ifndef UM_BUILD
  subroutine write_casa_output_nc(veg, casamet, casapool, casabal, casaflux, casaonly, ctime, lfinal)

    use casadimension,        only: mplant, mlitter, msoil, icycle
    use casavariable,         only: casa_met, casa_pool, casa_balance, casa_flux, &
         casafile, casa_timeunits
    use cable_common_module,  only: cable_user, filename, handle_err
    use cable_IO_vars_module, only: patch
    use cable_def_types_mod,  only: veg_parameter_type, mp
    ! use cable_def_types_mod,  only: sp => r_2
    use netcdf,               only: nf90_noerr, &
         nf90_put_var, nf90_clobber, nf90_create, nf90_global, nf90_put_att, &
#ifdef __NETCDF3__
         nf90_64bit_offset, &
#else
         nf90_netcdf4, nf90_classic_model, &
#endif
         nf90_def_dim, nf90_unlimited, nf90_int, nf90_def_var, nf90_float, nf90_enddef, nf90_put_var, nf90_close ! , nf90_double

    implicit none

    type(veg_parameter_type), intent(in) :: veg      ! vegetation parameters
    type(casa_met) ,          intent(in) :: casamet
    type(casa_pool),          intent(in) :: casapool
    type(casa_balance),       intent(in) :: casabal
    type(casa_flux),          intent(in) :: casaflux
    logical,                  intent(in) :: casaonly
    integer,                  intent(in) :: ctime
    logical,                  intent(in) :: lfinal

    integer :: status
    integer :: land_id, plant_id, litter_id, soil_id, t_id, i
    character(len=99) :: fname
    character(len=50) :: dum
    logical, save :: call1 = .true.
    integer, parameter :: sp = kind(1.0)

    ! 1 dim arrays (mp)
    character(len=20), dimension(3)  :: a0
    ! 2 dim arrays (mp,t)
    character(len=20), dimension(52) :: a1
    ! 3 dim arrays (mp,mplant,t)
    character(len=20), dimension(9)  :: a2
    ! 3 dim arrays (mp,mlitter,t)
    character(len=20), dimension(9)  :: a3
    ! 3 dim arrays (mp,msoil,t)
    character(len=20), dimension(8)  :: a4
    ! 4 dim arrays (mp,mlitter,mplant,t)
    character(len=20), dimension(1)  :: a5
    ! 4 dim arrays (mp,msoil,mlitter,t)
    character(len=20), dimension(1)  :: a6
    ! 4 dim arrays (mp,msoil,msoil,t)
    character(len=20), dimension(1)  :: a7

    integer, dimension(size(a0)), save :: vid0
    integer, dimension(size(a1)), save :: vid1
    integer, dimension(size(a2)), save :: vid2
    integer, dimension(size(a3)), save :: vid3
    integer, dimension(size(a4)), save :: vid4
    integer, dimension(size(a5)), save :: vid5
    integer, dimension(size(a6)), save :: vid6
    integer, dimension(size(a7)), save :: vid7
    integer, save :: vidtime, file_id, cnt
    integer :: na0, na1, na2, na3, na4, na5, na6, na7 ! actual size to write depending on icycle


    if (icycle < 1) return

    a0(1) = 'latitude'
    a0(2) = 'longitude'
    a0(3) = 'area_gridcell'
    na0 = 3

    ! C
    a1(1)  = 'glai'
    a1(2)  = 'clabile'
    a1(3)  = 'sumcbal'
    a1(4)  = 'Cgpp'
    a1(5)  = 'Cnpp'
    a1(6)  = 'stemnpp'
    a1(7)  = 'Crp'
    a1(8)  = 'Crgplant'
    a1(9)  = 'Clabloss'
    a1(10) = 'fraclabile'
    a1(11) = 'Cnep'
    a1(12) = 'Crsoil'
    a1(13) = 'FluxCtoco2'
    a1(14) = 'FCgppyear'
    a1(15) = 'FCrpyear'
    a1(16) = 'FCnppyear'
    a1(17) = 'FCrsyear'
    a1(18) = 'FCNeeyear'
    a1(19) = 'vcmax'
    na1 = 19
    ! N
    a1(20) = 'sumnbal'
    a1(21) = 'Nminfix'
    a1(22) = 'Nmindep'
    a1(23) = 'Nminloss'
    a1(24) = 'Nminleach'
    a1(25) = 'Nupland'
    a1(26) = 'Nlittermin'
    a1(27) = 'Nsmin'
    a1(28) = 'Nsimm'
    a1(29) = 'Nsnet'
    a1(30) = 'fNMinloss'
    a1(31) = 'Nsoilmin'
    if (icycle==2) na1 = 31
    ! P
    a1(32) = 'psoillab'
    a1(33) = 'psoilsorb'
    a1(34) = 'psoilocc'
    a1(35) = 'sumpbal'
    a1(36) = 'Plabuptake'
    a1(37) = 'Pdep'
    a1(38) = 'pwea'
    a1(39) = 'Pleach'
    a1(40) = 'Ploss'
    a1(41) = 'Pupland'
    a1(42) = 'Plittermin'
    a1(43) = 'Psmin'
    a1(44) = 'Psimm'
    a1(45) = 'Psnet'
    a1(46) = 'fPleach'
    a1(47) = 'kPlab'
    a1(48) = 'kPsorb'
    a1(49) = 'kpocc'
    a1(50) = 'kmlabP'
    a1(51) = 'Psorbmax'
    ! patchfrac
    a1(52) = 'patchfrac'
    if (icycle==3) na1 = 52

    ! C
    a2(1) = 'cplant'
    a2(2) = 'fracCalloc'
    a2(3) = 'kplant'
    a2(4) = 'Crmplant'
    a2(5) = 'kplant_fire'
    na2 = 5
    ! N
    a2(6) = 'nplant'
    a2(7) = 'fracNalloc'
    if (icycle==2) na2 = 7
    ! P
    a2(8) = 'pplant'
    a2(9) = 'fracPalloc'
    if (icycle==3) na2 = 9

    ! C
    a3(1) = 'clitter'
    a3(2) = 'klitter'
    a3(3) = 'fromLtoCO2'
    a3(4) = 'FluxCtolitter'
    a3(5) = 'klitter_fire'
    na3 = 5
    ! N
    a3(6) = 'nlitter'
    a3(7) = 'FluxNtolitter'
    if (icycle==2) na3 = 7
    ! P
    a3(8) = 'plitter'
    a3(9) = 'FluxPtolitter'
    if (icycle==3) na3 = 9

    ! C
    a4(1) = 'csoil'
    a4(2) = 'ksoil'
    a4(3) = 'fromStoCO2'
    a4(4) = 'FluxCtosoil'
    na4 = 4
    ! N
    a4(5) = 'nsoil'
    a4(6) = 'FluxNtosoil'
    if (icycle==2) na4 = 6
    ! P
    a4(7) = 'psoil'
    a4(8) = 'FluxPxtosoil'
    if (icycle==3) na4 = 8

    ! C
    a5(1) = 'fromPtoL'
    na5 = 1

    ! C
    a6(1) = 'fromLtoS'
    na6 = 1

    ! C
    a7(1) = 'fromStoS'
    na7 = 1

    ! Get File-Name
    if (len_trim(casafile%out) > 0) then
       fname = trim(casafile%out)
    else
       if (len_trim(cable_user%mettype) > 0) then
          if (cable_user%yearstart < 1000) then
             write(dum, fmt="(i3)") cable_user%yearstart
          else
             write(dum, fmt="(i4)") cable_user%yearstart
          endif
          if (cable_user%yearend < 1000) then
             write(dum, fmt="(a,a,i3)") trim(dum), '_', cable_user%yearend
          else
             write(dum, fmt="(a,a,i4)") trim(dum), '_', cable_user%yearend
          endif
          fname = trim(filename%path)//'/'//trim(cable_user%RunIden)//'_'//trim(dum)//'_casa_out.nc'
       else
          ! site data
          fname = trim(filename%path)//'/'//trim(cable_user%RunIden)//'_casa_out.nc'
       endif
    endif

    if ( call1 ) then
       cnt = 0

       ! create netcdf file:
#ifdef __NETCDF3__
       status = nf90_create(trim(fname), ior(nf90_clobber,nf90_64bit_offset), file_id)
#else
       status = nf90_create(trim(fname), ior(nf90_clobber,ior(nf90_netcdf4,nf90_classic_model)), file_id)
#endif
       if (status /= nf90_noerr) call handle_err(status)

       status = nf90_put_att(file_id, nf90_global, "icycle"   , icycle)
       status = nf90_put_att(file_id, nf90_global, "startyear", cable_user%yearstart)
       status = nf90_put_att(file_id, nf90_global, "endyear"  , cable_user%yearend)
       status = nf90_put_att(file_id, nf90_global, "runiden"  , cable_user%runiden)
       if (casaonly) then
          dum = 'casa-only run'
       else
          dum = 'cable-casa coupled run'
       endif
       status = nf90_put_att(file_id, nf90_global, "run-type", trim(dum))

       ! define dimensions:
       ! land (number of points)
       status = nf90_def_dim(file_id, 'land',    mp,             land_id)
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_def_dim(file_id, 'mplant',  mplant,         plant_id)
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_def_dim(file_id, 'mlitter', mlitter,        litter_id)
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_def_dim(file_id, 'msoil',   msoil,          soil_id)
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_def_dim(file_id, 'time',    nf90_unlimited, t_id)
       if (status /= nf90_noerr) call handle_err(status)

       ! define variables
       status = nf90_def_var(file_id, 'time', nf90_int, (/t_id/), vidtime)
       if (status /= nf90_noerr) call handle_err(status)

       STATUS = NF90_PUT_ATT(FILE_ID, VIDtime, 'units', TRIM(casa_timeunits))
       IF (STATUS /= NF90_NOERR)  CALL handle_err(STATUS)

       ! STATUS = NF90_PUT_ATT(FILE_ID, VIDtime, 'calendar', trim(casa_calendar))
       ! IF (STATUS /= NF90_NOERR)  CALL handle_err(STATUS)

       do i=1, na0
          status = nf90_def_var(file_id, trim(a0(i)), nf90_float, (/land_id/), vid0(i))
          if (status /= nf90_noerr) call handle_err(status)
       end do

       do i=1, na1
          status = nf90_def_var(file_id, trim(a1(i)), nf90_float, (/land_id,t_id/), vid1(i) &
#ifndef __NETCDF3__
               , deflate_level=4 &
#endif
               )
          if (status /= nf90_noerr) call handle_err(status)
       end do

       do i=1, na2
          status = nf90_def_var(file_id, trim(a2(i)), nf90_float, (/land_id,plant_id,t_id/), vid2(i) &
#ifndef __NETCDF3__
               , deflate_level=4 &
#endif
               )
          if (status /= nf90_noerr) call handle_err(status)
       end do

       do i=1, na3
          status = nf90_def_var(file_id, trim(a3(i)), nf90_float, (/land_id,litter_id,t_id/), vid3(i) &
#ifndef __NETCDF3__
               , deflate_level=4 &
#endif
               )
          if (status /= nf90_noerr) call handle_err(status)
       end do

       do i=1, na4
          status = nf90_def_var(file_id, trim(a4(i)), nf90_float, (/land_id,soil_id,t_id/), vid4(i) &
#ifndef __NETCDF3__
               , deflate_level=4 &
#endif
               )
          if (status /= nf90_noerr) call handle_err(status)
       end do

       do i=1, na5
          status = nf90_def_var(file_id, trim(a5(i)), nf90_float, (/land_id,litter_id,plant_id,t_id/), vid5(i) &
#ifndef __NETCDF3__
               , deflate_level=4 &
#endif
               )
          if (status /= nf90_noerr) call handle_err(status)
       end do

       do i=1, na6
          status = nf90_def_var(file_id, trim(a6(i)), nf90_float, (/land_id,soil_id,litter_id,t_id/), vid6(i) &
#ifndef __NETCDF3__
               , deflate_level=4 &
#endif
               )
          if (status /= nf90_noerr) call handle_err(status)
       end do

       do i=1, na7
          status = nf90_def_var(file_id, trim(a7(i)), nf90_float, (/land_id,soil_id,soil_id,t_id/), vid7(i) &
#ifndef __NETCDF3__
               , deflate_level=4 &
#endif
               )
          if (status /= nf90_noerr) call handle_err(status)
       end do

       ! end define mode:
       status = nf90_enddef(file_id)
       if (status /= nf90_noerr) call handle_err(status)

       ! put lat / lon ( mp )
       status = nf90_put_var(file_id, vid0(1), real(casamet%lat,sp))
       if(status /= nf90_noerr) call handle_err(status)

       status = nf90_put_var(file_id, vid0(2), real(casamet%lon,sp))
       if(status /= nf90_noerr) call handle_err(status)

       status = nf90_put_var(file_id, vid0(3), real(casamet%areacell,sp))
       if(status /= nf90_noerr) call handle_err(status)

       call1 = .false.
    endif ! call1

    cnt = cnt + 1

    ! time  ( t )
    status = nf90_put_var(file_id, vidtime, ctime, start=(/cnt/))
    if(status /= nf90_noerr) call handle_err(status)

    ! C
    status = nf90_put_var(file_id, vid1(1),  real(casamet%glai,sp),         start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(2),  real(casapool%clabile,sp),     start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(3),  real(casabal%sumcbal,sp),      start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(4),  real(casaflux%cgpp,sp),        start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(5),  real(casaflux%cnpp,sp),        start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(6),  real(casaflux%stemnpp,sp),     start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(7),  real(casaflux%crp,sp),         start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(8),  real(casaflux%crgplant,sp),    start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(9),  real(casaflux%clabloss,sp),    start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(10), real(casaflux%fracclabile,sp), start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(11), real(casaflux%cnep,sp),        start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(12), real(casaflux%crsoil,sp),      start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(13), real(casaflux%fluxctoco2,sp),  start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(14), real(casabal%fcgppyear,sp),    start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(15), real(casabal%fcrpyear,sp),     start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(16), real(casabal%fcnppyear,sp),    start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(17), real(casabal%fcrsyear,sp),     start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(18), real(casabal%fcneeyear,sp),    start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid1(19), real(veg%vcmax,sp),            start=(/1,cnt/), count=(/mp,1/) )
    if(status /= nf90_noerr) call handle_err(status)
    ! N
    if (icycle > 1) then
       status = nf90_put_var(file_id, vid1(20), real(casabal%sumnbal,sp),      start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(21), real(casaflux%nminfix,sp),     start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(22), real(casaflux%nmindep,sp),     start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(23), real(casaflux%nminloss,sp),    start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(24), real(casaflux%nminleach,sp),   start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(25), real(casaflux%nupland,sp),     start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(26), real(casaflux%nlittermin,sp),  start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(27), real(casaflux%nsmin,sp),       start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(28), real(casaflux%nsimm,sp),       start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(29), real(casaflux%nsnet,sp),       start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(30), real(casaflux%fnminloss,sp),   start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(31), real(casapool%nsoilmin,sp),    start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
    endif
    ! P
    if (icycle > 2) then
       status = nf90_put_var(file_id, vid1(32), real(casapool%psoillab,sp),    start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(33), real(casapool%psoilsorb,sp),   start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(34), real(casapool%psoilocc,sp),    start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(35), real(casabal%sumpbal,sp),      start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(36), real(casaflux%plabuptake,sp),  start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(37), real(casaflux%pdep,sp),        start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(38), real(casaflux%pwea,sp),        start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(39), real(casaflux%pleach,sp),      start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(40), real(casaflux%ploss,sp),       start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(41), real(casaflux%pupland,sp),     start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(42), real(casaflux%plittermin,sp),  start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(43), real(casaflux%psmin,sp),       start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(44), real(casaflux%psimm,sp),       start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(45), real(casaflux%psnet,sp),       start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(46), real(casaflux%fpleach,sp),     start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(47), real(casaflux%kplab,sp),       start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(48), real(casaflux%kpsorb,sp),      start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(49), real(casaflux%kpocc,sp),       start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(50), real(casaflux%kmlabp,sp),      start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid1(51), real(casaflux%psorbmax,sp),    start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
       ! patchfrac
       status = nf90_put_var(file_id, vid1(52), real(patch(:)%frac,sp),         start=(/1,cnt/), count=(/mp,1/) )
       if(status /= nf90_noerr) call handle_err(status)
    endif

    ! put 3d vars ( mp, mplant, t )
    ! C
    status = nf90_put_var(file_id, vid2(1), real(casapool%cplant,sp),      start=(/1,1,cnt/), count=(/mp,mplant,1/))
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid2(2), real(casaflux%fracCalloc,sp),  start=(/1,1,cnt/), count=(/mp,mplant,1/))
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid2(3), real(casaflux%kplant,sp),      start=(/1,1,cnt/), count=(/mp,mplant,1/))
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid2(4), real(casaflux%crmplant,sp),    start=(/1,1,cnt/), count=(/mp,mplant,1/))
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid2(5), real(casaflux%kplant_fire,sp), start=(/1,1,cnt/), count=(/mp,mplant,1/))
    if(status /= nf90_noerr) call handle_err(status)
    ! N
    if (icycle > 1) then
       status = nf90_put_var(file_id, vid2(6), real(casapool%nplant,sp),     start=(/1,1,cnt/), count=(/mp,mplant,1/))
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid2(7), real(casaflux%fracnalloc,sp), start=(/1,1,cnt/), count=(/mp,mplant,1/))
       if(status /= nf90_noerr) call handle_err(status)
    endif
    ! P
    if (icycle > 2) then
       status = nf90_put_var(file_id, vid2(8), real(casapool%pplant,sp),     start=(/1,1,cnt/), count=(/mp,mplant,1/))
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid2(9), real(casaflux%fracpalloc,sp), start=(/1,1,cnt/), count=(/mp,mplant,1/))
       if(status /= nf90_noerr) call handle_err(status)
    endif

    ! put 3d vars ( mp, mlitter, t )
    ! C
    status = nf90_put_var(file_id, vid3(1), real(casapool%clitter,sp),       start=(/1,1,cnt/), count=(/mp,mlitter,1/))
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid3(2), real(casaflux%klitter,sp),       start=(/1,1,cnt/), count=(/mp,mlitter,1/))
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid3(3), real(casaflux%fromltoco2,sp),    start=(/1,1,cnt/), count=(/mp,mlitter,1/))
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid3(4), real(casaflux%fluxctolitter,sp), start=(/1,1,cnt/), count=(/mp,mlitter,1/))
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid3(5), real(casaflux%klitter_fire,sp),  start=(/1,1,cnt/), count=(/mp,mlitter,1/))
    if(status /= nf90_noerr) call handle_err(status)
    ! N
    if (icycle > 1) then
       status = nf90_put_var(file_id, vid3(6), real(casapool%nlitter,sp),       start=(/1,1,cnt/), count=(/mp,mlitter,1/))
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid3(7), real(casaflux%fluxntolitter,sp), start=(/1,1,cnt/), count=(/mp,mlitter,1/))
       if(status /= nf90_noerr) call handle_err(status)
    endif
    ! P
    if (icycle > 2) then
       status = nf90_put_var(file_id, vid3(8), real(casapool%plitter,sp),       start=(/1,1,cnt/), count=(/mp,mlitter,1/))
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid3(9), real(casaflux%fluxptolitter,sp), start=(/1,1,cnt/), count=(/mp,mlitter,1/))
       if(status /= nf90_noerr) call handle_err(status)
    endif

    ! put 3d vars ( mp, msoil, t )
    ! C
    status = nf90_put_var(file_id, vid4(1), real(casapool%csoil,sp),       start=(/1,1,cnt/), count=(/mp,msoil,1/))
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid4(2), real(casaflux%ksoil,sp),       start=(/1,1,cnt/), count=(/mp,msoil,1/))
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid4(3), real(casaflux%fromstoco2,sp),  start=(/1,1,cnt/), count=(/mp,msoil,1/))
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_put_var(file_id, vid4(4), real(casaflux%fluxctosoil,sp), start=(/1,1,cnt/), count=(/mp,msoil,1/))
    if(status /= nf90_noerr) call handle_err(status)
    ! N
    if (icycle > 1) then
       status = nf90_put_var(file_id, vid4(5), real(casapool%nsoil,sp),       start=(/1,1,cnt/), count=(/mp,msoil,1/))
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid4(6), real(casaflux%fluxntosoil,sp), start=(/1,1,cnt/), count=(/mp,msoil,1/))
       if(status /= nf90_noerr) call handle_err(status)
    endif
    ! P
    if (icycle > 2) then
       status = nf90_put_var(file_id, vid4(7), real(casapool%psoil,sp),       start=(/1,1,cnt/), count=(/mp,msoil,1/))
       if(status /= nf90_noerr) call handle_err(status)
       status = nf90_put_var(file_id, vid4(8), real(casaflux%fluxptosoil,sp), start=(/1,1,cnt/), count=(/mp,msoil,1/))
       if(status /= nf90_noerr) call handle_err(status)
    endif

    ! put 4d vars ( mp, mlitter,mplant, t )
    ! C
    status = nf90_put_var(file_id, vid5(1), real(casaflux%fromptol,sp), start=(/1,1,1,cnt/), count=(/mp,mlitter,mplant,1/))
    if(status /= nf90_noerr) call handle_err(status)

    ! put 4d vars ( mp, msoil, mlitter, t )
    ! C
    status = nf90_put_var(file_id, vid6(1), real(casaflux%fromltos,sp), start=(/1,1,1,cnt/), count=(/mp,msoil,mlitter,1/))
    if(status /= nf90_noerr) call handle_err(status)

    ! put 4d vars ( mp, msoil, msoil, t )
    ! C
    status = nf90_put_var(file_id, vid7(1), real(casaflux%fromstos,sp), start=(/1,1,1,cnt/), count=(/mp,msoil,msoil,1/))
    if(status /= nf90_noerr) call handle_err(status)

    if ( lfinal ) then
       ! close netcdf file:
       status = nf90_close(file_id)
       file_id = -1
       if (status /= nf90_noerr) call handle_err(status)
       write(*,*) " casa output written to ", trim(fname)
    endif

  end subroutine write_casa_output_nc
#endif

end module casa_inout
