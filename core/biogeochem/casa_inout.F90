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
!
! ==============================================================================
! casa_inout.f90
!
! the following routines are used when "casacnp" is coupled to "cable"
!   casa_readbiome
!   casa_readphen
!   casa_readpoint   (removed, now done in parameter_module)
!   casa_init
!   casa_poolout
!   casa_cnpflux  (zeros casabal quantites on doy 1 and updates casabal at end of biogeochem)
!   biogeochem

!#define UM_BUILD YES
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
  REAL(r_2), DIMENSION(mvtype)       :: totroot
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
  REAL(r_2), DIMENSION(mvtype)       :: xkuplabp
  REAL(r_2), DIMENSION(mvtype,ms)    :: fracroot
  REAL(r_2) ::  xratioNPleafmin,xratioNPleafmax,         &
                xratioNPwoodmin,xratioNPwoodmax,         &
                xratioNPfrootmin,xratioNPfrootmax
  INTEGER :: i,iv1,nv,ns,npt,iv,is,iso
  INTEGER :: nv0,nv1,nv2,nv3,nv4,nv5,nv6,nv7,nv8,nv9,nv10,nv11,nv12,nv13
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
!     PRINT *, nv,nv0,casabiome%ivt2(nv)
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
!     PRINT *, 'nv1',nv,nv1
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv2, &
                casabiome%fracnpptoP(nv,leaf),casabiome%fracnpptoP(nv,wood), &
                casabiome%fracnpptoP(nv,froot),casabiome%rmplant(nv,leaf),   &
                casabiome%rmplant(nv,wood),casabiome%rmplant(nv,froot)
!     PRINT *, 'nv2', nv2
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
!         xfherbivore(nv),casabiome%ratiofrootleaf(nv),                  &
         casabiome%glaimax(nv),casabiome%glaimin(nv)
!     PRINT *, 'nv22',nv2
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv3, cleaf(nv),cwood(nv),cfroot(nv),cmet(nv),   &
                cstr(nv),ccwd(nv), cmic(nv), cslow(nv),cpass(nv)
!     PRINT *, 'nv3',nv3
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv4, &
         phen%TKshed(nv),xxkleafcoldmax(nv),casabiome%xkleafcoldexp(nv), &
         xxkleafdrymax(nv),casabiome%xkleafdryexp(nv)
!     PRINT *, 'nv4',nv4
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
!     PRINT *, 'nv6',nv6
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv7,nleaf(nv),nwood(nv),nfroot(nv), &
                nmet(nv),nstr(nv), ncwd(nv), &
                nmic(nv),nslow(nv),npass(nv),xnsoilmin(nv)
!     PRINT *, 'nv7',nv7
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
!     PRINT *, 'nv9',nv9
  ENDDO
  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv10, &
         xpleaf(nv),xpwood(nv),xpfroot(nv),xpmet(nv),xpstr(nv),xpcwd(nv), &
         xpmic(nv),xpslow(nv),xppass(nv),xplab(nv),xpsorb(nv),xpocc(nv)
!     PRINT *, 'nv10',nv10
  ENDDO

 !@@@@@@@@@@@@@@@@@@@@@@@@@
  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv11, &
         xxnpmax(nv),xq10soil(nv),xxkoptlitter(nv),xxkoptsoil(nv),xprodptase(nv), &
         xcostnpup(nv),xmaxfinelitter(nv),xmaxcwd(nv),xnintercept(nv),xnslope(nv)
  ENDDO
  !@@@@@@@@@@@@@@@@@@@@@

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

  fracroot   = 0.0
  depthsoila = 0.0
  depthsoilb = 0.0
  DO ns=1,ms
    depthsoilb(ns) = depthsoilb(ns) + soil%zse(ns)
    IF (ns==1) THEN
      depthsoila(ns) = 0.0
    ELSE
      depthsoila(ns) = depthsoilb(ns-1)
    ENDIF
  ENDDO

  DO nv=1,mvtype
    casabiome%sla(nv)             = slax(nv)
    casabiome%fraclabile(nv,leaf) = deltcasa*0.6    !1/day
    casabiome%fraclabile(nv,froot)= deltcasa*0.4    !1/day
    casabiome%fraclabile(nv,wood) = deltcasa*0.0
    casabiome%plantrate(nv,leaf)  = deltcasa/(leafage(nv)*(1.0-xfherbivore(nv)))
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
!    casabiome%kuplabp(nv)         = xkuplabp(nv)
    casabiome%rmplant(nv,:)       = casabiome%rmplant(nv,:)*deltcasa
    casabiome%kclabrate(nv)       = deltcasa/clabileage(nv)

!@@@@@@@@@@@@@@@@@
    casabiome%xnpmax(nv)          = xxnpmax(nv)
    casabiome%q10soil(nv)         = xq10soil(nv)
    casabiome%xkoptlitter(nv)     = xxkoptlitter(nv)
    casabiome%xkoptsoil(nv)       = xxkoptsoil(nv)
    casabiome%prodptase(nv)       = xprodptase(nv)/365.0   ! convert from yearly to daily
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
!@@@@@@@@@@@@@@
  ENDDO

!@@@@@@@@@@@@@@
  DO ns=1,mso
    casabiome%xkplab(ns)          =  xxkplab(ns)
    casabiome%xkpsorb(ns)         =  xxkpsorb(ns)
    casabiome%xkpocc(ns)          =  xxkpocc(ns)
  ENDDO

!@@@@@@@@@@@@@@

 ! PRINT *, 'casabiome%xkoptsoil = ', casabiome%xkoptsoil(2)

  DO npt = 1, mp
    iv1=veg%iveg(npt)
    iso=casamet%isorder(npt)
    ! The following to be commented out when coupled to CABLE
!    veg%froot(npt,:) =fracroot(iv1,:)
!    PRINT *, 'npt,iv1,iso ', npt,iv1, iso
    casamet%iveg2(npt) =casabiome%ivt2(iv1)
    casamet%lnonwood(npt) = 1
    casapool%cplant(npt,wood)  = 0.0
    casapool%clitter(npt,cwd)  = 0.0
    casapool%nplant(npt,wood)  = 0.0
    casapool%nlitter(npt,cwd)  = 0.0
    casapool%pplant(npt,wood)  = 0.0
    casapool%plitter(npt,cwd)  = 0.0
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
         casapool%cplant(npt,wood) = 0.01
         casapool%nplant(npt,wood)= casabiome%ratioNCplantmin(iv1,wood)* casapool%cplant(npt,wood)
         casapool%pplant(npt,wood)= casabiome%ratioPCplantmin(iv1,wood)* casapool%cplant(npt,wood)
      ENDIF
      !! vh_js

    ENDIF
    casapool%cplant(npt,leaf)     = cleaf(iv1)
    casapool%cplant(npt,froot)    = cfroot(iv1)
    casapool%clabile(npt)         = 0.0
    casapool%clitter(npt,metb)     = cmet(iv1)
    casapool%clitter(npt,str)     = cstr(iv1)
    casapool%csoil(npt,mic)       = cmic(iv1)
    casapool%csoil(npt,slow)      = cslow(iv1)
    casapool%csoil(npt,pass)      = cpass(iv1)
    IF (icycle==1) THEN
      casapool%ratioNCplant(npt,:)  = 1.0/ratioCNplant(iv1,:)
    ENDIF

    ! initializing glai in case not reading pool file (eg. during spin)
    casamet%glai(npt) = MAX(casabiome%glaimin(iv1), &
                        casabiome%sla(iv1) * casapool%cplant(npt,leaf))

    casaflux%fNminloss(npt)   = xfNminloss(iv1)
    ! comment out by ypw 12/07/2009
    casaflux%fNminleach(npt)  = 10.0*xfNminleach(iv1) * deltcasa
!    casaflux%fNminleach(npt)  = xfNminleach(iv1)
    casapool%nplant(npt,leaf) = nleaf(iv1)
    casapool%nplant(npt,froot)= nfroot(iv1)
    casapool%nlitter(npt,metb) = nmet(iv1)
!    casapool%nlitter(npt,str) = nstr(iv1)
    casapool%nlitter(npt,str) = cstr(iv1)*ratioNCstrfix
    casapool%nsoil(npt,mic)   = nmic(iv1)
    casapool%nsoil(npt,slow)  = nslow(iv1)
    casapool%nsoil(npt,pass)  = npass(iv1)
    casapool%nsoilmin(npt)    = xnsoilmin(iv1)
    casapool%pplant(npt,leaf) = xpleaf(iv1)
    casapool%pplant(npt,froot)= xpfroot(iv1)
    casapool%plitter(npt,metb) = xpmet(iv1)
!    casapool%plitter(npt,str) = xpstr(iv1)
    casapool%plitter(npt,str) = casapool%nlitter(npt,str)/ratioNPstrfix
    casapool%psoil(npt,mic)   = xpmic(iv1)
    casapool%psoil(npt,slow)  = xpslow(iv1)
    casapool%psoil(npt,pass)  = xppass(iv1)
    casapool%psoillab(npt)    = xplab(iv1)
    casapool%psoilsorb(npt)   = xpsorb(iv1)
    casapool%psoilocc(npt)    = xpocc(iv1)
    casapool%Psoillab(npt)    = 0.0
    casaflux%kmlabp(npt)      = xkmlabp(iso)
    casaflux%psorbmax(npt)    = xpsorbmax(iso)
    casaflux%fpleach(npt)     = xfPleach(iso) /(365.0)    ! convert from 1/year to 1/day
!   we used the spatially explicit estimate N fixation by Wang and Houlton (GRL)
!    casaflux%Nminfix(npt)     = xnfixrate(iv1)/365.0

    casapool%ratioNCplant(npt,:)  = 1.0/ratioCNplant(iv1,:)
    casapool%ratioNPplant(npt,:)  = casabiome%ratioNPplantmin(iv1,:)
    casapool%ratioNClitter(npt,:) = casapool%nlitter(npt,:)/(casapool%clitter(npt,:)+1.0e-10)
    casapool%ratioNPlitter(npt,:) = casapool%nlitter(npt,:)/(casapool%plitter(npt,:)+1.0e-10)
    casapool%ratioNCsoil(npt,:)   = 1.0/ratioCNsoil(iv1,:)
    casapool%ratioNPsoil(npt,:)   = ratioNPsoil(iso,:)
    casapool%ratioNCsoilmin(npt,:)   = 1.0/ratioCNsoilmax(iv1,:)
    casapool%ratioNCsoilmax(npt,:)   = 1.0/ratioCNsoilmin(iv1,:)
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

!  DO npt=1,mp
!    IF (veg%iveg(npt)==12) PRINT *, npt, veg%iveg(npt), &
!         casapool%Psoil(npt,:),casapool%psoilsorb(npt), &
!         casaflux%psorbmax(npt),casapool%psoillab(npt),casaflux%kmlabp(npt)
!  ENDDO

END SUBROUTINE casa_readbiome

SUBROUTINE casa_readphen(veg,casamet,phen)
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
  TYPE (veg_parameter_type), INTENT(IN)    :: veg  ! vegetation parameters
  TYPE (casa_met),           INTENT(IN)    :: casamet
  TYPE (phen_variable),      INTENT(INOUT) :: phen

  ! local variables
  INTEGER, PARAMETER            :: nphen=8! was 10(IGBP). changed by Q.Zhang @01/12/2011
  INTEGER np,nx,ilat
  INTEGER, DIMENSION(271,mvtype) :: greenup, fall,  phendoy1
  INTEGER, DIMENSION(nphen)     :: greenupx,fallx,xphendoy1
  INTEGER, DIMENSION(nphen)     :: ivtx
  REAL(r_2), DIMENSION(271)     :: xlat

  ! initilize for evergreen PFTs
  greenup(:,:) = -50
  fall(:,:)    = 367
  phendoy1(:,:)= 2

  OPEN(101,file=casafile%phen)
  READ(101,*)
  READ(101,*) (ivtx(nx),nx=1,nphen) ! fixed at 10, as only 10 of 17 IGBP PFT
                                    ! have seasonal leaf phenology
  DO ilat=271,1,-1
    READ(101,*) xlat(ilat),(greenupx(nx),nx=1,nphen), &
                (fallx(nx),nx=1,nphen),(xphendoy1(nx),nx=1,nphen)
    DO nx=1,nphen
      greenup(ilat,ivtx(nx)) = greenupx(nx)
      fall(ilat,ivtx(nx))    = fallx(nx)
      phendoy1(ilat,ivtx(nx))= xphendoy1(nx)
    ENDDO
  ENDDO

  DO np=1,mp
    ilat=(casamet%lat(np)+55.25)/0.5+1
    ilat= MIN(271,MAX(1,ilat))
    phen%phase(np) = phendoy1(ilat,veg%iveg(np))
    phen%doyphase(np,1) = greenup(ilat,veg%iveg(np)) ! DOY for greenup
    phen%doyphase(np,2) = phen%doyphase(np,1) +14    ! DOY for steady LAI
    phen%doyphase(np,3) = fall(ilat,veg%iveg(np))    ! DOY for leaf senescence
    phen%doyphase(np,4) = phen%doyphase(np,3) +14    ! DOY for minimal LAI season
    IF (phen%doyphase(np,2) > 365) phen%doyphase(np,2)=phen%doyphase(np,2)-365
    IF (phen%doyphase(np,4) > 365) phen%doyphase(np,4)=phen%doyphase(np,4)-365

 ENDDO

END SUBROUTINE casa_readphen

!SUBROUTINE casa_readpoint(veg,soil,casaflux,casamet,rad)
!! Transfer grid information from CABLE internally, read N&P input from
!! integral NETCDF file "cnpdata_r21.nc" (Q.Zhang 01/08/2011)
!
!!SUBROUTINE casa_readpoint(mvt,veg,soil,casaflux,casamet,patch,rad)
!  USE netcdf
!  USE cable_def_types_mod
!  USE abort_module
!  USE io_variables, ONLY: landpt,patch        ! add landpt, Q.Zhang 05/08/2011
!  USE casaparm
!  USE casadimension
!  USE casavariable
!  IMPLICIT NONE
!!  INTEGER,               INTENT(IN)    :: mvt
!  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
!  TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
!  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
!  TYPE (casa_met),            INTENT(INOUT) :: casamet
!  TYPE (radiation_type),      INTENT(IN)    :: rad
!
!  ! local variables
!  INTEGER, DIMENSION(:,:), ALLOCATABLE :: iso
!  REAL,DIMENSION(:,:), ALLOCATABLE :: annNdep,annNfix,annPwea,annPdust
!  REAL,DIMENSION(:), ALLOCATABLE:: latx, lonx
!  INTEGER :: nlat, nlon, ii, jj, g, p
!  INTEGER :: ncid, ok, varid
!
!  ok = NF90_OPEN(casafile%cnppoint,0,ncid)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error opening cnpdata_r21.nc')
!
!  ok = NF90_INQ_DIMID(ncid,'lon',varid)
!  ok = NF90_INQUIRE_DIMENSION(ncid,varid,LEN=nlon)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error getting longitude')
!  ok = NF90_INQ_DIMID(ncid,'lat',varid)
!  ok = NF90_INQUIRE_DIMENSION(ncid,varid,LEN=nlat)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error getting latitude')
!
!  ALLOCATE( iso(nlon,nlat) )
!  ALLOCATE( annNdep(nlon,nlat) )
!  ALLOCATE( annNfix(nlon,nlat) )
!  ALLOCATE( annPwea(nlon,nlat) )
!  ALLOCATE( annPdust(nlon,nlat) )
!  ALLOCATE( latx(nlat) )
!  ALLOCATE( lonx(nlon) )
!
!  ! Read temporary variables
!  ! soil order
!  ok = NF90_INQ_VARID(ncid,'sorder',varid)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable sorder')
!  ok = NF90_GET_VAR(ncid,varid,iso)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable sorder')
!  ! N deposition
!  ok = NF90_INQ_VARID(ncid,'ndep',varid)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable ndep')
!  ok = NF90_GET_VAR(ncid,varid,annNdep)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable ndep')
!  ! N fixation rate
!  ok = NF90_INQ_VARID(ncid,'nfix',varid)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable nfix')
!  ok = NF90_GET_VAR(ncid,varid,annNfix)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable nfix')
!  ! P dust deposition
!  ok = NF90_INQ_VARID(ncid,'pdust',varid)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable pdust')
!  ok = NF90_GET_VAR(ncid,varid,annPdust)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable pdust')
!  ! P weathering rate
!  ok = NF90_INQ_VARID(ncid,'pweather',varid)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding variable pweather')
!  ok = NF90_GET_VAR(ncid,varid,annPwea)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable pweather')
!  ! lat and lon
!  ok = NF90_INQ_VARID(ncid,'lat',varid)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding lat')
!  ok = NF90_GET_VAR(ncid,varid,latx)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable lat')
!  ok = NF90_INQ_VARID(ncid,'lon',varid)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding lon')
!  ok = NF90_GET_VAR(ncid,varid,lonx)
!  IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading variable lon')
!
!  ! reorder lon from (0,360) to (-180,180)
!  where (lonx > 180.) lonx = lonx - 360.
!
!  do g = 1,mland
!    ii = landpt(g)%ilon
!    jj = landpt(g)%ilat
!
!   do p = landpt(g)%cstart,landpt(g)%cend
!     casamet%lon(p) = lonx(ii)
!     casamet%lat(p) = latx(jj)
!     if (ABS(casamet%lat(p) - patch(p)%latitude) > 0.1 .or. &
!        ABS(casamet%lon(p) - patch(p)%longitude) > 0.1) then
!       print*, "check nutrient input, coordinate unmatch"
!       print*, p, casamet%lon(p), patch(p)%longitude
!       print*, p, casamet%lat(p), patch(p)%latitude
!       stop
!     end if
!
!     casamet%isorder(p)  = iso(ii,jj)
!     casaflux%Nmindep(p) = annNdep(ii,jj)/365.0*1.e-3! gN/m2/day
!     casaflux%Nminfix(p) = annNfix(ii,jj)/365.0      ! gN/m2/day
!     casaflux%Pdep(p)    = annPdust(ii,jj)/365.0     ! gP/m2/day
!     casaflux%Pwea(p)    = annPwea(ii,jj)/365.0      ! gP/m2/day
!
!     if(veg%iveg(p)==9 .or. veg%iveg(p)==10) then
!     ! P fertilizer =13 Mt P globally in 1994
!       casaflux%Pdep(p) = casaflux%Pdep(p)+0.7/365.0
!     ! N fertilizer =86 Mt N globally in 1994
!       casaflux%Nminfix(p) = casaflux%Nminfix(p)+4.3/365.0
!     endif
!   end do
!  end do
!
!  ok = NF90_CLOSE(ncid)
!  if (ok /= NF90_NOERR) CALL nc_abort(ok,'error closing cnpdata_r21.nc')
!
!  DEALLOCATE(latx,lonx,iso,annNdep,annNfix,annPwea,annPdust)
!
!!  ! local variables
!!  INTEGER :: np,nland
!!  REAL(r_2) :: annNdep,annNfix,annPwea,annPdust
!!  REAL(r_2) :: annNfert,annPfert   ! not really used yet
!!  INTEGER, DIMENSION(mp) :: vtypex,stypex
!!  INTEGER :: nlandx,ivtigbp,inPatch,ilat,ilon
!!  REAL    :: frac,ssat,swilt,sfc   ! used in offline version, Q.Zhang @ 25/02/2011
!!
!!  OPEN(101,file=casafile%cnppoint,FORM='FORMATTED')
!!  READ(101,*)
!!!  READ(101,*) inPatch
!!  PRINT * ,'Within casa_readpoint, mp = ', mp
!!!  PRINT * ,'Input file has ', inPatch, ' patches.'
!!
!!  np = 0
!!  DO nland=1,mp
!!    np = np + 1
!!!    READ(101,*) &
!!!        nlandx,ivtigbp,stypex(np),casamet%isorder(np), &
!!!               casamet%lat(np),casamet%lon(np),casamet%areacell(np), &
!!!               annNfix,annNdep,annNfert,annPwea,annPdust,annPfert
!!
!!    ! 'ijgcm,j,i,lat,lon,frac,iveg,isoil,ist,parea,ssat,swilt,sfc,ndep,nfix,pwea,pdust'
!!    read(101,*) nlandx,ilat,ilon,casamet%lat(np),casamet%lon(np),&
!!                frac,vtypex(np),stypex(np),casamet%isorder(np),&
!!                casamet%areacell(np),ssat,swilt,&
!!                sfc,annNdep,annNfix,annPwea,annPdust
!!
!!!    PRINT * , nlandx,ivtigbp,stypex(np),veg%iveg(np),soil%isoilm(np), &
!!!              patch(np)%frac,patch(np)%latitude,patch(np)%longitude, &
!!!              casamet%lat(np),casamet%lon(np)
!!
!!!    IF (ivtigbp == 0) ivtigbp = iceland
!!
!!    IF (ABS(casamet%lat(np) - patch(np)%latitude) < 0.1 .AND. &
!!        ABS(casamet%lon(np) - patch(np)%longitude) < 0.1) THEN
!!      IF (vtypex(np) /= veg%iveg(np) .OR. stypex(np) /= soil%isoilm(np)) THEN
!!        PRINT * ,'Check why iveg, isoil do not match'
!!        STOP
!!      ELSE
!!        casaflux%Nmindep(np) = annNdep/365.0
!!        casaflux%Nminfix(np) = annNfix/365.0
!!        casaflux%Pdep(np)    = annPdust/365.0     ! gP/m2/day
!!        casaflux%Pwea(np)    = annPwea/365.0      ! gP/m2/day
!!!        IF (mvtype==17) THEN
!!!          vtypex(np)  = ivtigbp  ! for running IGBP veg type only
!!!        END IF
!!      END IF
!!    ELSE
!!      PRINT * ,'Check why lat, lon do not match'
!!      print * ,'casamet',casamet%lat(np),casamet%lon(np)
!!      print * ,'cable  ',patch(np)%latitude,patch(np)%longitude
!!      STOP
!!    END IF
!!
!!    if(veg%iveg(np)==9 .or. veg%iveg(np)==10) then
!!    ! P fertilizer =13 Mt P globally in 1994
!!      casaflux%Pdep(np) = casaflux%Pdep(np)+0.7/365.0
!!    ! N fertilizer =86 Mt N globally in 1994
!!      casaflux%Nminfix(np) = casaflux%Nminfix(np)+4.3/365.0
!!    endif
!!!    IF (veg%iveg(np)==12 .OR. veg%iveg(np)==14) casaflux%Pdep(np)= &
!!!       casaflux%Pdep(np)+0.7/365.0    ! P fertilizer =13 Mt P globally in 1994
!!
!!  ENDDO
!!  CLOSE(101)
!
!END SUBROUTINE casa_readpoint

SUBROUTINE casa_init(casabiome,casamet,casaflux,casapool,casabal,veg,phen)
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
  USE cable_io_vars_module, ONLY: landpt, patch
  USE cable_common_module, only: cable_user

! end addition (BP may2010)
  IMPLICIT NONE
!  INTEGER,        INTENT(IN)    :: mst
  TYPE (casa_biome),   INTENT(IN)    :: casabiome
  TYPE (casa_met),     INTENT(INOUT) :: casamet
  TYPE (casa_flux),    INTENT(INOUT) :: casaflux
  TYPE (casa_pool),    INTENT(INOUT) :: casapool
  TYPE (casa_balance), INTENT(INOUT) :: casabal
! for first time reading file *_1220.csv  (BP may2010)
  TYPE (veg_parameter_type), INTENT(IN) :: veg
  TYPE (phen_variable),   INTENT(INOUT) :: phen
  REAL(r_2) :: clabile,cplant(3),clitter(3),csoil(3)
  REAL(r_2) :: nplant(3),nlitter(3),nsoil(3),nsoilmin,pplant(3)
  REAL(r_2) :: plitter(3),psoil(3),psoillab,psoilsorb,psoilocc
! end addition (BP may2010)

  ! local variables
  INTEGER   :: np,npt,npz
  INTEGER   :: nyearz,ivtz,istz,isoz
  REAL(r_2) :: latz,lonz,areacellz,glaiz,slaz
  LOGICAL   :: EXRST

  
if (.NOT.cable_user%casa_fromzero) THEN
   PRINT *, 'initial pool from restart file'
ENDIF
  PRINT *, 'icycle,initcasa,mp ', icycle,initcasa,mp
  !phen%phase = 2

  !CLN initialise all !!!!! THIS NEEDS FIXING because of e.g. ICE-WATER
  casaflux%Cgpp         = 0.
  casaflux%Cnpp         = 0.
  casaflux%Crp          = 0.
  casaflux%Crgplant     = 0.
 ! casaflux%Nminfix      = 0.
  casaflux%Nminuptake   = 0.
  casaflux%Plabuptake   = 0.
  casaflux%Clabloss     = 0.
  casaflux%fracClabile  = 0.
  casaflux%stemnpp      = 0.
  casaflux%frac_sapwood = 0.
  casaflux%sapwood_area = 0.
  casaflux%FluxCtohwp = 0.
  casaflux%FluxCtoClear = 0.
  casaflux%fracCalloc   = 0.
  casaflux%fracNalloc   = 0.
  casaflux%fracPalloc   = 0.
  casaflux%Crmplant     = 0.
  casaflux%kplant       = 0.

  casaflux%fromPtoL     = 0.

  casaflux%Cnep         = 0.
  casaflux%Crsoil       = 0.
  casapool%dClabiledt = 0.0
  !casaflux%Nmindep      =  casaflux%Nmindep /2.0
 !casaflux%Nmindep      = 0.
  casaflux%Nminloss     = 0.
  casaflux%Nminleach    = 0.
  casaflux%Nupland      = 0.
  casaflux%Nlittermin   = 0.
  casaflux%Nsmin        = 0.
  casaflux%Nsimm        = 0.
  casaflux%Nsnet        = 0.
  !casaflux%fNminloss    = 0.
  !casaflux%fNminleach   = 0.
  !casaflux%Pdep         = 0.
  !casaflux%Pwea         = 0.
  casaflux%Pleach       = 0.
  casaflux%Ploss        = 0.
  casaflux%Pupland      = 0.
  casaflux%Plittermin   = 0.
  casaflux%Psmin        = 0.
  casaflux%Psimm        = 0.
  casaflux%Psnet        = 0.
!  casaflux%fPleach      = 0. !vh ! this should be a parameter, not a flux variable
  casaflux%kplab        = 0.
  casaflux%kpsorb       = 0.
  casaflux%kpocc        = 0.
!  casaflux%kmlabp       = 0.  !vh ! this should be a paramter, not a flux variable
!  casaflux%Psorbmax     = 0. !vh ! this should be a paramter, not a flux variable

  casaflux%klitter      = 0.
  casaflux%ksoil        = 0.
  casaflux%fromLtoS     = 0.
  casaflux%fromStoS     = 0.
  casaflux%fromLtoCO2   = 0.
  casaflux%fromStoCO2   = 0.
  casaflux%FluxCtolitter= 0.
  casaflux%FluxNtolitter= 0.
  casaflux%FluxPtolitter= 0.
  casaflux%FluxCtosoil  = 0.
  casaflux%FluxNtosoil  = 0.
  casaflux%FluxPtosoil  = 0.
  casaflux%FluxCtoCO2   = 0.

  casaflux%FluxFromPtoL = 0.
  casaflux%FluxFromLtoS = 0.
  casaflux%FluxFromStoS = 0.
  casaflux%FluxFromPtoCO2 = 0.
  casaflux%FluxFromLtoCO2 = 0.
  casaflux%FluxFromStoCO2 = 0.
  casaflux%FluxFromPtoHarvest = 0.

  casaflux%Cplant_turnover = 0.
  casaflux%fHarvest = 0.0
  casaflux%Charvest = 0.0
  casaflux%Nharvest = 0.0
  casaflux%fcrop    = 0.0

  phen%doyphase(:,1) = -50
  phen%doyphase(:,2) = phen%doyphase(:,1) +14
  phen%doyphase(:,3) = 367
  phen%doyphase(:,4) = phen%doyphase(:,3) + 14
  phen%phase(:) = 2
  phen%phen(:) = 1
  phen%aphen(:) = 0
  !CLN add more if necessary

  IF (initcasa==1) THEN
     if (.NOT.cable_user%casa_fromzero) THEN
#ifndef UM_BUILD
        CALL READ_CASA_RESTART_NC (  casamet, casapool, casaflux, phen )
#endif
     ELSE
        WRITE(*,*)'casa_init: not using restart file!'
        WRITE(*,*)'Using input from readbiome.!!!'
        WRITE(*,*) 'initialising frac_sapwood=1 and sapwood_area = 0)'
        casaflux%frac_sapwood(:) = 1.0
        casaflux%sapwood_area(:) = 0.0
     ENDIF
  ENDIF
 WHERE(casamet%lnonwood==1) casapool%cplant(:,WOOD) = 0.0
 WHERE(casamet%lnonwood==1) casapool%nplant(:,WOOD) = 0.0
 WHERE(casamet%lnonwood==1) casapool%pplant(:,WOOD) = 0.0
!!$IF (initcasa==1) THEN
!!$     INQUIRE( FILE=TRIM(casafile%cnpipool), EXIST=EXRST )
!!$!! vh_js!!
!!$     IF ( EXRST ) THEN
!!$
!!$           PRINT*, ' Reading cnppoolOutfile as input: ,',casafile%cnpipool
!!$
!!$    OPEN(99,file=casafile%cnpipool)
!!$
!!$    DO npt =1, mp
!!$       SELECT CASE(icycle)
!!$       CASE(1)
!!$          !! vh_js !!
!!$          IF (cable_user%CALL_POP) THEN
!!$
!!$             READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!!$                  casamet%glai(npt),slaz,phen%phase(npt) , &
!!$                  phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
!!$                  casapool%clabile(npt) ,casapool%cplant(npt,:) ,  &
!!$                  casapool%clitter(npt,:),casapool%csoil(npt,:), &
!!$                  casaflux%frac_sapwood(npt), casaflux%sapwood_area(npt)
!!$
!!$
!!$             ELSE
!!$              READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!!$                  casamet%glai(npt),slaz,phen%phase(npt) , &
!!$                  phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
!!$                  casapool%clabile(npt) ,casapool%cplant(npt,:) ,  &
!!$                  casapool%clitter(npt,:),casapool%csoil(npt,:)
!!$             casaflux%frac_sapwood(:) = 1.0
!!$             casaflux%sapwood_area(:) = 0.0
!!$          ENDIF
!!$
!!$
!!$       CASE(2)
!!$!! vh_js !!
!!$          IF (cable_user%CALL_POP) THEN
!!$             READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!!$                  casamet%glai(npt),slaz,phen%phase(npt), &
!!$                  phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
!!$                  casapool%clabile(npt),casapool%cplant(npt,:),   &
!!$                  casapool%clitter(npt,:),casapool%csoil(npt,:),       &
!!$                  casaflux%frac_sapwood(npt), casaflux%sapwood_area(npt), &
!!$                  casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
!!$                  casapool%nsoil(npt,:),casapool%nsoilmin(npt)
!!$
!!$          ELSE
!!$             READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!!$                  casamet%glai(npt),slaz,phen%phase(npt), &
!!$                  phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
!!$                  casapool%clabile(npt),casapool%cplant(npt,:),   &
!!$                  casapool%clitter(npt,:),casapool%csoil(npt,:),       &
!!$                  casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
!!$                  casapool%nsoil(npt,:),casapool%nsoilmin(npt)
!!$             casaflux%frac_sapwood(:) = 1.0
!!$             casaflux%sapwood_area(:) = 0.0
!!$
!!$          ENDIF
!!$       CASE(3)
!!$!! vh_js !!
!!$          IF (cable_user%CALL_POP) THEN
!!$             READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!!$                  casamet%glai(npt),slaz,phen%phase(npt), &
!!$                  phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
!!$                  casapool%clabile(npt),casapool%cplant(npt,:),   &
!!$                  casapool%clitter(npt,:),casapool%csoil(npt,:),       &
!!$                  casaflux%frac_sapwood(npt), casaflux%sapwood_area(npt), &
!!$                  casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
!!$                  casapool%nsoil(npt,:),casapool%nsoilmin(npt),        &
!!$                  casapool%pplant(npt,:),casapool%plitter(npt,:),      &
!!$                  casapool%psoil(npt,:),casapool%psoillab(npt),        &
!!$                  casapool%psoilsorb(npt),casapool%psoilocc(npt)
!!$          ELSE
!!$             READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!!$                  casamet%glai(npt),slaz,phen%phase(npt), &
!!$                  phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
!!$                  casapool%clabile(npt),casapool%cplant(npt,:),   &
!!$                  casapool%clitter(npt,:),casapool%csoil(npt,:),       &
!!$                  casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
!!$                  casapool%nsoil(npt,:),casapool%nsoilmin(npt),        &
!!$                  casapool%pplant(npt,:),casapool%plitter(npt,:),      &
!!$                  casapool%psoil(npt,:),casapool%psoillab(npt),        &
!!$                  casapool%psoilsorb(npt),casapool%psoilocc(npt)
!!$             casaflux%frac_sapwood(:) = 1.0
!!$             casaflux%sapwood_area(:) = 0.0
!!$
!!$
!!$          ENDIF
!!$       END SELECT
!!$       IF (ABS(patch(npt)%longitude - lonz) > 0.9 .OR. &
!!$            ABS(patch(npt)%latitude  - latz) > 0.9) THEN
!!$          PRINT *, 'patch(npt)%longitude, lonz:', patch(npt)%longitude, lonz
!!$          PRINT *, 'patch(npt)%latitude,  latz:', patch(npt)%latitude,  latz
!!$          PRINT *, 'npt = ', npt
!!$          STOP
!!$       ENDIF
!!$    ENDDO
!!$    CLOSE(99)
!!$
!!$
!!$ ELSE
!!$ !! vh_js !!
!!$    WRITE(*,*)'No valid restart file for casa_init found.'
!!$    WRITE(*,*)'Using input from readbiome.!!!'
!!$    WRITE(*,*) 'initialising frac_sapwood=1 and sapwood_area = 0)'
!!$    casaflux%frac_sapwood(:) = 1.0
!!$    casaflux%sapwood_area(:) = 0.0
!!$
!!$
!!$ ENDIF  ! IF (EXRST)

!!$ENDIF
!92 format(5(i6,2x),5(f18.6,3x),2(i6,',',2x),',',2x,100(f18.6,3x))
92    format(5(i6,',',2x),5(f18.6,',',2x),2(i6,',',2x),',',2x,100(f18.6,',',2x))


if(initcasa==0) then
   nyearz = 1
   do npt=1,mp
      casamet%lon(npt) = patch(npt)%longitude
      casamet%lat(npt) = patch(npt)%latitude
   enddo
endif

  ! reset labile C pool,comment out by Q.Zhang 10/09/2011
  !  casapool%clabile    = 0.0
  ! check pool sizes
  casapool%Ctot_0 = 0.0
  casapool%Ctot = 0.0
  casapool%cplant     = MAX(0.0,casapool%cplant)
  casapool%clitter    = MAX(0.0,casapool%clitter)
  casapool%csoil      = MAX(0.0,casapool%csoil)
  casabal%cplantlast  = casapool%cplant
  casabal%clitterlast = casapool%clitter
  casabal%csoillast   = casapool%csoil
  casabal%clabilelast = casapool%clabile
  casabal%sumcbal     = 0.0
  casabal%FCgppyear=0.0;casabal%FCrpyear=0.0
  casabal%FCnppyear=0;casabal%FCrsyear=0.0;casabal%FCneeyear=0.0
  !vh !
  WHERE(casamet%lnonwood==1) casapool%cplant(:,WOOD) = 0.0
  IF (icycle==1) THEN
    casapool%Nplant(:,:) = casapool%cplant(:,:) * casapool%ratioNCplant(:,:)
    casapool%Nsoil(:,:)  = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
    casapool%Psoil(:,:)  = casapool%Nsoil(:,:)/casapool%ratioNPsoil(:,:)
    casapool%Nsoilmin(:) = 2.5
  ENDIF

  IF (icycle >=1) THEN
     casapool%nplant     = MAX(1.e-6,casapool%nplant)
     casapool%nlitter    = MAX(1.e-6,casapool%nlitter)
     casapool%nsoil      = MAX(1.e-6,casapool%nsoil)
     casapool%nsoilmin   = MAX(1.e-6,casapool%nsoilmin)
     casabal%nplantlast  = casapool%nplant
     casabal%nlitterlast = casapool%nlitter
     casabal%nsoillast   = casapool%nsoil
     casabal%nsoilminlast= casapool%nsoilmin
     casabal%sumnbal     = 0.0
     casabal%FNdepyear=0.0;casabal%FNfixyear=0.0;casabal%FNsnetyear=0.0
     casabal%FNupyear=0.0;casabal%FNleachyear=0.0;casabal%FNlossyear=0.0
     !vh !
     WHERE(casamet%lnonwood==1) casapool%nplant(:,WOOD) = 0.0
  ENDIF

  IF (icycle >=1) THEN
     casapool%pplant       = MAX(1.0e-7,casapool%pplant)
     casapool%plitter      = MAX(1.0e-7,casapool%plitter)
     casapool%psoil        = MAX(1.0e-7,casapool%psoil)
     casapool%Psoillab     = MAX(1.0e-7,casapool%psoillab)  ! was 2.0, changed according to  YP
     casapool%psoilsorb    = MAX(1.0e-7,casapool%psoilsorb) ! was 10.0, -
     casapool%psoilocc     = MAX(1.0e-7,casapool%psoilocc)  ! was 50.0, -
     casabal%pplantlast    = casapool%pplant
     casabal%plitterlast   = casapool%plitter
     casabal%psoillast     = casapool%psoil
     casabal%psoillablast  = casapool%psoillab
     casabal%psoilsorblast = casapool%psoilsorb
     casabal%psoilocclast  = casapool%psoilocc
     casabal%sumpbal       = 0.0
     casabal%FPweayear=0.0;casabal%FPdustyear=0.0; casabal%FPsnetyear=0.0
     casabal%FPupyear=0.0;casabal%FPleachyear=0.0;casabal%FPlossyear=0.0
     !vh !
     WHERE(casamet%lnonwood==1) casapool%pplant(:,WOOD) = 0.0
  EndIF


END SUBROUTINE casa_init


SUBROUTINE casa_poolout(ktau,veg,soil,casabiome,casapool,casaflux,casamet, &
                        casabal,phen)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  USE cable_common_module, only: cable_user
  IMPLICIT NONE
  INTEGER,               INTENT(IN)    :: ktau
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (casa_balance),        INTENT(INOUT) :: casabal
  TYPE (phen_variable),       INTENT(INOUT) :: phen

  ! local variables
  REAL(r_2), DIMENSION(mso) :: Psorder,pweasoil,xpsoil50
  REAL(r_2), DIMENSION(mso) :: fracPlab,fracPsorb,fracPocc,fracPorg
  REAL(r_2), DIMENSION(mp)  :: totpsoil
  INTEGER  npt,nout,nso

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
  DATA psorder/61.3,103.9,92.8,136.9,98.2,107.6,84.1,110.1,35.4,41.0,51.5,190.6/
  DATA pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/
  DATA fracpLab/0.08,0.08,0.10,0.02,0.08,0.08,0.08,0.06,0.02,0.05,0.09,0.05/
  DATA fracPsorb/0.32,0.37,0.57,0.67,0.37,0.37,0.37,0.32,0.24,0.22,0.21,0.38/
  DATA fracPocc/0.36,0.38,0.25,0.26,0.38,0.38,0.38,0.44,0.38,0.38,0.37,0.45/
  DATA fracPorg/0.25,0.17,0.08,0.05,0.17,0.17,0.17,0.18,0.36,0.35,0.34,0.12/
  DATA xpsoil50/7.6,4.1,4.2,3.4,4.1,4.1,4.8,4.1,6.9,6.9,6.9,1.7/
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

  PRINT *, 'Within casa_poolout, mp = ', mp
  nout=103
  OPEN(nout,file=casafile%cnpepool)
  PRINT *, 'Opened file ', casafile%cnpepool

  casabal%sumcbal=MIN(9999.0,MAX(-9999.0,casabal%sumcbal))
  casabal%sumnbal=MIN(9999.0,MAX(-9999.0,casabal%sumnbal))
  casabal%sumpbal=MIN(9999.0,MAX(-9999.0,casabal%sumpbal))

  DO npt =1, mp
    nso = casamet%isorder(npt)
    totpsoil(npt) = psorder(nso) *xpsoil50(nso)
  if(casamet%iveg2(npt)>0 ) then
    IF (icycle<2) THEN
      casapool%Nplant(npt,:) = casapool%ratioNCplant(npt,:)  &
                             * casapool%cplant(npt,:)
      casapool%Nlitter(npt,:)= casapool%ratioNClitter(npt,:) &
                             * casapool%clitter(npt,:)
      casapool%Nsoil(npt,:)  = casapool%ratioNCsoil(npt,:)   &
                             * casapool%Csoil(npt,:)
      casapool%nsoilmin(npt) = 2.0
      casabal%sumnbal(npt)   = 0.0
      if(casamet%iveg2(npt)==grass) then
         casapool%nplant(npt,wood) = 0.0
         casapool%nlitter(npt,cwd) = 0.0
      endif
    ENDIF

    IF (icycle<3) THEN
      casabal%sumpbal(npt)   = 0.0
      casapool%pplant(npt,:)  = casapool%Nplant(npt,:)/casapool%ratioNPplant(npt,:)
      casapool%plitter(npt,:) = casapool%Nlitter(npt,:)/(casapool%ratioNPlitter(npt,:)+1.0e-10)
      casapool%psoil(npt,:)   = casapool%Nsoil(npt,:)/casapool%ratioNPsoil(npt,:)
      casapool%psoillab(npt) = totpsoil(npt) *fracpLab(nso)
      casapool%psoilsorb(npt)= casaflux%psorbmax(npt) * casapool%psoillab(npt) &
                                /(casaflux%kmlabp(npt)+casapool%psoillab(npt))
      casapool%psoilocc(npt) = totpsoil(npt) *fracPocc(nso)
      if(casamet%iveg2(npt)==grass) then
         casapool%pplant(npt,wood) = 0.0
         casapool%plitter(npt,cwd) = 0.0
      endif
    ENDIF
  else
     casapool%cplant(npt,:)=0.0; casapool%clitter(npt,:)=0.0; casapool%csoil(npt,:) = 0.0; casapool%clabile(npt) = 0.0
     casapool%nplant(npt,:)=0.0; casapool%nlitter(npt,:)=0.0; casapool%nsoil(npt,:) = 0.0; casapool%nsoilmin(npt) = 0.0
     casapool%pplant(npt,:)=0.0; casapool%plitter(npt,:)=0.0; casapool%psoil(npt,:) = 0.0
     casapool%psoillab(npt) = 0.0; casapool%psoilsorb(npt) = 0.0; casapool%psoilocc(npt) = 0.0
     casabal%sumcbal(npt) =0.0; casabal%sumnbal(npt) =0.0; casabal%sumpbal(npt) = 0.0
  endif

!! vh_js  !! 
  IF (cable_user%CALL_POP) THEN
   
     WRITE(nout,92) ktau,npt,veg%iveg(npt),soil%isoilm(npt) ,     &
          casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
         casamet%areacell(npt)*(1.0e-9),casamet%glai(npt),       &
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
          casamet%areacell(npt)*(1.0e-9),casamet%glai(npt),       &
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

92    format(5(i6,',',2x),5(f18.6,',',2x),2(i6,',',2x),100(f18.6,',',2x))
END SUBROUTINE casa_poolout

! casa_fluxout output data for Julie Tang; comment out (BP apr2010)
SUBROUTINE casa_fluxout(myear,veg,soil,casabal,casamet)
!SUBROUTINE casa_fluxout(myear,clitterinput,csoilinput)
  USE cable_def_types_mod
!  USE cableDeclare, ONLY: veg, soil
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
!  USE casaDeclare
  IMPLICIT NONE
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (casa_balance),        INTENT(INOUT) :: casabal
  INTEGER,               INTENT(IN)    :: myear
!  REAL(r_2),    INTENT(IN) :: clitterinput(mp,3),csoilinput(mp,3)

  ! local variables
  INTEGER  npt,nout
  REAL(r_2) xyear, totGPP, totNPP

  totGPP =0.0
  totNPP =0.0
  nout=104
  xyear=1.0/FLOAT(myear)
  casabal%FCgppyear=casabal%FCgppyear * xyear
  casabal%FCnppyear=casabal%FCnppyear * xyear
  casabal%FCrmleafyear=casabal%FCrmleafyear * xyear
  casabal%FCrmwoodyear=casabal%FCrmwoodyear * xyear
  casabal%FCrmrootyear=casabal%FCrmrootyear * xyear
  casabal%FCrgrowyear=casabal%FCrgrowyear * xyear
  casabal%FCrsyear=casabal%FCrsyear * xyear
  casabal%FCneeyear=casabal%FCneeyear * xyear
  casabal%FNdepyear=casabal%FNdepyear * xyear
  casabal%FNfixyear=casabal%FNfixyear * xyear
  casabal%FNsnetyear=casabal%FNsnetyear * xyear
  casabal%FNupyear=casabal%FNupyear * xyear
  casabal%FNleachyear=casabal%FNleachyear * xyear
  casabal%FNlossyear=casabal%FNlossyear * xyear
  casabal%FPweayear=casabal%FPweayear * xyear
  casabal%FPdustyear=casabal%FPdustyear * xyear
  casabal%FPsnetyear=casabal%FPsnetyear * xyear
  casabal%FPupyear=casabal%FPupyear * xyear
  casabal%FPleachyear=casabal%FPleachyear * xyear
  casabal%FPlossyear=casabal%FPlossyear * xyear
!  clitterinput = clitterinput * xyear
!  csoilinput   = csoilinput   * xyear

  print *, 'writing CNP fluxes out to file ', casafile%cnpflux
  OPEN(nout,file=casafile%cnpflux)
    DO npt =1,mp
      SELECT CASE(icycle)
      CASE(1)

        WRITE(nout,*) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
            casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
            casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt),  &
            casabal%Fcnppyear(npt),  &
            casabal%Fcrmleafyear(npt),casabal%Fcrmwoodyear(npt),     &
            casabal%Fcrmrootyear(npt),casabal%Fcrgrowyear(npt),     &
            casabal%Fcrsyear(npt),casabal%Fcneeyear(npt)  ! ,           &
!            clitterinput(npt,:),csoilinput(npt,:)

      CASE(2)
        WRITE(nout,*) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
            casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
            casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt),  &
            casabal%FCnppyear(npt),                                 &
            casabal%Fcrmleafyear(npt),casabal%Fcrmwoodyear(npt),     &
            casabal%Fcrmrootyear(npt),casabal%Fcrgrowyear(npt),     &
            casabal%FCrsyear(npt), casabal%FCneeyear(npt),          &
!        clitterinput(npt,:),csoilinput(npt,:), &
        casabal%FNdepyear(npt),casabal%FNfixyear(npt),casabal%FNsnetyear(npt), &
        casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt)

      CASE(3)
        WRITE(nout,*) myear,npt,veg%iveg(npt),soil%isoilm(npt), &
        casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt),  &
        casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt), &
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

    print *, 'totGPP global = ', totGPP*(1.0e-15)
    print *, 'totNPP global = ', totNPP*(1.0e-15)
  CLOSE(nout)
92    format(5(i6,',',2x),100(f15.6,',',2x))
END SUBROUTINE casa_fluxout

! clitterinput and csoilinput are for Julie Tang; comment out (BP apr2010)
!SUBROUTINE casa_cnpflux(clitterinput,csoilinput)
SUBROUTINE casa_cnpflux(casaflux,casapool,casabal,zeroflux)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  IMPLICIT NONE
  TYPE (casa_flux),    INTENT(INOUT) :: casaflux
  TYPE (casa_pool),    INTENT(IN)    :: casapool
  TYPE (casa_balance), INTENT(INOUT) :: casabal
  LOGICAL :: zeroflux
  !  REAL(r_2), INTENT(INOUT) :: clitterinput(mp,3),csoilinput(mp,3)
  INTEGER n

  IF(zeroflux) THEN
     casabal%FCgppyear    = 0.0
     casabal%FCrpyear     = 0.0   
     casabal%FCrmleafyear = 0.0
     casabal%FCrmwoodyear = 0.0
     casabal%FCrmrootyear = 0.0
     casabal%FCrgrowyear  = 0.0
     casabal%FCnppyear    = 0.0
     casabal%FCrsyear     = 0.0
     casabal%FCneeyear    = 0.0
     casabal%dCdtyear    = 0.0
    

     casabal%FNdepyear    = 0.0
     casabal%FNfixyear    = 0.0
     casabal%FNsnetyear   = 0.0
     casabal%FNupyear     = 0.0
     casabal%FNleachyear  = 0.0
     casabal%FNlossyear   = 0.0

     casabal%FPweayear   = 0.0
     casabal%FPdustyear  = 0.0
     casabal%FPsnetyear  = 0.0
     casabal%FPupyear    = 0.0
     casabal%FPleachyear = 0.0
     casabal%FPlossyear  = 0.0

     casaflux%FluxCtohwp = 0.0
     casaflux%FluxNtohwp = 0.0
     casaflux%FluxPtohwp = 0.0
     casaflux%FluxCtoclear = 0.0
     casaflux%FluxNtoclear = 0.0
     casaflux%FluxPtoclear = 0.0
     casaflux%CtransferLUC = 0.0

     casaflux%Cplant_turnover_disturbance = 0
     casaflux%Cplant_turnover_crowding = 0
     casaflux%Cplant_turnover_resource_limitation = 0

  ELSE

     casabal%FCgppyear = casabal%FCgppyear + casaflux%Cgpp   * deltpool
     casabal%FCrpyear  = casabal%FCrpyear  + casaflux%Crp    * deltpool
     casabal%FCrmleafyear(:)  = casabal%FCrmleafyear(:)  + casaflux%Crmplant(:,leaf)    * deltpool
     casabal%FCrmwoodyear(:)  = casabal%FCrmwoodyear(:)  + casaflux%Crmplant(:,wood)    * deltpool
     casabal%FCrmrootyear(:)  = casabal%FCrmrootyear(:)  + casaflux%Crmplant(:,froot)    * deltpool
     casabal%FCrgrowyear      = casabal%FCrgrowyear  + casaflux%Crgplant              * deltpool
     ! change made ypwang 17-nov-2013 to accoutn for change in labile carbon pool  size
     casabal%FCnppyear        = casabal%FCnppyear + (casaflux%Cnpp+casapool%dClabiledt)   * deltpool
     casabal%FCrsyear  = casabal%FCrsyear  + casaflux%Crsoil * deltpool
     !casabal%FCneeyear = casabal%FCneeyear &
     !     + (casaflux%Cnpp+casapool%dClabiledt-casaflux%Crsoil) * deltpool
     casabal%FCneeyear = casabal%FCneeyear &
          + (casaflux%Cnpp-casaflux%Crsoil) * deltpool
     casabal%dCdtyear =  casabal%dCdtyear + (casapool%Ctot-casapool%Ctot_0)*deltpool
   
     !  DO n=1,3
     !    clitterinput(:,n)= clitterinput(:,n) + casaflux%kplant(:,n) * casapool%cplant(:,n) * deltpool
     !    csoilinput(:,n) = csoilinput(:,n) + casaflux%fluxCtosoil(:,n) * deltpool
     !    !csoilinput(:,n) = csoilinput(:,n)+casaflux%fluxCtolitter(:,n)*deltpool
     !  ENDDO

     IF (icycle >1) THEN
        casabal%FNdepyear   = casabal%FNdepyear   + casaflux%Nmindep    * deltpool
        casabal%FNfixyear   = casabal%FNfixyear   + casaflux%Nminfix    * deltpool
        casabal%FNsnetyear  = casabal%FNsnetyear  + casaflux%Nsnet      * deltpool
        casabal%FNupyear    = casabal%FNupyear    + casaflux%Nminuptake * deltpool
        casabal%FNleachyear = casabal%FNleachyear + casaflux%Nminleach  * deltpool
        casabal%FNlossyear  = casabal%FNlossyear  + casaflux%Nminloss   * deltpool
     ENDIF

     IF (icycle >2) THEN
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
SUBROUTINE biogeochem(ktau,dels,idoY,LALLOC,veg,soil,casabiome,casapool,casaflux, &
     casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf,xkleafcold,xkleafdry,&
     cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
     nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
     pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
  
  USE cable_def_types_mod
  USE casadimension
  USE casa_cnp_module
  USE POP_TYPES,            ONLY: POP_TYPE
  USE cable_IO_vars_module, ONLY: wlogn
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)    :: ktau
  REAL,    INTENT(IN)    :: dels
  INTEGER, INTENT(IN)    :: idoy
  INTEGER, INTENT(IN)    :: LALLOC
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_balance),          INTENT(INOUT) :: casabal
  TYPE (phen_variable),         INTENT(INOUT) :: phen
  TYPE(POP_TYPE),               INTENT(IN)    :: POP
  TYPE(climate_TYPE),           INTENT(IN)    :: climate

  ! local variables added by ypwang following Chris Lu 5/nov/2012

  real, dimension(mp), INTENT(OUT)   :: cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
       nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
       pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd

  ! local variables
  REAL(r_2),    DIMENSION(mp) :: xnplimit,xNPuptake
  REAL(r_2),    DIMENSION(mp) :: xklitter,xksoil,xkNlimiting
  REAL(r_2),    DIMENSION(mp) :: xkleafcold,xkleafdry,xkleaf
  INTEGER  npt,j
  REAL, ALLOCATABLE :: tmp(:)

  xKNlimiting = 1.0

 ! zero annual sums
  if (idoy==1) CALL casa_cnpflux(casaflux,casapool,casabal,.true.)

  IF (cable_user%PHENOLOGY_SWITCH.eq.'MODIS') THEN
     call phenology(idoy,veg,phen)
  ENDIF
  call avgsoil(veg,soil,casamet)
  call casa_rplant(veg,casabiome,casapool,casaflux,casamet,climate)


   IF (.NOT.cable_user%CALL_POP) THEN
      call casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen,LALLOC)
   ENDIF

   call casa_xrateplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome, &
        casamet,phen)
   call casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
        casaflux,casamet,phen)

   call casa_xnp(xnplimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)

   IF (cable_user%CALL_POP) THEN

      call casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen,LALLOC)
      WHERE (pop%pop_grid(:)%cmass_sum_old.gt.0.001 .and. pop%pop_grid(:)%cmass_sum.gt.0.001 )
         
         casaflux%frac_sapwood(POP%Iwood) = POP%pop_grid(:)%csapwood_sum/ POP%pop_grid(:)%cmass_sum
         casaflux%sapwood_area(POP%Iwood) = max(POP%pop_grid(:)%sapwood_area/10000., 1e-6)
         veg%hc(POP%Iwood) = POP%pop_grid(:)%height_max
          
         WHERE (pop%pop_grid(:)%LU ==2)

            casaflux%kplant(POP%Iwood,2) =  1.0 -  &
              (1.0-  max( min((POP%pop_grid(:)%stress_mortality + &
              POP%pop_grid(:)%crowding_mortality + &
              POP%pop_grid(:)%fire_mortality ) &
              /(POP%pop_grid(:)%cmass_sum+POP%pop_grid(:)%growth) + &
              1.0/veg%disturbance_interval(POP%Iwood,1), 0.99), 0.0))**(1.0/365.0)

        ELSEWHERE
            casaflux%kplant(POP%Iwood,2) =  1.0 -  &
              (1.0-  max( min((POP%pop_grid(:)%stress_mortality + &
              POP%pop_grid(:)%crowding_mortality + &
              POP%pop_grid(:)%fire_mortality+POP%pop_grid(:)%cat_mortality  ) &
              /(POP%pop_grid(:)%cmass_sum+POP%pop_grid(:)%growth), 0.99), 0.0))**(1.0/365.0)

         ENDWHERE

         veg%hc(POP%Iwood) = POP%pop_grid(:)%height_max
      ELSEWHERE
         casaflux%frac_sapwood(POP%Iwood) = 1.0
         casaflux%sapwood_area(POP%Iwood) = max(POP%pop_grid(:)%sapwood_area/10000., 1e-6)
         casaflux%kplant(POP%Iwood,2) = 0.0
         veg%hc(POP%Iwood) = POP%pop_grid(:)%height_max
      ENDWHERE

   ENDIF
!!$if (idoy.eq.365) then
!!$ write(667,*) pop%LU
!!$ write(667,*) veg%ilu
!!$ write(667,991) casaflux%FluxCtohwp(POP%Iwood,1)
!!$ write(667,991) POP%pop_grid(:)%cat_mortality/POP%pop_grid(:)%cmass_sum_old
!!$   write(667,991)max(min((POP%pop_grid(:)%cat_mortality                &
!!$        /POP%pop_grid(:)%cmass_sum_old),0.99),0.0)**(1.0/365.0)
!!$   write(667,991) (1.0 - (1.0 -max( min((POP%pop_grid(:)%cat_mortality  &
!!$        /POP%pop_grid(:)%cmass_sum_old),0.99), 0.0))**(1.0/365.0))
!!$write(667,*)
!!$   endif
!!$  
!write(667,991) casaflux%cgpp(147),casaflux%cnpp(147),casaflux%kplant(147,2),casapool%cplant(147,:)
!  write(*,991)casaflux%cgpp(2058),casaflux%cnpp(2058),casaflux%fracClabile(2058), &
!            casaflux%fracCalloc(2058,:),casaflux%crmplant(2058,:),casaflux%crgplant(2058), casapool%Nsoilmin(2058), &
!            casaflux%cgpp(2058)-casaflux%cnpp(2058)-casaflux%fracClabile(2058)*casaflux%cgpp(2058)-sum(casaflux%crmplant(2058,:))-casaflux%crgplant(2058)
   !991  format('point 147',20(f10.4,2x))
   991  format(20(e12.4,2x))

  call casa_xratesoil(xklitter,xksoil,veg,soil,casamet,casabiome)
  call casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casaflux,casamet)

  IF (icycle>1) THEN
    call casa_xkN(xkNlimiting,casapool,casaflux,casamet,casabiome,veg)
   ! vh ! why is this here and not in casa_cnp.F90?

    DO j=1,mlitter
      casaflux%klitter(:,j) = casaflux%klitter(:,j)* xkNlimiting(:)
    ENDDO
    call casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
    IF (icycle >2) call casa_puptake(veg,xkNlimiting,casabiome, &
                                     casapool,casaflux,casamet)
  ENDIF

 !write(62,"(100e16.6)") xkNlimiting
!!$
!!$ ! reset base soil turnover rates for grass tiles to be the same 
!!$  ! as potential veg in the same gridcell
!!$  IF (cable_user%CALL_POP) THEN
!!$     DO j=1,msoil
!!$        where (veg%ilu ==3)
!!$           casaflux%ksoil(:,j) = casaflux%ksoil(:,j) * casabiome%soilrate(veg%ivegp(:),j)/ &
!!$                casabiome%soilrate(veg%iveg(:),j)
!!$        ENDWHERE
!!$     ENDDO
!!$  ENDIF

  ! changed by ypwang following Chris Lu on 5/nov/2012
  call casa_delplant(veg,casabiome,casapool,casaflux,casamet,                &
       cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
       nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
       pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

  casaflux%Cplant_turnover_disturbance = 0
  casaflux%Cplant_turnover_crowding = 0
  casaflux%Cplant_turnover_resource_limitation = 0

  if (cable_user%CALL_POP) THEN
     if (.not.allocated(tmp)) allocate(tmp(size(POP%pop_grid)))
     tmp = (POP%pop_grid(:)%stress_mortality + POP%pop_grid(:)%crowding_mortality &
          +POP%pop_grid(:)%cat_mortality &
          + POP%pop_grid(:)%fire_mortality  )
     where (tmp.gt. 1.e-12)
        casaflux%Cplant_turnover_disturbance(POP%Iwood) =  &
             casaflux%Cplant_turnover(POP%Iwood,2)*(POP%pop_grid(:)%cat_mortality &
             + POP%pop_grid(:)%fire_mortality  )/tmp
        casaflux%Cplant_turnover_crowding(POP%Iwood) =  &
             casaflux%Cplant_turnover(POP%Iwood,2)*POP%pop_grid(:)%crowding_mortality/tmp
        casaflux%Cplant_turnover_resource_limitation(POP%Iwood) = &
             casaflux%Cplant_turnover(POP%Iwood,2)*POP%pop_grid(:)%stress_mortality/tmp
     endwhere
  endif

  call casa_delsoil(veg,casapool,casaflux,casamet,casabiome)

  call casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet, LALLOC)
   !! vh_js !!
  !CLN ndummy must be before pdummy!!!!
  IF (icycle<3) then
      IF (icycle<2) call casa_ndummy(casapool)
      call casa_pdummy(casapool)
  ENDIF

  call casa_cnpbal(casapool,casaflux,casabal)

  call casa_cnpflux(casaflux,casapool,casabal,.false.)

  ! for spinning up only
  IF (cable_user%limit_labile) THEN
     casapool%Nsoilmin = max(casapool%Nsoilmin,5.0)
     casapool%Psoillab = max(casapool%Psoillab,1.0)
  ENDIF




END SUBROUTINE biogeochem

#ifndef UM_BUILD
SUBROUTINE WRITE_CASA_RESTART_NC ( casamet, casapool, casaflux, phen, CASAONLY )

  USE CASAVARIABLE, ONLY : casa_met, casa_pool, casa_flux, icycle, mplant, mlitter, msoil
  USE CABLE_COMMON_MODULE
  USE CABLE_DEF_TYPES_MOD, ONLY: MET_TYPE, mp
  USE phenvariable
  USE casavariable
  USE netcdf

  IMPLICIT NONE

 
  TYPE (casa_met),  INTENT(IN) :: casamet
  TYPE (casa_pool),  INTENT(IN) :: casapool
  TYPE (casa_flux),           INTENT(IN) :: casaflux
  TYPE (phen_variable),       INTENT(IN) :: phen

  INTEGER*4 :: mp4
  INTEGER*4, parameter   :: pmp4 =0
  INTEGER, parameter   :: fmp4 = kind(pmp4)
  INTEGER*4   :: STATUS
  INTEGER*4   :: FILE_ID, land_ID, plnt_ID, litt_ID, soil_ID, i
  LOGICAL   :: CASAONLY
  CHARACTER :: CYEAR*4, FNAME*99,dum*50

  ! ! 1 dim arrays (npt )
  ! CHARACTER(len=20),DIMENSION(7), PARAMETER :: A1 = (/ 'latitude', 'longitude', 'glai', &
  !      'clabile', 'psoillab','psoilsorb','psoilocc' /)
  ! ! 2 dim arrays (npt,mplant)
  ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A2 = (/ 'cplant' , 'nplant' , 'pplantc' /)
  ! ! 2 dim arrays (npt,mlitter)
  ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A3 = (/ 'clitter', 'nlitter', 'plitter' /)
  ! ! 2 dim arrays (npt,msoil)
  ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A4 = (/ 'csoil', 'nsoil', 'psoil' /)

  ! 1 dim arrays (npt )
  CHARACTER(len=20),DIMENSION(14) :: A1
  CHARACTER(len=20),DIMENSION(2) :: AI1
  ! 2 dim arrays (npt,mplant)
  CHARACTER(len=20),DIMENSION(3) :: A2
  ! 2 dim arrays (npt,mlitter)
  CHARACTER(len=20),DIMENSION(3) :: A3
  ! 2 dim arrays (npt,msoil)
  CHARACTER(len=20),DIMENSION(3) :: A4
  INTEGER*4 :: VID1(SIZE(A1)), VIDI1(SIZE(AI1)), VID2(SIZE(A2)), VID3(SIZE(A3)), VID4(SIZE(A4))

  mp4=int(mp,fmp4)
  A1(1) = 'latitude'
  A1(2) = 'longitude'
  A1(3) = 'glai'
  A1(4) = 'clabile'
  A1(5) = 'psoillab'
  A1(6) = 'psoilsorb'
  A1(7) = 'psoilocc'
  A1(8) = 'frac_sapwood'
  A1(9) = 'sapwood_area'
  A1(10) = 'phen'
  A1(11) = 'aphen'
  A1(12) = 'nsoilmin'
  A1(13) = 'fHarvest'
  A1(14) = 'fCrop'

  AI1(1) = 'phase'
  AI1(2) = 'doyphase3'


  A2(1) = 'cplant'
  A2(2) = 'nplant'
  A2(3) = 'pplant'
  A3(1) = 'clitter'
  A3(2) = 'nlitter'
  A3(3) = 'plitter'
  A4(1) = 'csoil'
  A4(2) = 'nsoil'
  A4(3) = 'psoil'

  ! Get File-Name
  WRITE(CYEAR, FMT='(I4)') CurYear + 1

    IF (( LEN_TRIM(casafile%cnpepool) ) .gt. 0) THEN
        fname=TRIM(casafile%cnpepool)
    ELSE
        fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
            '_casa_rst.nc'
    ENDIF

  ! Create NetCDF file:
  STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
write(*,*) 'writing casa restart', fname 
  ! Put the file in define mode:
  STATUS = NF90_redef(FILE_ID)

  STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Valid restart date", "01/01/"//CYEAR  )
  STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Icycle", icycle  )
  IF ( CASAONLY ) THEN
     dum = 'CASA-ONLY run'
  ELSE
     dum = 'CABLE-CASA coupled run'
  ENDIF
  STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Run-Type", TRIM(dum) )

  ! Define dimensions:
  ! Land (number of points)
  STATUS = NF90_def_dim(FILE_ID, 'land'   , mp4     , land_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_def_dim(FILE_ID, 'mplant' , mplant , plnt_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_def_dim(FILE_ID, 'mlitter', mlitter, litt_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_def_dim(FILE_ID, 'msoil'  , msoil  , soil_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  DO i = 1, SIZE(A1)
     STATUS = NF90_def_var(FILE_ID,TRIM(A1(i)) ,NF90_FLOAT,(/land_ID/),VID1(i))
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  END DO

  DO i = 1, SIZE(AI1)
     STATUS = NF90_def_var(FILE_ID,TRIM(AI1(i)) ,NF90_INT,(/land_ID/),VIDI1(i))
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  END DO

  DO i = 1, SIZE(A2)
     STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/land_ID,plnt_ID/),VID2(i))
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  END DO

  DO i = 1, SIZE(A3)
     STATUS = NF90_def_var(FILE_ID,TRIM(A3(i)) ,NF90_FLOAT,(/land_ID,litt_ID/),VID3(i))
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  END DO

  DO i = 1, SIZE(A4)
     STATUS = NF90_def_var(FILE_ID,TRIM(A4(i)) ,NF90_FLOAT,(/land_ID,soil_ID/),VID4(i))
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  END DO

  ! End define mode:
  STATUS = NF90_enddef(FILE_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  ! PUT LAT / LON
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(1), casamet%lat )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(2), casamet%lon )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT VARS
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(3), casamet%glai )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(4), casapool%clabile )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


  STATUS = NF90_PUT_VAR(FILE_ID, VID1(8), casaflux%frac_sapwood )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(9), casaflux%sapwood_area )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), phen%phen )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), phen%aphen )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), casapool%Nsoilmin )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(13), casaflux%fHarvest )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(14), casaflux%fCrop )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(1), phen%phase )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(2), phen%doyphase(:,3) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), casapool%cplant  )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID2(2), casapool%nplant  )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), casapool%clitter  )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(2), casapool%nlitter )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  
  STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), casapool%csoil )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), casapool%nsoil )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  IF (icycle ==3) then
     STATUS = NF90_PUT_VAR(FILE_ID, VID1(5), casapool%psoillab )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     
     STATUS = NF90_PUT_VAR(FILE_ID, VID1(6), casapool%psoilsorb )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     
     STATUS = NF90_PUT_VAR(FILE_ID, VID1(7), casapool%psoilocc )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     

     STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), casapool%psoil )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     
     STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), casapool%pplant  )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     
     STATUS = NF90_PUT_VAR(FILE_ID, VID3(3), casapool%plitter )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
      
  ENDIF
  ! Close NetCDF file:
  STATUS = NF90_close(FILE_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

END SUBROUTINE WRITE_CASA_RESTART_NC

#ifndef UM_BUILD
SUBROUTINE READ_CASA_RESTART_NC (  casamet, casapool, casaflux,phen )

  USE CASAVARIABLE
  USE phenvariable
  USE CABLE_COMMON_MODULE
  USE CABLE_DEF_TYPES_MOD, ONLY: MET_TYPE, r_2, mp
  USE netcdf

  IMPLICIT NONE

  !INTEGER, INTENT(in)    :: YEAR
  TYPE (casa_met) , INTENT(inout) :: casamet
  TYPE (casa_pool), INTENT(inout) :: casapool
  TYPE (casa_flux), INTENT(inout) :: casaflux
  TYPE (phen_variable),       INTENT(INOUT) :: phen

  INTEGER*4 :: mp4
  INTEGER*4, parameter   :: pmp4 =0
  INTEGER, parameter   :: fmp4 = kind(pmp4)
  INTEGER*4   :: STATUS, i
  INTEGER*4   :: FILE_ID, dID, land_dim, mp_dim, ml_dim, ms_dim
  CHARACTER :: FRST_IN*99, CYEAR*4, CDATE*12, RSTDATE*12, FNAME*99

  ! ! 1 dim arrays (npt )
  ! CHARACTER(len=20),DIMENSION(7), PARAMETER :: A1 = (/ 'latitude', 'longitude', 'glai', &
  !      'clabile', 'psoillab','psoilsorb','psoilocc' /)
  ! ! 2 dim arrays (npt,mplant)
  ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A2 = (/ 'cplant' , 'nplant' , 'pplantc' /)
  ! ! 2 dim arrays (npt,mlitter)
  ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A3 = (/ 'clitter', 'nlitter', 'plitter' /)
  ! ! 2 dim arrays (npt,msoil)
  ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A4 = (/ 'csoil', 'nsoil', 'psoil' /)
  REAL(r_2), DIMENSION(mp)          :: LAT, LON, TMP
  REAL(r_2)                         :: TMP2(mp,mplant),TMP3(mp,mlitter),TMP4(mp,msoil)

  ! 1 dim arrays (npt )
  CHARACTER(len=20),DIMENSION(14) :: A1
  CHARACTER(len=20),DIMENSION(2) :: AI1
  ! 2 dim arrays (npt,mplant)
  CHARACTER(len=20),DIMENSION(3) :: A2
  ! 2 dim arrays (npt,mlitter)
  CHARACTER(len=20),DIMENSION(3) :: A3
  ! 2 dim arrays (npt,msoil)
  CHARACTER(len=20),DIMENSION(3) :: A4
  INTEGER :: VID1(SIZE(A1)), VID2(SIZE(A2)), VID3(SIZE(A3)), VID4(SIZE(A4))
  LOGICAL            ::  EXISTFILE, EXISTFILE1
  mp4=int(mp,fmp4)
  A1(1) = 'latitude'
  A1(2) = 'longitude'
  A1(3) = 'glai'
  A1(4) = 'clabile'
  A1(5) = 'psoillab'
  A1(6) = 'psoilsorb'
  A1(7) = 'psoilocc'
  A1(8) = 'frac_sapwood'
  A1(9) = 'sapwood_area'
  A1(10) = 'phen'
  A1(11) = 'aphen'
  A1(12) = 'nsoilmin'
  A1(13) = 'fHarvest'
  A1(14) = 'fCrop'

  AI1(1) = 'phase'
  AI1(2) = 'doyphase3'

  A2(1) = 'cplant'
  A2(2) = 'nplant'
  A2(3) = 'pplant'
  A3(1) = 'clitter'
  A3(2) = 'nlitter'
  A3(3) = 'plitter'
  A4(1) = 'csoil'
  A4(2) = 'nsoil'
  A4(3) = 'psoil'

!fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
!       '_casa_rst.nc'
  fname =  TRIM(casafile%cnpipool)
  INQUIRE( FILE=TRIM(fname), EXIST=EXISTFILE )
  IF (EXISTFILE) THEN
     STATUS = NF90_OPEN( TRIM(fname), NF90_NOWRITE, FILE_ID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     PRINT *, 'initial pool from restart file: ', fname
  ELSE
     write(*,*) 'CASA restart file:', TRIM(fname), ' does not exist'
     fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
          '_casa_rst.nc'   
     INQUIRE( FILE=TRIM(fname), EXIST=EXISTFILE1 )
     IF (EXISTFILE1) THEN
        STATUS = NF90_OPEN( TRIM(fname), NF90_NOWRITE, FILE_ID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        PRINT *, 'initial pool from restart file: ', fname
     ELSE
        write(*,*) 'CASA restart file:', TRIM(fname), ' does not exist either'
        write(*,*) 'Set cable_user%CASA_fromZero to true to initialise without restart file.'
        write(*,*) 'Otherwise set casafile%cnpipool to netcdf restart file name in cable.nml'
        stop
     ENDIF
  ENDIF

  ! TIME
  STATUS = NF90_GET_ATT( FILE_ID, NF90_GLOBAL, "Valid restart date", RSTDATE )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!$
  WRITE(CYEAR, FMT="(I4)") CurYear
  CDATE = '01/01/'//CYEAR
  ! compare current year with restart year (only for non-site type met data)
  IF ( CDATE .NE. RSTDATE .and. &
      TRIM(cable_user%MetType).NE.'' .and. TRIM(cable_user%MetType).NE.'site' ) THEN
     WRITE(*,*)"Restart Date in rst file doesn't match start date of Run!"
     WRITE(*,*)"File: "//RSTDATE//' Run: '//CDATE
    ! STOP
  ENDIF

  ! DIMS
  STATUS = NF90_INQ_DIMID( FILE_ID, 'land', dID )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=land_dim )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  STATUS = NF90_INQ_DIMID( FILE_ID, 'mplant', dID )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=mp_dim )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  STATUS = NF90_INQ_DIMID( FILE_ID, 'mlitter', dID )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=ml_dim )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  STATUS = NF90_INQ_DIMID( FILE_ID, 'msoil', dID )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=ms_dim )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  IF ( land_dim .NE. SIZE(casamet%lon) .OR. mp_dim .NE. mplant .OR. &
       ml_dim   .NE. mlitter             .OR. ms_dim .NE. msoil ) THEN
     WRITE(*,*)"Dimension misfit!"
     WRITE(*,*)"Restart file      Run"
     WRITE(*,*)"# points  ",land_dim,"     ",SIZE(casamet%lon)
     WRITE(*,*)"# mplant  ",mp_dim,"     ",mplant
     WRITE(*,*)"# mlitter ",ml_dim,"     ",mlitter
     WRITE(*,*)"# msoil   ",ms_dim,"     ",msoil
     STOP
  ENDIF

  ! LAT & LON
  STATUS = NF90_INQ_VARID( FILE_ID, A1(1), dID )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_GET_VAR( FILE_ID, dID, LAT )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  STATUS = NF90_INQ_VARID( FILE_ID, A1(2), dID )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_GET_VAR( FILE_ID, dID, LON )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  ! CHECK FOR VALID LONS

  ! READ 1-dimensional fields
  DO i = 3, SIZE(A1)
     STATUS = NF90_INQ_VARID( FILE_ID, A1(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A1(i)))
     CASE ('glai'      ) ; casamet%glai       = TMP
     CASE ('clabile'   ) ; casapool%clabile   = TMP
     CASE ('frac_sapwood' ) ; casaflux%frac_sapwood  = TMP
     CASE ( 'sapwood_area' ) ; casaflux%sapwood_area  = TMP
     CASE ( 'phen' ) ; phen%phen  = TMP
     CASE ( 'aphen' ) ; phen%aphen  = TMP
     CASE ( 'nsoilmin' ) ; casapool%Nsoilmin  = TMP
     CASE ( 'fHarvest' ) ; casaflux%fHarvest  = TMP
     CASE ( 'fCrop' ) ; casaflux%fCrop  = TMP
     END SELECT
  END DO
IF (icycle==3) then
  DO i = 3, SIZE(A1)
     STATUS = NF90_INQ_VARID( FILE_ID, A1(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A1(i)))
     CASE ('psoillab'  ) ; casapool%psoillab  = TMP
     CASE ('psoilsorb' ) ; casapool%psoilsorb = TMP
     CASE ('psoilocc'  ) ; casapool%psoilocc  = TMP
     END SELECT
  END DO
ENDIF

  DO i = 1, SIZE(AI1)
     STATUS = NF90_INQ_VARID( FILE_ID, AI1(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(AI1(i)))
     CASE ( 'phase' ) ; phen%phase  = TMP
     CASE ( 'doyphase3' ) ; phen%doyphase(:,3)  = TMP
     END SELECT
  END DO

  ! READ 2-dimensional fields (mplant)
  DO i = 1, SIZE(A2)
     STATUS = NF90_INQ_VARID( FILE_ID, A2(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP2 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A2(i)))
     CASE ('cplant' ) ; casapool%cplant = TMP2
     CASE ('nplant' ) ; casapool%nplant = TMP2
     END SELECT
  END DO


IF (icycle==3) then
   DO i = 1, SIZE(A2)
     STATUS = NF90_INQ_VARID( FILE_ID, A2(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP2 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A2(i)))
     CASE ('pplant' ) ; casapool%pplant = TMP2
     END SELECT
  END DO
ENDIF

  ! READ 2-dimensional fields (mlitter)
  DO i = 1, SIZE(A3)
     STATUS = NF90_INQ_VARID( FILE_ID, A3(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP3 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A3(i)))
     CASE ('clitter' ) ; casapool%clitter = TMP3
     CASE ('nlitter' ) ; casapool%nlitter = TMP3
     END SELECT
  END DO

IF (icycle==3) then

  DO i = 1, SIZE(A3)
     STATUS = NF90_INQ_VARID( FILE_ID, A3(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP3 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A3(i)))
     CASE ('plitter' ) ; casapool%plitter = TMP3
     END SELECT
  END DO


ENDIF

  ! READ 2-dimensional fields (msoil)
  DO i = 1, SIZE(A4)
     STATUS = NF90_INQ_VARID( FILE_ID, A4(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP4 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     SELECT CASE ( TRIM(A4(i)))
     CASE ('csoil' ) ; casapool%csoil = TMP4
     CASE ('nsoil' ) ; casapool%nsoil = TMP4
     END SELECT
  END DO
IF (icycle==3) then
 DO i = 1, SIZE(A4)
     STATUS = NF90_INQ_VARID( FILE_ID, A4(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP4 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A4(i)))
     CASE ('psoil' ) ; casapool%psoil = TMP4
     END SELECT
  END DO
ENDIF

  STATUS = NF90_CLOSE( FILE_ID )

END SUBROUTINE READ_CASA_RESTART_NC
#endif
SUBROUTINE WRITE_CASA_OUTPUT_NC ( veg, casamet, casapool, casabal, casaflux, &
     CASAONLY, ctime, FINAL )

  USE CASAVARIABLE
  USE CABLE_COMMON_MODULE

  
  USE cable_def_types_mod, ONLY: veg_parameter_type
  
  USE netcdf

  IMPLICIT NONE

  TYPE (casa_met) ,   INTENT(in) :: casamet
  TYPE (casa_pool),   INTENT(in) :: casapool
  TYPE (casa_balance),INTENT(in) :: casabal
  TYPE (casa_flux),   INTENT(in) :: casaflux
  TYPE (veg_parameter_type),  INTENT(IN) :: veg  ! vegetation parameters

  INTEGER   :: STATUS, ctime
  INTEGER   :: land_ID, plnt_ID, litt_ID, soil_ID, t_ID, i
  LOGICAL   :: CASAONLY, FINAL
  CHARACTER :: CYEAR*4, FNAME*99,dum*50
  LOGICAL, SAVE :: CALL1 = .TRUE.

  ! ! 1 dim arrays (mp )
  ! CHARACTER(len=20),DIMENSION(2), PARAMETER :: A0 = (/ 'latitude', 'longitude' /)
  ! ! 2 dim arrays (mp,t)
  ! CHARACTER(len=20),DIMENSION(44),PARAMETER :: A1 = (/ 'glai', 'clabile',      &
  !      'psoillab','psoilsorb','psoilocc', 'sumcbal','sumnbal','sumpbal','Cgpp',&
  !      'Cnpp','stemnpp','Crp','Crgplant','Nminfix','Plabuptake','Clabloss',    &
  !      'fraclabile','Cnep','Crsoil','Nmindep','Nminloss','Nminleach',          &
  !      'Nupland','Nlittermin','Nsmin','Nsnet','fNMinloss','fNMinleach','Pdep', &
  !      'pwea','Pleach','Ploss','Pupland','Plittermin','Psmin','Psimm','Psnet', &
  !      'fPleach','kPlab','kPsorb','kpocc','kmlabP','Psorbmax','FluxCtoco2'/)
  ! ! 3 dim arrays (mp,mplant,t)
  ! CHARACTER(len=20),DIMENSION(8), PARAMETER :: A2 = (/ 'cplant' , 'nplant' ,   &
  !      'pplantc','fracCalloc','fracNalloc','fracPalloc','kplant','Crmplant'/)
  ! ! 3 dim arrays (mp,mlitter,t)
  ! CHARACTER(len=20),DIMENSION(8), PARAMETER :: A3 = (/ 'clitter', 'nlitter',   &
  !      'plitter','klitter','fromL2CO2','FluxCtolitter','FluxNtolitter',        &
  !      'FluxPtolitter' /)
  ! ! 3 dim arrays (mp,msoil,t)
  ! CHARACTER(len=20),DIMENSION(8), PARAMETER :: A4 = (/ 'csoil','nsoil','psoil',&
  !      'ksoil','fromStoCO2','FluxCtosoil','FluxNtosoil','FluxPxtosoil'/)

  ! 1 dim arrays (mp )
  CHARACTER(len=20),DIMENSION(2) :: A0
  ! 2 dim arrays (mp,t)
  CHARACTER(len=20),DIMENSION(51):: A1
  ! 3 dim arrays (mp,mplant,t)
  CHARACTER(len=20),DIMENSION(8) :: A2
  ! 3 dim arrays (mp,mlitter,t)
  CHARACTER(len=20),DIMENSION(8) :: A3
  ! 3 dim arrays (mp,msoil,t)
  CHARACTER(len=20),DIMENSION(8) :: A4

  ! 4 dim arrays (mp,mlitter,mplant,t)
  CHARACTER(len=20),DIMENSION(1), PARAMETER :: A5 = (/ 'fromPtoL'/)
  ! 4 dim arrays (mp,msoil,mlitter,t)
  CHARACTER(len=20),DIMENSION(1), PARAMETER :: A6 = (/ 'fromLtoS'/)
  ! 4 dim arrays (mp,msoil,msoil,t)
  CHARACTER(len=20),DIMENSION(1), PARAMETER :: A7 = (/ 'fromStoS'/)

  INTEGER, SAVE :: VIDtime, VID0(SIZE(A0)),VID1(SIZE(A1)),VID2(SIZE(A2)),VID3(SIZE(A3))
  INTEGER, SAVE :: VID4(SIZE(A4)),VID5(SIZE(A5)),VID6(SIZE(A6)),VID7(SIZE(A7))
  INTEGER, SAVE :: FILE_ID, CNT = 0
  LOGICAL   :: EXRST
  CHARACTER(len=50) :: RecordDimName


  A0(1) = 'latitude'
  A0(2) = 'longitude'

  A1(1) = 'glai'
  A1(2) = 'clabile'
  A1(3) = 'psoillab'
  A1(4) = 'psoilsorb'
  A1(5) = 'psoilocc'
  A1(6) = 'sumcbal'
  A1(7) = 'sumnbal'
  A1(8) = 'sumpbal'
  A1(9) = 'Cgpp'
  A1(10) = 'Cnpp'
  A1(11) = 'stemnpp'
  A1(12) = 'Crp'
  A1(13) = 'Crgplant'
  A1(14) = 'Nminfix'
  A1(15) = 'Plabuptake'
  A1(16) = 'Clabloss'
  A1(17) = 'fraclabile'
  A1(18) = 'Cnep'
  A1(19) = 'Crsoil'
  A1(20) = 'Nmindep'
  A1(21) = 'Nminloss'
  A1(22) = 'Nminleach'
  A1(23) = 'Nupland'
  A1(24) = 'Nlittermin'
  A1(25) = 'Nsmin'
  A1(26) = 'Nsimm'
  A1(27) = 'Nsnet'
  A1(28) = 'fNMinloss'
  A1(29) = 'Pdep'
  A1(30) = 'pwea'
  A1(31) = 'Pleach'
  A1(32) = 'Ploss'
  A1(33) = 'Pupland'
  A1(34) = 'Plittermin'
  A1(35) = 'Psmin'
  A1(36) = 'Psimm'
  A1(37) = 'Psnet'
  A1(38) = 'fPleach'
  A1(39) = 'kPlab'
  A1(40) = 'kPsorb'
  A1(41) = 'kpocc'
  A1(42) = 'kmlabP'
  A1(43) = 'Psorbmax'
  A1(44) = 'FluxCtoco2'
  A1(45) = 'FCgppyear'
  A1(46) = 'FCrpyear'
  A1(47) = 'FCnppyear'
  A1(48) = 'FCrsyear'
  A1(49) = 'FCNeeyear'
  A1(50) = 'vcmax'
  A1(51) = 'Nsoilmin'

  A2(1) = 'cplant'
  A2(2) = 'nplant'
  A2(3) = 'pplant'
  A2(4) = 'fracCalloc'
  A2(5) = 'fracNalloc'
  A2(6) = 'fracPalloc'
  A2(7) = 'kplant'
  A2(8) = 'Crmplant'

  A3(1) = 'clitter'
  A3(2) = 'nlitter'
  A3(3) = 'plitter'
  A3(4) = 'klitter'
  A3(5) = 'fromL2CO2'
  A3(6) = 'FluxCtolitter'
  A3(7) = 'FluxNtolitter'
  A3(8) = 'FluxPtolitter'

  A4(1) = 'csoil'
  A4(2) = 'nsoil'
  A4(3) = 'psoil'
  A4(4) = 'ksoil'
  A4(5) = 'fromStoCO2'
  A4(6) = 'FluxCtosoil'
  A4(7) = 'FluxNtosoil'
  A4(8) = 'FluxPxtosoil'


  CNT = CNT + 1

  IF ( CALL1 ) THEN
     ! Get File-Name
     IF (TRIM(casafile%out).NE.'' ) THEN
        fname = TRIM(casafile%out)
     ELSE
        IF (TRIM(cable_user%MetType).NE.'' ) THEN

           WRITE( dum, FMT="(I4,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
           IF (CABLE_USER%YEARSTART.lt.1000.and.CABLE_USER%YEAREND.lt.1000) THEN
              WRITE( dum, FMT="(I3,'_',I3)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
           ELSEIF (CABLE_USER%YEARSTART.lt.1000) THEN
              WRITE( dum, FMT="(I3,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
           ENDIF
           fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//&
                TRIM(dum)//'_casa_out.nc'
        ELSE
           ! site data
           fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_casa_out.nc'
        ENDIF
     ENDIF

     INQUIRE( FILE=TRIM( fname ), EXIST=EXRST )
     EXRST = .FALSE.
     IF ( EXRST ) THEN
        STATUS = NF90_open(fname, mode=nf90_write, ncid=FILE_ID)
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        CALL1 = .FALSE.

        STATUS = nf90_inq_dimid(FILE_ID, 'time', t_id)
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

!CRM        status = nf90_inquire_dimension(FILE_ID, t_id,name = RecordDimName, len = CNT)
!CRM        if (status /= nf90_noerr) call handle_err(status)
!CRM        CNT = CNT+1

        STATUS = nf90_inq_varid(FILE_ID, 'time', VIDTime)
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

        DO i = 1, SIZE(A0)
           STATUS = nf90_inq_varid(FILE_ID,TRIM(A0(i)),VID0(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        DO i = 1, SIZE(A1)
           STATUS = nf90_inq_varid(FILE_ID,TRIM(A1(i)), VID1(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        DO i = 1, SIZE(A2)
           STATUS = nf90_inq_varid(FILE_ID,TRIM(A2(i)) , VID2(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        DO i = 1, SIZE(A3)
           STATUS = nf90_inq_varid(FILE_ID,TRIM(A3(i)) ,VID3(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        DO i = 1, SIZE(A4)
           STATUS = nf90_inq_varid(FILE_ID,TRIM(A4(i)) ,VID4(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        DO i = 1, SIZE(A5)
           STATUS = nf90_inq_varid(FILE_ID,TRIM(A5(i)), VID5(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        DO i = 1, SIZE(A6)
           STATUS = nf90_inq_varid(FILE_ID,TRIM(A6(i)), VID6(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        DO i = 1, SIZE(A7)
           STATUS = nf90_inq_varid(FILE_ID,TRIM(A7(i)),VID7(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

     ELSE
        ! Create NetCDF file:
        STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

        ! Put the file in define mode:
        STATUS = NF90_redef(FILE_ID)

        STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Icycle"   , icycle  )
        STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "StartYear", CABLE_USER%YEARSTART )
        STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "EndYear"  , CABLE_USER%YEAREND   )
        STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "RunIden"  , CABLE_USER%RunIden   )
        IF ( CASAONLY ) THEN
           dum = 'CASA-ONLY run'
        ELSE
           dum = 'CABLE-CASA coupled run'
        ENDIF
        STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Run-Type", TRIM(dum) )

        ! Define dimensions:
        ! Land (number of points)
        STATUS = NF90_def_dim(FILE_ID, 'land'   , mp     , land_ID)
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_def_dim(FILE_ID, 'mplant' , mplant , plnt_ID)
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_def_dim(FILE_ID, 'mlitter', mlitter, litt_ID)
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_def_dim(FILE_ID, 'msoil'  , msoil  , soil_ID)
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_def_dim(FILE_ID, 'time'   , NF90_UNLIMITED, t_ID)
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

        ! Define variables
        STATUS = NF90_def_var(FILE_ID,'time' ,NF90_INT,(/t_ID/),VIDtime )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

        DO i = 1, SIZE(A0)
           STATUS = NF90_def_var(FILE_ID,TRIM(A0(i)) ,NF90_FLOAT,(/land_ID/),VID0(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        DO i = 1, SIZE(A1)
           STATUS = NF90_def_var(FILE_ID,TRIM(A1(i)) ,NF90_FLOAT,(/land_ID,t_ID/),VID1(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        DO i = 1, SIZE(A2)
           STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/land_ID,plnt_ID,t_ID/),VID2(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        DO i = 1, SIZE(A3)
           STATUS = NF90_def_var(FILE_ID,TRIM(A3(i)) ,NF90_FLOAT,(/land_ID,litt_ID,t_ID/),VID3(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        DO i = 1, SIZE(A4)
           STATUS = NF90_def_var(FILE_ID,TRIM(A4(i)) ,NF90_FLOAT,(/land_ID,soil_ID,t_ID/),VID4(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        DO i = 1, SIZE(A5)
           STATUS = NF90_def_var(FILE_ID,TRIM(A5(i)) ,NF90_FLOAT, &
                (/land_ID,litt_ID,plnt_ID,t_ID/),VID5(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        DO i = 1, SIZE(A6)
           STATUS = NF90_def_var(FILE_ID,TRIM(A6(i)) ,NF90_FLOAT, &
                (/land_ID,soil_ID,litt_ID,t_ID/),VID6(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        DO i = 1, SIZE(A7)
           STATUS = NF90_def_var(FILE_ID,TRIM(A7(i)) ,NF90_FLOAT, &
                (/land_ID,soil_ID,soil_ID,t_ID/),VID7(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        ! End define mode:
        STATUS = NF90_enddef(FILE_ID)
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


        ! PUT LAT / LON ( mp )
        STATUS = NF90_PUT_VAR(FILE_ID, VID0(1), casamet%lat )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

        STATUS = NF90_PUT_VAR(FILE_ID, VID0(2), casamet%lon )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

        CALL1 = .FALSE.
     ENDIF !( EXRST )
  ENDIF

  ! TIME  ( t )
  STATUS = NF90_PUT_VAR(FILE_ID, VIDtime, ctime, start=(/ CNT /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

SELECT CASE(icycle)
CASE(1)
  ! PUT 2D VARS ( mp, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), casamet%glai,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), casapool%clabile,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 6), casabal%sumcbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 7), casabal%sumnbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 9), casaflux%Cgpp,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), REAL(casaflux%Cnpp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), casaflux%stemnpp,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), casaflux%Crp,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(13), casaflux%Crgplant,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(16), casaflux%Clabloss,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(17), casaflux%fracClabile,start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(18), casaflux%Cnep,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(19), casaflux%Crsoil,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  
  
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(45), casabal%FCgppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(46), casabal%FCrpyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(47), casabal%FCnppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(48), casabal%FCrsyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(49), casabal%FCneeyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(50), veg%vcmax,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

! PUT 3D VARS ( mp, mplant, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), casapool%cplant,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(4), casaflux%fracCalloc,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(7), casaflux%kplant,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(8), casaflux%Crmplant,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

! PUT 3D VARS ( mp, mlitter, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), casapool%clitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(4), casaflux%klitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(5), casaflux%fromLtoCO2,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(6), casaflux%FluxCtolitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


  ! PUT 3D VARS ( mp, msoil, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), casapool%csoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(4), casaflux%ksoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(5), casaflux%fromStoCO2,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(6), casaflux%FluxCtosoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


! PUT 4D VARS ( mp, mlitter,mplant, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), casaflux%fromPtoL,   &
       start=(/ 1,1,1,CNT /), count=(/ mp,mlitter,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT 4D VARS ( mp, msoil, mlitter, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID6(1), casaflux%fromLtoS,   &
       start=(/ 1,1,1,CNT /), count=(/ mp,msoil,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT 4D VARS ( mp, msoil, msoil, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID7(1), casaflux%fromStoS,   &
       start=(/ 1,1,1,CNT /), count=(/ mp,msoil,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

CASE(2)


  ! PUT 2D VARS ( mp, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), casamet%glai,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), casapool%clabile,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 3), casapool%psoillab,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 4), casapool%psoilsorb,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 5), casapool%psoilocc,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 6), casabal%sumcbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 7), casabal%sumnbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 8), casabal%sumpbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 9), casaflux%Cgpp,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), casaflux%Cnpp,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), casaflux%stemnpp,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), casaflux%Crp,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(13), casaflux%Crgplant,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(14), casaflux%Nminfix,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
 ! STATUS = NF90_PUT_VAR(FILE_ID, VID1(15), casaflux%Nminuptake, start=(/ 1, CNT /), count=(/ mp, 1 /) )
 ! IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(15), casaflux%Plabuptake, start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(16), casaflux%Clabloss,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(17), casaflux%fracClabile,start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(18), casaflux%Cnep,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(19), casaflux%Crsoil,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(20), casaflux%Nmindep,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(21), casaflux%Nminloss,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(22), casaflux%Nminleach,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(23), casaflux%Nupland,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(24), casaflux%Nlittermin, start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(25), casaflux%Nsmin,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(26), casaflux%Nsimm,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(27), casaflux%Nsnet,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(28), casaflux%fNminloss,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(29), casaflux%Pdep,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(30), casaflux%Pwea,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(31), casaflux%Pleach,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(32), casaflux%Ploss,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(33), casaflux%Pupland,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(34), casaflux%Plittermin, start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(35), casaflux%Psmin,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(36), casaflux%Psimm,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(37), casaflux%Psnet,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(38), casaflux%fPleach,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(39), casaflux%kplab,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(40), casaflux%kpsorb,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(41), casaflux%kpocc,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(42), casaflux%kmlabP,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(43), casaflux%Psorbmax,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(44), casaflux%FluxCtoco2,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(45), casabal%FCgppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(46), casabal%FCrpyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(47), casabal%FCnppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(48), casabal%FCrsyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(49), casabal%FCneeyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(50), veg%vcmax,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(51), casapool%Nsoilmin,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT 3D VARS ( mp, mplant, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), casapool%cplant,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(2), casapool%nplant,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), casapool%pplant,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(4), casaflux%fracCalloc,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(5), casaflux%fracNalloc,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(6), casaflux%fracPalloc,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(7), casaflux%kplant,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(8), casaflux%Crmplant,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT 3D VARS ( mp, mlitter, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), casapool%clitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(2), casapool%nlitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(3), casapool%plitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(4), casaflux%klitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(5), casaflux%fromLtoCO2,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(6), casaflux%FluxCtolitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(7), casaflux%FluxNtolitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(8), casaflux%FluxPtolitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT 3D VARS ( mp, msoil, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), casapool%csoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), casapool%nsoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), casapool%psoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(4), casaflux%ksoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(5), casaflux%fromStoCO2,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(6), casaflux%FluxCtosoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(7), casaflux%FluxNtosoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(8), casaflux%FluxPtosoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT 4D VARS ( mp, mlitter,mplant, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), casaflux%fromPtoL,   &
       start=(/ 1,1,1,CNT /), count=(/ mp,mlitter,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT 4D VARS ( mp, msoil, mlitter, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID6(1), casaflux%fromLtoS,   &
       start=(/ 1,1,1,CNT /), count=(/ mp,msoil,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT 4D VARS ( mp, msoil, msoil, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID7(1), casaflux%fromStoS,   &
       start=(/ 1,1,1,CNT /), count=(/ mp,msoil,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

CASE(3)
  ! PUT 2D VARS ( mp, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), casamet%glai,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), casapool%clabile,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 3), casapool%psoillab,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 4), casapool%psoilsorb,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 5), casapool%psoilocc,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 6), casabal%sumcbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 7), casabal%sumnbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 8), casabal%sumpbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1( 9), casaflux%Cgpp,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), casaflux%Cnpp,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), casaflux%stemnpp,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), casaflux%Crp,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(13), casaflux%Crgplant,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(14), casaflux%Nminfix,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
!  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!  STATUS = NF90_PUT_VAR(FILE_ID, VID1(15), casaflux%Nminuptake, start=(/ 1, CNT /), count=(/ mp, 1 !/) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(15), casaflux%Plabuptake, start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(16), casaflux%Clabloss,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(17), casaflux%fracClabile,start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(18), casaflux%Cnep,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(19), casaflux%Crsoil,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(20), casaflux%Nmindep,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(21), casaflux%Nminloss,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(22), casaflux%Nminleach,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(23), casaflux%Nupland,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(24), casaflux%Nlittermin, start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(25), casaflux%Nsmin,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(26), casaflux%Nsimm,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(27), casaflux%Nsnet,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(28), casaflux%fNminloss,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(29), casaflux%Pdep,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(30), casaflux%Pwea,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(31), casaflux%Pleach,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(32), casaflux%Ploss,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(33), casaflux%Pupland,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(34), casaflux%Plittermin, start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(35), casaflux%Psmin,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(36), casaflux%Psimm,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(37), casaflux%Psnet,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(38), casaflux%fPleach,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(39), casaflux%kplab,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(40), casaflux%kpsorb,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(41), casaflux%kpocc,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(42), casaflux%kmlabP,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(43), casaflux%Psorbmax,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(44), casaflux%FluxCtoCo2,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(45), casabal%FCgppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(46), casabal%FCrpyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(47), casabal%FCnppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(48), casabal%FCrsyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(49), casabal%FCneeyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(50), veg%vcmax,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(51), casapool%Nsoilmin,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


  ! PUT 3D VARS ( mp, mplant, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), casapool%cplant,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(2), casapool%nplant,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), casapool%pplant,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(4), casaflux%fracCalloc,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(5), casaflux%fracNalloc,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(6), casaflux%fracPalloc,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(7), casaflux%kplant,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  STATUS = NF90_PUT_VAR(FILE_ID, VID2(8), casaflux%Crmplant,   &
       start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT 3D VARS ( mp, mlitter, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), casapool%clitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(2), casapool%nlitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(3), casapool%plitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(4), casaflux%klitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(5), casaflux%fromLtoCO2,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(6), casaflux%FluxCtolitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(7), casaflux%FluxNtolitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(8), casaflux%FluxPtolitter,   &
       start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT 3D VARS ( mp, msoil, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), casapool%csoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), casapool%nsoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), casapool%psoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(4), casaflux%ksoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(5), casaflux%fromStoCO2,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(6), casaflux%FluxCtosoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(7), casaflux%FluxNtosoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(8), casaflux%FluxPtosoil,   &
       start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT 4D VARS ( mp, mlitter,mplant, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), casaflux%fromPtoL,   &
       start=(/ 1,1,1,CNT /), count=(/ mp,mlitter,mplant,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT 4D VARS ( mp, msoil, mlitter, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID6(1), casaflux%fromLtoS,   &
       start=(/ 1,1,1,CNT /), count=(/ mp,msoil,mlitter,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT 4D VARS ( mp, msoil, msoil, t )
  STATUS = NF90_PUT_VAR(FILE_ID, VID7(1), casaflux%fromStoS,   &
       start=(/ 1,1,1,CNT /), count=(/ mp,msoil,msoil,1 /) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

END SELECT
  IF ( FINAL ) THEN
     ! Close NetCDF file:
     STATUS = NF90_close(FILE_ID)
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     WRITE(*,*) " Casa Output written to ", TRIM(fname)
  ENDIF

END SUBROUTINE WRITE_CASA_OUTPUT_NC
#endif
