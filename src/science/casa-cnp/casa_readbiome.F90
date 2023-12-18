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

MODULE casa_readbiome_module

CONTAINS

  SUBROUTINE casa_readbiome(veg,soil,casabiome,casapool,casaflux,casamet,phen)
    ! mst actually not used in this routine (BP sep2010)
    !SUBROUTINE casa_readbiome(mvt,mst,veg,soil, &
    !                          casabiome,casapool,casaflux,casamet,phen)
    USE cable_def_types_mod
    USE casadimension
    USE casaparm
USE casavariable, ONLY : casafile
USE casavariable, ONLY :  casa_biome
USE casavariable, ONLY :  casa_pool
USE casavariable, ONLY :  casa_flux
USE casavariable, ONLY :  casa_met
   
    USE phenvariable
    !! vh_js !!
    USE cable_common_module, ONLY: cable_user
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
    INTEGER :: nv0,nv1,nv2,nv3,nv4,nv5,nv6,nv7,nv8,nv9,nv10,nv11,nv12
    REAL(r_2), DIMENSION(mvtype)       :: xxnpmax,xq10soil,xxkoptlitter,xxkoptsoil,xprodptase, &
         xcostnpup,xmaxfinelitter,xmaxcwd,xnintercept,xnslope
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
       IF (casamet%iveg2(npt)==forest.OR.casamet%iveg2(npt)==shrub) THEN
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
       casaflux%kmlabp(npt)      = xkmlabp(iso)
       casaflux%psorbmax(npt)    = xpsorbmax(iso)
       casaflux%fpleach(npt)     = xfPleach(iso) 
       !   we used the spatially explicit estimate N fixation by Wang and Houlton (GRL)
       !    casaflux%Nminfix(npt)     = xnfixrate(iv1)/365.0

       ! use the PFT-specific C:N:P stoichiometry
       casapool%ratioNCplant(npt,:)  = 1.0/ratioCNplant(iv1,:)
       casapool%ratioNPplant(npt,:)  = casabiome%ratioNPplantmin(iv1,:)
       casapool%ratioPCplant(npt,:)  = 1.0/(ratioCNplant(iv1,:) *casabiome%ratioNPplantmin(iv1,:) )

       casapool%ratioNClitter(npt,metb) = casapool%ratioNCplant(npt,leaf)  * casabiome%ftransNPtoL(iv1,leaf)
       casapool%ratioNClitter(npt,str)  = casapool%ratioNCplant(npt,froot) * casabiome%ftransNPtoL(iv1,froot)
       casapool%ratioNClitter(npt,cwd)  = casapool%ratioNCplant(npt,wood)  * casabiome%ftransNPtoL(iv1,wood)

       casapool%ratioPClitter(npt,metb) = casapool%ratioPCplant(npt,leaf)  * casabiome%ftransPPtoL(iv1,leaf)
       casapool%ratioPClitter(npt,str)  = casapool%ratioPCplant(npt,froot) * casabiome%ftransPPtoL(iv1,froot)
       casapool%ratioPClitter(npt,cwd)  = casapool%ratioPCplant(npt,wood)  * casabiome%ftransPPtoL(iv1,wood)

       casapool%ratioNPlitter(npt,metb) = casapool%ratioNClitter(npt,metb)/(casapool%ratioPClitter(npt,metb) +1.0e-10)
       casapool%ratioNPlitter(npt,str)  = casapool%ratioNClitter(npt,str)/(casapool%ratioPClitter(npt,str) +1.0e-10)
       casapool%ratioNPlitter(npt,cwd)  = casapool%ratioNClitter(npt,cwd)/(casapool%ratioPClitter(npt,cwd) +1.0e-10)

       casapool%ratioNCsoil(npt,:)   = 1.0/ratioCNsoil(iv1,:)
       casapool%ratioNPsoil(npt,:)   = ratioNPsoil(iso,:)
       casapool%ratioPCsoil(npt,:)   = 1.0/(ratioCNsoil(iv1,:)*ratioNPsoil(iso,:))

       casapool%ratioNCsoilmin(npt,:)   = 1.0/ratioCNsoilmax(iv1,:)
       casapool%ratioNCsoilmax(npt,:)   = 1.0/ratioCNsoilmin(iv1,:)
       casapool%ratioNCsoilnew(npt,:)   = casapool%ratioNCsoilmax(npt,:)
    ENDDO

    IF(icycle<2) THEN
       casapool%Nplant(:,:)  = casapool%Cplant(:,:)  * casapool%ratioNCplant(:,:)
       casapool%Nlitter(:,:) = casapool%Clitter(:,:) * casapool%ratioNClitter(:,:)
       casapool%Nsoil(:,:)   = casapool%Csoil(:,:)   * casapool%ratioNCsoil(:,:) 
    ENDIF
    IF(icycle<3) THEN
       casapool%Pplant(:,:)  = casapool%Cplant(:,:)  * casapool%ratioPCplant(:,:)
       casapool%Plitter(:,:) = casapool%Clitter(:,:) * casapool%ratioPClitter(:,:)
       casapool%Psoil(:,:)   = casapool%Csoil(:,:)   * casapool%ratioPCsoil(:,:) 
       casapool%Psoilsorb(:) = casaflux%psorbmax(:) * casapool%psoillab(:) &
            /(casaflux%kmlabp(:)+casapool%psoillab(:))
    ENDIF


  END SUBROUTINE casa_readbiome

END MODULE casa_readbiome_module
