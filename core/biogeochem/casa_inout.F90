!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.accessimulator.org.au/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
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
!   casa_cnpflux  (not used?)
!   biogeochem

SUBROUTINE casa_readbiome(veg,soil,casabiome,casapool,casaflux,casamet,phen)
! mst actually not used in this routine (BP sep2010)
!SUBROUTINE casa_readbiome(mvt,mst,veg,soil, &
!                          casabiome,casapool,casaflux,casamet,phen)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
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
  REAL(r_2), DIMENSION(mvtype)       :: micage,slowage,passage,clabileage
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
  INTEGER :: nv0,nv1,nv2,nv3,nv4,nv5,nv6,nv7,nv8,nv9,nv10

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
                micage(nv),slowage(nv),passage(nv),clabileage(nv) 
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
!     PRINT *, 'nv8',nv8
  ENDDO

  READ(101,*)
  READ(101,*)
  DO iso=1,mso
    READ(101,*) nv9,xkmlabp(iso),xpsorbmax(iso),xfPleach(iso), &
                ratioNPsoil(iso,mic),ratioNPsoil(iso,slow),ratioNPsoil(iso,pass)
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
    casabiome%sla(nv)             = 0.025 * (leafage(nv)**(-0.5)) ! see eqn A1 of Arora and Boer, GCB, 2005
!    casabiome%sla(nv)             = 2.0E-4 * exp(6.15)/((12*leafage(nv))**0.46) ! see eqn 6 of Sitch, GCB, 2003
!    casabiome%fherbivore(nv)     = deltcasa*xfherbivore(nv)
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
  ENDDO

!  PRINT *, 'casabiome%ivt2 = ', casabiome%ivt2

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
      casapool%cplant(npt,wood)  = cwood(iv1) 
      casapool%clitter(npt,cwd)  = ccwd(iv1)
      casapool%nplant(npt,wood)  = nwood(iv1) 
      casapool%nlitter(npt,cwd)  = ncwd(iv1)
      casapool%pplant(npt,wood)  = xpwood(iv1)
      casapool%plitter(npt,cwd)  = xpcwd(iv1)
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
      casapool%rationcplant(npt,:)  = 1.0/ratiocnplant(iv1,:)
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
    casapool%plitter(npt,str) = cstr(iv1)*ratioPCstrfix
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

    casapool%rationcplant(npt,:)  = 1.0/ratioCNplant(iv1,:)
    casapool%ratiopcplant(npt,:)  = casabiome%ratioPcplantmax(iv1,:)
    casapool%rationclitter(npt,:) = casapool%nlitter(npt,:)/(casapool%clitter(npt,:)+1.0e-10)
    casapool%ratiopclitter(npt,:) = casapool%plitter(npt,:)/(casapool%clitter(npt,:)+1.0e-10)
    casapool%ratioNCsoil(npt,:)   = 1.0/ratioCNsoil(iv1,:)
    casapool%ratioPCsoil(npt,:)   = 1.0/(ratioCNsoil(iv1,:)*ratioNPsoil(iso,:))
    casapool%ratioNCsoilmin(npt,:)   = 1.0/ratioCNsoilmax(iv1,:)
    casapool%ratioNCsoilmax(npt,:)   = 1.0/ratioCNsoilmin(iv1,:)
    casapool%ratioNCsoilnew(npt,:)   = casapool%ratioNCsoilmax(npt,:)

  ENDDO

  IF (icycle==1) THEN
    casapool%nplant(:,:)  = casapool%cplant(:,:) * casapool%rationcplant(:,:)
  ELSE
    casapool%Nsoil(:,:)   = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
    casapool%Psoil(:,:)   = casapool%ratioPCsoil(:,:) * casapool%Csoil(:,:)
    casapool%psoilsorb(:) = casaflux%psorbmax(:) * casapool%psoillab(:) &
                            /(casaflux%kmlabp(:)+casapool%psoillab(:))
  ENDIF 
      
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

SUBROUTINE casa_init(casabiome,casamet,casapool,casabal,veg,phen)
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
  
! end addition (BP may2010)
  IMPLICIT NONE
!  INTEGER,        INTENT(IN)    :: mst
  TYPE (casa_biome),   INTENT(IN)    :: casabiome
  TYPE (casa_met),     INTENT(INOUT) :: casamet
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

  PRINT *, 'initial pool from ',TRIM(casafile%cnpipool)
  PRINT *, 'icycle,initcasa,mp ', icycle,initcasa,mp
  !phen%phase = 2
  IF (initcasa==1) THEN
    OPEN(99,file=casafile%cnpipool)

    DO npt =1, mp
      SELECT CASE(icycle)
      CASE(1)
        READ(99,92) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
                   casamet%glai(npt),slaz,phen%phase(npt), &
                   casapool%clabile(npt),casapool%cplant(npt,:),  &
                   casapool%clitter(npt,:),casapool%csoil(npt,:)
      CASE(2)
        READ(99,92) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
                   casamet%glai(npt),slaz,phen%phase(npt), &
                   casapool%clabile(npt),casapool%cplant(npt,:),   &
                   casapool%clitter(npt,:),casapool%csoil(npt,:),       &
                   casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
                   casapool%nsoil(npt,:),casapool%nsoilmin(npt)
      CASE(3)
        READ(99,92) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
                   casamet%glai(npt),slaz,phen%phase(npt), &
                   casapool%clabile(npt),casapool%cplant(npt,:),   &
                   casapool%clitter(npt,:),casapool%csoil(npt,:),       &
                   casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
                   casapool%nsoil(npt,:),casapool%nsoilmin(npt),        &
                   casapool%pplant(npt,:),casapool%plitter(npt,:),      &
                   casapool%psoil(npt,:),casapool%psoillab(npt),        &
                   casapool%psoilsorb(npt),casapool%psoilocc(npt)
      END SELECT 
      IF (ABS(patch(npt)%longitude - lonz) > 0.9 .OR. &
          ABS(patch(npt)%latitude  - latz) > 0.9) THEN
        PRINT *, 'patch(npt)%longitude, lonz:', patch(npt)%longitude, lonz
        PRINT *, 'patch(npt)%latitude,  latz:', patch(npt)%latitude,  latz
        PRINT *, 'npt = ', npt
        STOP
      ENDIF
    ENDDO
    CLOSE(99)
  ENDIF 
92  format(5(i6,3x),5(f15.6,3x),i6,3x,100(f15.6,3x))

  ! reset labile C pool,comment out by Q.Zhang 10/09/2011
!  casapool%clabile    = 0.0    
  ! check pool sizes
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

  IF (icycle==1) THEN
    casapool%nplant(:,:) = casapool%cplant(:,:) * casapool%rationcplant(:,:)
    casapool%Nsoil(:,:)  = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
    casapool%Psoil(:,:)  = casapool%ratioPCsoil(:,:) * casapool%Csoil(:,:)
    casapool%Nsoilmin(:) = 2.5
  ENDIF 
    
  IF (icycle >1) THEN
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
  ENDIF 

  IF (icycle >2) THEN
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
  EndIF 

END SUBROUTINE casa_init


SUBROUTINE casa_poolout(ktau,veg,soil,casabiome,casapool,casaflux,casamet, &
                        casabal,phen)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
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
  ! Alfisol	1	61.3
  ! Andisol	2	103.9
  ! Aridisol	3	92.8
  ! Entisol	4	136.9
  ! Gellisol	5	98.2
  ! Histosol	6	107.6
  ! Inceptisol	7	84.1
  ! Mollisol	8	110.1
  ! Oxisol	9	35.4	
  ! Spodosol	10	41.0	
  ! Ultisol	11	51.5	
  ! Vertisol	12	190.6
  DATA psorder/61.3,103.9,92.8,136.9,98.2,107.6,84.1,110.1,35.4,41.0,51.5,190.6/
  DATA pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/
  DATA fracpLab/0.08,0.08,0.10,0.02,0.08,0.08,0.08,0.06,0.02,0.05,0.09,0.05/
  DATA fracPsorb/0.32,0.37,0.57,0.67,0.37,0.37,0.37,0.32,0.24,0.22,0.21,0.38/
  DATA fracPocc/0.36,0.38,0.25,0.26,0.38,0.38,0.38,0.44,0.38,0.38,0.37,0.45/
  DATA fracPorg/0.25,0.17,0.08,0.05,0.17,0.17,0.17,0.18,0.36,0.35,0.34,0.12/
  DATA xpsoil50/7.6,4.1,4.2,3.4,4.1,4.1,4.8,4.1,6.9,6.9,6.9,1.7/

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

    IF (icycle<2) THEN
      casapool%nplant(npt,:) = casapool%rationcplant(npt,:)  &
                             * casapool%cplant(npt,:)
      casapool%nlitter(npt,:)= casapool%rationclitter(npt,:) &
                             * casapool%clitter(npt,:)
      casapool%nsoil(npt,:)  = casapool%ratioNCsoil(npt,:)   &
                             * casapool%Csoil(npt,:)
      casapool%nsoilmin(npt) = 2.0
      casabal%sumnbal(npt)   = 0.0 
    ENDIF 

    IF (icycle<3) THEN
      casabal%sumpbal(npt)   = 0.0
      casapool%pplant(npt,:) = casapool%ratiopcplant(npt,:)  &
                             * casapool%cplant(npt,:)
      casapool%plitter(npt,:)= casapool%ratiopclitter(npt,:) &
                             * casapool%clitter(npt,:)
      casapool%psoil(npt,:)  = casapool%ratioPCsoil(npt,:)   &
                             * casapool%Csoil(npt,:)
      casapool%psoillab(npt) = totpsoil(npt) *fracpLab(nso)
      casapool%psoilsorb(npt)= casaflux%psorbmax(npt) * casapool%psoillab(npt) &
                                /(casaflux%kmlabp(npt)+casapool%psoillab(npt))
      casapool%psoilocc(npt) = totpsoil(npt) *fracPocc(nso)
    ENDIF 

    WRITE(nout,92) ktau,npt,veg%iveg(npt),soil%isoilm(npt),     &
        casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
        casamet%areacell(npt)*(1.0e-9),casamet%glai(npt),       &
        casabiome%sla(veg%iveg(npt)), phen%phase(npt), casapool%clabile(npt), &
        casapool%cplant(npt,:),casapool%clitter(npt,:),casapool%csoil(npt,:), &
        casapool%nplant(npt,:),casapool%nlitter(npt,:),casapool%nsoil(npt,:), &
        casapool%nsoilmin(npt),casapool%pplant(npt,:),          &
        casapool%plitter(npt,:), casapool%psoil(npt,:),         &
        casapool%psoillab(npt),casapool%psoilsorb(npt),casapool%psoilocc(npt), &
        casabal%sumcbal(npt),casabal%sumnbal(npt),casabal%sumpbal(npt)
  ENDDO

  CLOSE(nout)

92    format(5(i6,',',2x),5(f15.6,',',2x),i6,',',2x,100(f15.6,',',2x))
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
        WRITE(nout,92) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
            casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
            casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt),  &
            casabal%Fcnppyear(npt),  &
            casabal%Fcrsyear(npt),casabal%Fcneeyear(npt)  ! ,           &
!            clitterinput(npt,:),csoilinput(npt,:)

      CASE(2)
        WRITE(nout,92) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
            casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
            casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt),  &
        casabal%FCnppyear(npt),casabal%FCrsyear(npt), casabal%FCneeyear(npt), &
!        clitterinput(npt,:),csoilinput(npt,:), &
        casabal%FNdepyear(npt),casabal%FNfixyear(npt),casabal%FNsnetyear(npt), &
        casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt)

      CASE(3)
        WRITE(nout,92) myear,npt,veg%iveg(npt),soil%isoilm(npt), &
        casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt),  &
        casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt), &
       casabal%FCnppyear(npt),casabal%FCrsyear(npt),   casabal%FCneeyear(npt),&
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
SUBROUTINE casa_cnpflux(casaflux,casabal)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  IMPLICIT NONE
  TYPE (casa_flux),    INTENT(INOUT) :: casaflux
  TYPE (casa_balance), INTENT(INOUT) :: casabal
!  REAL(r_2), INTENT(INOUT) :: clitterinput(mp,3),csoilinput(mp,3)
  INTEGER n

  casabal%FCgppyear = casabal%FCgppyear + casaflux%Cgpp   * deltpool
  casabal%FCrpyear  = casabal%FCrpyear  + casaflux%Crp    * deltpool
  casabal%FCnppyear = casabal%FCnppyear + casaflux%Cnpp   * deltpool
  casabal%FCrsyear  = casabal%FCrsyear  + casaflux%Crsoil * deltpool
  casabal%FCneeyear = casabal%FCneeyear &
                    + (casaflux%Cnpp-casaflux%Crsoil) * deltpool
 
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

END SUBROUTINE casa_cnpflux

SUBROUTINE biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
                    casamet,casabal,phen)
  USE cable_def_types_mod
  USE casadimension
  USE casa_cnp_module
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: ktau
  REAL,    INTENT(IN)    :: dels
  INTEGER, INTENT(IN)    :: idoy
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters  
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_balance),          INTENT(INOUT) :: casabal
  TYPE (phen_variable),         INTENT(INOUT) :: phen

  ! local variables
  REAL(r_2),    DIMENSION(mp) :: xnplimit,xNPuptake
  REAL(r_2),    DIMENSION(mp) :: xklitter,xksoil,xkNlimiting
  REAL(r_2),    DIMENSION(mp) :: xkleafcold,xkleafdry,xkleaf
  INTEGER  npt,j

  xKNlimiting = 1.0
  call phenology(idoy,veg,phen)
  call avgsoil(veg,soil,casamet)
  call casa_rplant(veg,casabiome,casapool,casaflux,casamet)

  call casa_allocation(veg,soil,casabiome,casaflux,casamet,phen)

  call casa_xrateplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome, &
                       casamet,phen)
  call casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
                       casaflux,casamet)

  call casa_xnp(xnplimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)

  call casa_xratesoil(xklitter,xksoil,veg,soil,casamet)
  call casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casaflux,casamet)

  IF (icycle>1) THEN
    call casa_xkN(xkNlimiting,casapool,casaflux,casamet,veg)
    DO j=1,mlitter
      casaflux%klitter(:,j) = casaflux%klitter(:,j)* xkNlimiting(:)
    ENDDO
    call casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
    IF (icycle >2) call casa_puptake(veg,xkNlimiting,casabiome, &
                                     casapool,casaflux,casamet)
  ENDIF 

  call casa_delplant(veg,casabiome,casapool,casaflux,casamet)

  call casa_delsoil(veg,casapool,casaflux,casamet)

  call casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet)

  IF (icycle==1) call casa_ndummy(casapool)

  call casa_cnpbal(casapool,casaflux,casabal)
  call casa_cnpflux(casaflux,casabal)


END SUBROUTINE biogeochem

