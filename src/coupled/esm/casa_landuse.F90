module landuse_mod

contains

SUBROUTINE newlitter( casabiome,frac_x,ifpre_x,frac_y,ifpre_y, & 
                      cplant_x,nplant_x,pplant_x,cplant_y,nplant_y,pplant_y, &
                      clitter_x,nlitter_x,plitter_x,clitter_y,nlitter_y,plitter_y )
! Used for LAND USE CHANGE SIMULATION
! Call by casa_reinit
! Transfer the deforest C to litter, and re-allocate litter pools. 
! Q.Zhang @ 29/05/2011
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable

  implicit none

  TYPE (casa_biome),        INTENT(IN) :: casabiome
  logical,DIMENSION(mvtype),INTENT(in) :: ifpre_x,ifpre_y
  real,DIMENSION(mvtype),INTENT(in) :: frac_x,frac_y
  real(r_2),DIMENSION(mvtype,mplant),INTENT(in) :: cplant_x,nplant_x,pplant_x
  real(r_2),DIMENSION(mvtype,mlitter),INTENT(in) :: clitter_x,nlitter_x,plitter_x
  real(r_2),DIMENSION(mvtype,mplant),INTENT(inout) :: cplant_y,nplant_y,pplant_y
  real(r_2),DIMENSION(mvtype,mlitter),INTENT(inout) :: clitter_y,nlitter_y,plitter_y 

  ! local variable
  real(r_2),DIMENSION(mvtype,mlitter,mplant) :: fromPtoL
  real(r_2),DIMENSION(mvtype,mplant) :: dcplant,dnplant,dpplant,ratioLignintoN
  real(r_2),DIMENSION(mvtype,mlitter) :: dclitter,dnlitter,dplitter,clitter_g,nlitter_g,plitter_g
  real(r_2),DIMENSION(mlitter) :: dcY, dnY, dpY
  real,DIMENSION(mvtype) :: dfrac
  real :: sum_y                   
  integer nL, nP, nv
  
  dcplant = 0.
  dnplant = 0.
  dpplant = 0.
  ratioLignintoN = 0.
  fromPtoL = 0.
  dclitter = 0.
  dnlitter = 0.
  dplitter = 0.
  dcY = 0.
  dnY = 0.
  dpY = 0.
  sum_y=0.
  dfrac = frac_y - frac_x

! I. transfer removed plant to litter
  DO nP =1,mplant
                    dcplant(:,nP) = cplant_x(:,nP) * frac_x(:) - cplant_y(:,nP) * frac_y(:)
    IF (icycle > 1) dnplant(:,nP) = nplant_x(:,nP) * frac_x(:) - nplant_y(:,nP) * frac_y(:)
    IF (icycle > 2) dpplant(:,nP) = pplant_x(:,nP) * frac_x(:) - pplant_y(:,nP) * frac_y(:)
  END DO
! NB: logged wood should not be transfered to litter
                  dcplant(1:mlogmax,wood) = 0.
  IF (icycle > 1) dnplant(1:mlogmax,wood) = 0.
  IF (icycle > 2) dpplant(1:mlogmax,wood) = 0.

  WHERE(sum(dcplant,2) > 0.)
  ! In land use, all plant nutient is allocated to litter pools without re-asorbsion.Q.Zhang 11/08/2011 
    ratioLignintoN(:,leaf) = cplant_x(:,leaf)/max(1.0e-10,nplant_x(:,leaf)) &
                             * casabiome%fracLigninplant(:,leaf)
    ratioLignintoN(:,froot)= cplant_x(:,froot)/max(1.0e-10,nplant_x(:,froot)) &
                             * casabiome%fracLigninplant(:,froot)

    fromPtoL(:,metb,leaf)   = max(0.001, 0.85 - 0.018 *ratioLignintoN(:,leaf))
    fromPtoL(:,metb,froot)  = max(0.001, 0.85 - 0.018 *ratioLignintoN(:,froot))
    fromPtoL(:,str,leaf)    = 1.0 - fromPtoL(:,metb,leaf)
    fromPtoL(:,str,froot)   = 1.0 - fromPtoL(:,metb,froot)
    fromPtoL(:,cwd,wood)    = 1.0
  ENDWHERE

  DO nv=1, mvtype
  ! transfer removed C,N,P pools from plant to litter 
    IF(ifpre_x(nv) .and. frac_x(nv)>frac_y(nv))THEN

      DO nL=1,mlitter
        DO nP=1,mplant
           dclitter(nv,nL) = dclitter(nv,nL)+ fromPtoL(nv,nL,nP) * dcplant(nv,nP)
        ENDDO
      ENDDO
  
      IF(icycle > 1) THEN
         dnlitter(nv,str) = (fromPtoL(nv,str,leaf) * dcplant(nv,leaf) &
                           + fromPtoL(nv,str,froot) * dcplant(nv,froot)) * ratioNCstrfix
         dnlitter(nv,metb) = dnplant(nv,leaf) + dnplant(nv,froot) - dnlitter(nv,str)
         dnlitter(nv,CWD) = dnplant(nv,wood)
      ENDIF !end "icycle >1"
    
      IF(icycle > 2) THEN
         dplitter(nv,str) = (fromPtoL(nv,str,leaf) * dcplant(nv,leaf) &
                           + fromPtoL(nv,str,froot)* dcplant(nv,froot)) * ratioPCstrfix
         dplitter(nv,metb) = dpplant(nv,leaf) + dpplant(nv,froot) -dplitter(nv,str)
         dplitter(nv,CWD) = dpplant(nv,wood)
      ENDIF  !of "icycle >2"
    ENDIF
  END DO

! II. re-allocate litter pools according to patch weights.
! average pool variables from gridcell to new patches.
  DO nv=1,mvtype
    IF (ifpre_x(nv) .and. dfrac(nv)<0.0) THEN
                      dcY(:) = dcY(:) + dclitter(nv,:) + clitter_x(nv,:)*abs(dfrac(nv))
      IF (icycle > 1) dnY(:) = dnY(:) + dnlitter(nv,:) + nlitter_x(nv,:)*abs(dfrac(nv))
      IF (icycle > 2) dpY(:) = dpY(:) + dplitter(nv,:) + plitter_x(nv,:)*abs(dfrac(nv))
    ENDIF
  END DO

  DO nv=1,mvtype
     IF ((frac_x(nv)-frac_y(nv))<0.) THEN
       sum_y = sum_y + abs(dfrac(nv)) 
     ENDIF
  END DO

  DO nv=1,mvtype
    IF (ifpre_y(nv)) THEN   ! pft exist in the 2nd year
      IF ((frac_x(nv)-frac_y(nv))>0.) THEN  ! patch weight decrease 
                        clitter_y(nv,:) = clitter_x(nv,:)
        IF (icycle > 1) nlitter_y(nv,:) = nlitter_x(nv,:)
        IF (icycle > 2) plitter_y(nv,:) = plitter_x(nv,:)
      ELSE IF ((frac_x(nv)-frac_y(nv))<0.) THEN ! patch increase
                        clitter_y(nv,:) = (clitter_x(nv,:)*frac_x(nv) + dcY(:)*dfrac(nv)/sum_y) / frac_y(nv)
        IF (icycle > 1) nlitter_y(nv,:) = (nlitter_x(nv,:)*frac_x(nv) + dnY(:)*dfrac(nv)/sum_y)	/ frac_y(nv)
        IF (icycle > 2) plitter_y(nv,:) = (plitter_x(nv,:)*frac_x(nv) + dpY(:)*dfrac(nv)/sum_y)	/ frac_y(nv)
      ELSE ! no change 
                        clitter_y(nv,:) = clitter_x(nv,:)
        IF (icycle > 1) nlitter_y(nv,:) = nlitter_x(nv,:)
	IF (icycle > 2) plitter_y(nv,:) = plitter_x(nv,:)
      ENDIF
    ENDIF
  END DO

END SUBROUTINE newlitter


SUBROUTINE newlitter_thin( &
      casabiome, &
      tile_exists, &
      cplant_x, &
      nplant_x, &
      pplant_x, &
      cplant_y, &
      nplant_y, &
      pplant_y, &
      clitter, &
      nlitter, &
      plitter, &
      thinning)
   !* Transfer the thinned forest leaf and root biomass to litter pools.
   !
   ! ## Procedure
   !
   ! 1. Find the difference between the plant biomass pools before and after
   !    thinning.
   ! 2. Ignore wood biomass (because it's already transferred to wood harvest)
   ! 3. Calculate the plant to litter ratio matrix based on C/N
   ! 4. Transfer biomass to litter pools
   USE cable_def_types_mod
   USE casadimension
   USE casaparm
   USE casavariable

   IMPLICIT NONE

   TYPE (casa_biome), INTENT (IN) :: casabiome
   LOGICAL, INTENT (IN) :: tile_exists(mvtype)
   REAL, INTENT (IN) :: thinning(mvtype) !! Thinning fraction (1 = no thinning)
   REAL (r_2), INTENT (IN) :: cplant_x(mvtype,mplant) !! Plant C before thinning
   REAL (r_2), INTENT (IN) :: nplant_x(mvtype,mplant) !! Plant N before thinning
   REAL (r_2), INTENT (IN) :: pplant_x(mvtype,mplant) !! Plant P before thinning
   REAL (r_2), INTENT (IN) :: cplant_y(mvtype,mplant) !! Plant C after thinning
   REAL (r_2), INTENT (IN) :: nplant_y(mvtype,mplant) !! Plant N after thinning
   REAL (r_2), INTENT (IN) :: pplant_y(mvtype,mplant) !! Plant P after thinning
   REAL (r_2), INTENT (INOUT) :: clitter(mvtype,mlitter) !! Litter C
   REAL (r_2), INTENT (INOUT) :: nlitter(mvtype,mlitter) !! Litter N
   REAL (r_2), INTENT (INOUT) :: plitter(mvtype,mlitter) !! Litter P

   ! Local variables
   REAL (r_2) :: fromPtoL(mvtype,mlitter,mplant)
   REAL (r_2) :: dcplant(mvtype,mplant)
   REAL (r_2) :: dnplant(mvtype,mplant)
   REAL (r_2) :: dpplant(mvtype,mplant)
   REAL (r_2) :: ratioLignintoN(mvtype,mplant)
   REAL (r_2) :: dclitter(mvtype,mlitter)
   REAL (r_2) :: dnlitter(mvtype,mlitter)
   REAL (r_2) :: dplitter(mvtype,mlitter)
   REAL (r_2) :: imbalance
   INTEGER :: nl, np, nv

   dcplant = 0.0
   dnplant = 0.0
   dpplant = 0.0
   ratioLignintoN = 0.0
   fromPtoL = 0.0
   dclitter = 0.0
   dnlitter = 0.0
   dplitter = 0.0

   ! Find the change in the plant pools
   DO np=1,mplant
                      dcplant(:,np) = cplant_x(:,np) - cplant_y(:,np)
      IF (icycle > 1) dnplant(:,np) = nplant_x(:,np) - nplant_y(:,np)
      IF (icycle > 2) dpplant(:,np) = pplant_x(:,np) - pplant_y(:,np)
   END DO

   ! Wood should not be transfered to litter, it has already gone to products
                   dcplant(:,wood) = 0.0
   IF (icycle > 1) dnplant(:,wood) = 0.0
   IF (icycle > 2) dpplant(:,wood) = 0.0

   ! Calculate plant->litter allocation ratios from C/N.
   ! All plant nutients are allocated to litter pools without re-asorpsion.
   WHERE (SUM(dcplant, 2) > 0.0)
      ratioLignintoN(:,leaf) = &
            cplant_x(:,leaf)/MAX(1.0e-10, nplant_x(:,leaf)) &
            *casabiome%fracLigninplant(:,leaf)
      ratioLignintoN(:,froot) = &
            cplant_x(:,froot)/MAX(1.0e-10, nplant_x(:,froot)) &
            *casabiome%fracLigninplant(:,froot)

      fromPtoL(:,metb,leaf) = MAX(0.001, 0.85 - 0.018*ratioLignintoN(:,leaf))
      fromPtoL(:,metb,froot) = MAX(0.001, 0.85 - 0.018*ratioLignintoN(:,froot))
      fromPtoL(:,str,leaf) = 1.0 - fromPtoL(:,metb,leaf)
      fromPtoL(:,str,froot) = 1.0 - fromPtoL(:,metb,froot)
      fromPtoL(:,cwd,wood) = 0.0
   END WHERE

   DO nv=1,mlogmax
      IF (tile_exists(nv) .AND. thinning(nv)<1.0) THEN
         ! Caluclate the change in each litter pools.
         DO nl=1,mlitter
            DO np=1,mplant
               dclitter(nv,nl) = &
                     dclitter(nv,nl) + fromPtoL(nv,nl,np)*dcplant(nv,np)
            END DO
         END DO

         IF (icycle > 1) THEN
            dnlitter(nv,str) = &
                    (fromPtoL(nv,str,leaf)*dcplant(nv,leaf) &
                  + fromPtoL(nv,str,froot)*dcplant(nv,froot)) &
                  *ratioNCstrfix
            dnlitter(nv,metb) = &
                    dnplant(nv,leaf) &
                  + dnplant(nv,froot) &
                  - dnlitter(nv,str)
            dnlitter(nv,CWD) = dnplant(nv,wood)
         END IF

         IF (icycle > 2) THEN
            dplitter(nv,str) = &
                    (fromPtoL(nv,str,leaf)*dcplant(nv,leaf) &
                  + fromPtoL(nv,str,froot)*dcplant(nv,froot)) &
                  *ratioPCstrfix
            dplitter(nv,metb) = &
                    dpplant(nv,leaf) &
                  + dpplant(nv,froot) &
                  - dplitter(nv,str)
            dplitter(nv,CWD) = dpplant(nv,wood)
         END IF

         ! Modify the litter pools.
                         clitter(nv,:) = clitter(nv,:) + dclitter(nv,:)
         IF (icycle > 1) nlitter(nv,:) = nlitter(nv,:) + dnlitter(nv,:)
         IF (icycle > 2) plitter(nv,:) = plitter(nv,:) + dplitter(nv,:)
      END IF
   END DO

   ! Check for conservation of mass
   imbalance = ABS(SUM(dcplant(1:mlogmax,:)) - SUM(dclitter(1:mlogmax,:)))
   IF (imbalance > 1.0E-10) THEN
      WRITE (6,*) "Violation of carbon conservation in newlitter_thin"
      WRITE (6,*) "difference", &
            SUM(dcplant(1:mlogmax,:)) - SUM(dclitter(1:mlogmax,:))
      WRITE (6,*) "dcplant", dcplant
      WRITE (6,*) "dclitter", dclitter
      WRITE (6,*) "tile_exists", tile_exists
   END IF
END SUBROUTINE newlitter_thin


SUBROUTINE newplant(cplant_x,frac_x,ifpre_x, &
                    cplant_y,frac_y,ifpre_y,logc)
! Used for LAND USE CHANGE SIMULATION
! Call by casa_reinit
! Re-allcate plant C,N,P pools to new patch array. 
! Q.Zhang @ 29/05/2011
  USE cable_def_types_mod
  USE casadimension
  USE casaparm

  implicit none

  logical,DIMENSION(mvtype),INTENT(in) :: ifpre_x,ifpre_y
  real,DIMENSION(mvtype),INTENT(in) :: frac_x,frac_y
  real(r_2),DIMENSION(mvtype,mplant),INTENT(inout) :: cplant_x,cplant_y
  real(r_2),DIMENSION(mvtype),INTENT(inout) :: logc
  ! local variable
  integer p 

  DO p = 1, mvtype
     ! exist in both years    
     IF (ifpre_x(p) .and. ifpre_y(p)) THEN
        IF (abs(frac_x(p)-0.0)<1.e-8 .or. abs(frac_y(p)-0.0)<1.e-8) THEN
print *, 'Lest Veg0', p,ifpre_x(p),ifpre_y(p),frac_x(p),frac_y(p)
          STOP "vegetation fraction .eq. 0"
        END IF
        IF ((frac_x(p)-frac_y(p))>0.) THEN  ! patch weight decrease 
          ! New pools
          cplant_y(p,:) = cplant_x(p,:)
          ! Save wood log
          IF (p<=mlogmax) logc(p) = cplant_x(p,wood) * (frac_x(p) - frac_y(p))
        ELSE ! patch weight incease
          cplant_y(p,:) = cplant_x(p,:)*frac_x(p)/frac_y(p)
        END IF
   ! plant clear in the second year
     ELSEIF (ifpre_x(p) .and. .not.ifpre_y(p)) THEN
        IF (p<=mlogmax) logc(p) = cplant_x(p,wood) * (frac_x(p) - frac_y(p))
    ! does not exist in both years
     ELSE 

     END IF

  END DO  ! end pft loop

END SUBROUTINE newplant


SUBROUTINE newsoil(nd,csoil_x,frac_x,ifpre_x,csoil_y,frac_y,ifpre_y)
! Used for LAND USE CHANGE SIMULATION
! Call by casa_reinit
! Re-allocate soil C,N and P pools
! Q.Zhang @ 29/05/2011
! L.stevens @ 19/01/2018
  USE cable_def_types_mod
  USE casadimension
  USE casaparm

  implicit none

  integer,INTENT(in) :: nd  ! dimension of soil pool
  logical,DIMENSION(mvtype),INTENT(in) :: ifpre_x,ifpre_y
  real,DIMENSION(mvtype),INTENT(in) :: frac_x,frac_y
  real(r_2),DIMENSION(mvtype,nd),INTENT(inout) :: csoil_x,csoil_y
! local variable
  real,DIMENSION(nd)     :: tmpVar
  real                   :: Rcount
  real,DIMENSION(mvtype) :: dfrac
  integer nv

  dfrac = frac_y - frac_x ! current minus previous
  tmpVar = 0.0
  Rcount = 0.0

  DO nv=1,mvtype
    IF (dfrac(nv)<0.0) THEN
      tmpVar = tmpVar + ABS(dfrac(nv))*csoil_x(nv,:)
      Rcount = Rcount + ABS(dfrac(nv))
    ENDIF
  END DO
  if (Rcount > 0) tmpVar = tmpVar/Rcount
  
  DO nv=1,mvtype
    IF (dfrac(nv) > 0.0) THEN   ! pft exist in the 2nd year
      csoil_y(nv,:) = (dfrac(nv)*tmpVar + &
                 csoil_x(nv,:) * frac_x(nv)) / frac_y(nv)
    ELSE
      csoil_y(nv,:) = csoil_x(nv,:)
    ENDIF
  END DO

END SUBROUTINE newsoil

End module landuse_mod
