MODULE BLAZE_MPI

  USE MPI
  USE cable_mpicommon
  USE cable_def_types_mod, ONLY: ncp
  USE BLAZE_MOD
  USE SIMFIRE_MOD

  ! Total number of restart parameters for BLAZE
  INTEGER, PARAMETER :: n_blaze_restart = 9

  ! Total number of output parameters for BLAZE
  ! for BLAZE%OUTMODE == "std"
  INTEGER, PARAMETER :: n_blaze_output_std   = 10
  ! add for BLAZE%OUTMODE == "full"
  INTEGER, PARAMETER :: n_blaze_output_extra = 13

  ! Total number of restart parameters for SIMFIRE
  INTEGER, PARAMETER :: n_simfire_restart = 4

  ! Total number of input parameters for SIMFIRE
  INTEGER, PARAMETER :: n_simfire_input = 1

  ! Total number of output parameters for SIMFIRE
  INTEGER, PARAMETER :: n_simfire_output = 4

  ! Constant number of days = 366 
  INTEGER, PARAMETER :: ndoy = 366
  
CONTAINS

SUBROUTINE master_blaze_types (comm, wland, mp, BLAZE, blaze_restart_ts, blaze_out_ts)
  
  ! Send blaze restart data to workers  
  
  IMPLICIT NONE
  
  INTEGER              , INTENT(IN)  :: comm ! MPI communicator to talk to the workers
  INTEGER              , INTENT(IN)  :: mp   ! Number of gridcells
  TYPE(TYPE_BLAZE)     , INTENT(IN)  :: BLAZE
  TYPE(lpdecomp_t), DIMENSION(:) , INTENT(IN)  :: wland
  INTEGER, DIMENSION(:), INTENT(OUT) :: blaze_restart_ts,blaze_out_ts 

  ! MPI: temp arrays for marshalling all types into a struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types
  INTEGER :: ntyp ! number of worker's types

  INTEGER :: last2d, i

  ! MPI: block lenghts for hindexed representing all vectors
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blen

  ! MPI: block lengths and strides for hvector representing matrices
  INTEGER :: r1len, r2len, i1len
  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride, istride

  INTEGER :: tsize, totalrecv, totalsend
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

  INTEGER :: rank, off, cnt, wnp
  INTEGER :: bidx, midx, vidx, ierr

  ! Restart value handles (restart_blaze_ts()) 
  
  ntyp = n_blaze_restart

  ALLOCATE (blocks(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))

  istride  = mp * extid ! short integer 
  r1stride = mp * extr1 ! single precision
  r2stride = mp * extr2 ! double precision

  DO rank = 1, wnp
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     i1len = cnt * extid  
     r1len = cnt * extr1
     r2len = cnt * extr2

     bidx = 0

     ! ------------- 2D arrays -------------

     ! Annual (daily) rainfall (ncells,366)
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AnnRainf(off,1), displs(bidx), ierr) ! 1
     CALL MPI_Type_create_hvector (ndoy, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1

!CLN     ! Above ground life woody biomass 
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (BLAZE%AGLB_w(off,1), displs(bidx), ierr) ! 2
!CLN     CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
!CLN     &                             types(bidx), ierr)
!CLN     blocks(bidx) = 1
!CLN
!CLN     ! Above ground life grassy biomass 
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (BLAZE%AGLB_g(off,1), displs(bidx), ierr) ! 3
!CLN     CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
!CLN     &                             types(bidx), ierr)
!CLN     blocks(bidx) = 1

     ! Above ground woody litter
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AGLit_w(off,1), displs(bidx), ierr) ! 4
     CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1

     ! Above ground grassy litter 
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AGLit_g(off,1), displs(bidx), ierr) ! 5
     CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1

     last2d = bidx
     
     ! ------------- 1D vectors -------------

     ! Integer days since last rainfall
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%DSLR(off), displs(bidx), ierr)
     blocks(bidx) = i1len

     ! Real(sp) Last rainfall
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%LR(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! current KBDI
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%KBDI(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! DEADWOOD
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%DEADWOOD(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! ------------- Wrap up -------------

     types(last2d+1:bidx) = MPI_BYTE

     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE (*,*) 'invalid blaze ntyp constant, fix it!'
        CALL MPI_Abort (comm, 1, ierr)
     END IF

     CALL MPI_Type_create_struct (bidx, blocks, displs, types, blaze_restart_ts(rank), ierr)
     CALL MPI_Type_commit (blaze_restart_ts(rank), ierr)

     CALL MPI_Type_size (blaze_restart_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent (blaze_restart_ts(rank), tmplb, text, ierr)

     WRITE (*,*) 'restart results recv from worker, size, extent, lb: ', &
   &       rank,tsize,text,tmplb

     totalrecv = totalrecv + tsize

     ! free the partial types used for matrices
     DO i = 1, last2d
        CALL MPI_Type_free (types(i), ierr)
     END DO

  END DO

  WRITE (*,*) 'total size of restart fields received from all workers: ', totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  CALL MPI_Reduce (MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, &
    &     0, comm, ierr)

  WRITE (*,*) 'total size of restart fields sent by all workers: ', totalsend

  IF (totalrecv /= totalsend) THEN
          WRITE (*,*) 'error: restart fields totalsend and totalrecv differ'
          CALL MPI_Abort (comm, 0, ierr)
  END IF

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  !=============================================================================
  ! Standard desired Output (blaze_ts()) 
  !=============================================================================

  ntyp = n_blaze_output_std
  IF( TRIM(BLAZE%OUTMODE) == "full") &
       ntyp = ntyp + n_blaze_output_extra
  
  ALLOCATE (blocks(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))
  
  istride  = mp * extid ! short integer 
  r1stride = mp * extr1 ! single precision
  r2stride = mp * extr2 ! double precision

  DO rank = 1, wnp
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     i1len = cnt * extid  
     r1len = cnt * extr1
     r2len = cnt * extr2

     bidx   = 0
     last2d = 0
     
     ! ------------- 2D arrays -------------

     IF ( TRIM(BLAZE%OUTMODE) == "full" ) THEN 
        
        ! Annual (daily) rainfall (ncells,366)
        bidx = bidx + 1
        CALL MPI_Get_address (BLAZE%TO(off,1), displs(bidx), ierr) ! 1
        CALL MPI_Type_create_hvector (NTO, r1len, r1stride, MPI_BYTE, &
             &                             types(bidx), ierr)
        blocks(bidx) = 1
        
!CLN        ! Above ground life woody biomass 
!CLN        bidx = bidx + 1
!CLN        CALL MPI_Get_address (BLAZE%AGLB_w(off,1), displs(bidx), ierr) ! 2
!CLN        CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
!CLN             &                             types(bidx), ierr)
!CLN        blocks(bidx) = 1
!CLN        
!CLN        ! Above ground life grassy biomass 
!CLN        bidx = bidx + 1
!CLN        CALL MPI_Get_address (BLAZE%AGLB_g(off,1), displs(bidx), ierr) ! 3
!CLN        CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
!CLN             &                             types(bidx), ierr)
!CLN        blocks(bidx) = 1
        
        ! Above ground woody litter
        bidx = bidx + 1
        CALL MPI_Get_address (BLAZE%AGLit_w(off,1), displs(bidx), ierr) ! 4
        CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
             &                             types(bidx), ierr)
        blocks(bidx) = 1
        
        ! Above ground grassy litter 
        bidx = bidx + 1
        CALL MPI_Get_address (BLAZE%AGLit_g(off,1), displs(bidx), ierr) ! 5
        CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
             &                             types(bidx), ierr)
        blocks(bidx) = 1

     ENDIF
     
     last2d = bidx
     
     ! ------------- 1D vectors -------------

     IF ( TRIM(BLAZE%OUTMODE) == "full" ) THEN 

        ! Integer days since last rainfall
        bidx = bidx + 1
        CALL MPI_Get_address (BLAZE%DSLR(off), displs(bidx), ierr)
        blocks(bidx) = i1len
        
        ! Real(sp) Last rainfall
        bidx = bidx + 1
        CALL MPI_Get_address (BLAZE%LR(off), displs(bidx), ierr)
        blocks(bidx) = r1len
        
        ! RAINF
        bidx = bidx + 1
        CALL MPI_Get_address (BLAZE%RAINF(off), displs(bidx), ierr)
        blocks(bidx) = r1len
        
        ! U10
        bidx = bidx + 1
        CALL MPI_Get_address (BLAZE%U10(off), displs(bidx), ierr)
        blocks(bidx) = r1len
        
        ! RH
        bidx = bidx + 1
        CALL MPI_Get_address (BLAZE%RH(off), displs(bidx), ierr)
        blocks(bidx) = r1len
        
        ! TMAX
        bidx = bidx + 1
        CALL MPI_Get_address (BLAZE%TMAX(off), displs(bidx), ierr)
        blocks(bidx) = r1len
        
        ! TMIN
        bidx = bidx + 1
        CALL MPI_Get_address (BLAZE%TMIN(off), displs(bidx), ierr)
        blocks(bidx) = r1len
        
        ! DEADWOOD
        bidx = bidx + 1
        CALL MPI_Get_address (BLAZE%DEADWOOD(off), displs(bidx), ierr)
        blocks(bidx) = r1len

     END IF
     
     ! current KBDI
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%KBDI(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! FLIx
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%FLIx(off), displs(bidx), ierr)
     blocks(bidx) = i1len

     ! MacArthur Forest Fire Dangeer index 
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%FFDI(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! FLI at burntime
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%FLI(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! daily FLI
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%DFLI(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! AB
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AB(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! available fuel w (CLN check unuit!!!)
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%w(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! MacArthur drought-factor D
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%D(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! Flame height Z
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%Z(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! Rate-of-Spread ROS
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%ROS(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     
     ! ------------- Wrap up -------------

     types(last2d+1:bidx) = MPI_BYTE

     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE (*,*) 'invalid blaze ntyp constant, fix it!'
        CALL MPI_Abort (comm, 1, ierr)
     END IF

     CALL MPI_Type_create_struct (bidx, blocks, displs, types, blaze_out_ts(rank), ierr)
     CALL MPI_Type_commit (blaze_out_ts(rank), ierr)

     CALL MPI_Type_size (blaze_out_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent (blaze_out_ts(rank), tmplb, text, ierr)

     WRITE (*,*) 'restart results recv from worker, size, extent, lb: ', &
   &       rank,tsize,text,tmplb

     totalrecv = totalrecv + tsize

     ! free the partial types used for matrices
     DO i = 1, last2d
        CALL MPI_Type_free (types(i), ierr)
     END DO

  END DO

  WRITE (*,*) 'total size of restart fields received from all workers: ', totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  CALL MPI_Reduce (MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, &
    &     0, comm, ierr)

  WRITE (*,*) 'total size of restart fields sent by all workers: ', totalsend

  IF (totalrecv /= totalsend) THEN
          WRITE (*,*) 'error: restart fields totalsend and totalrecv differ'
          CALL MPI_Abort (comm, 0, ierr)
  END IF

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)


END SUBROUTINE master_blaze_types

SUBROUTINE worker_blaze_types(comm, mp, BLAZE, blaze_restart_t, blaze_out_t)

  INTEGER         , INTENT(IN)  :: comm ! MPI communicator to talk to the master
  INTEGER         , INTENT(IN)  :: mp   ! Number of gridcells
  TYPE(TYPE_BLAZE), INTENT(IN)  :: BLAZE
  INTEGER         , INTENT(OUT) :: blaze_restart_t,blaze_out_t 
  
  ! MPI: temp arrays for marshalling all types into a struct
  INTEGER, ALLOCATABLE, DIMENSION(:)                        :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:)                        :: types
  INTEGER :: ntyp ! number of worker's types
  
  ! MPI: block lengths and strides for hvector representing matrices
  INTEGER :: r1len, r2len, i1len, llen
  
  INTEGER :: rank, off, cnt
  INTEGER :: bidx, midx, vidx, ierr, nd, ny

  INTEGER :: tsize
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

   
  CALL MPI_Comm_rank (comm, rank, ierr)

  !=============================================================================
  ! RESTART comm (blaze_restart_t) 
  !=============================================================================
  
  ntyp = n_blaze_restart
  
  ALLOCATE (blocks(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))
  
  cnt  = mp
  
  r1len = cnt * extr1
  r2len = cnt * extr2
  i1len = cnt * extid
  llen  = cnt * extl

  off  = 1
  bidx = 0

  ! ------------- 2D arrays -------------
  
  ! Annual (daily) rainfall (ncells,366)
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%AnnRainf(off,1), displs(bidx), ierr)
  blocks(bidx) = r1len * ndoy

!CLN  ! Above ground life woody biomass 
!CLN  bidx = bidx + 1
!CLN  CALL MPI_Get_address (BLAZE%AGLB_w(off,1), displs(bidx), ierr)
!CLN  blocks(bidx) = r1len * ncp
!CLN
!CLN  ! Above ground life grassy biomass 
!CLN  bidx = bidx + 1
!CLN  CALL MPI_Get_address (BLAZE%AGLB_g(off,1), displs(bidx), ierr)
!CLN  blocks(bidx) = r1len * ncp

  ! Above ground woody litter
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%AGLit_w(off,1), displs(bidx), ierr)
  blocks(bidx) = r1len * ncp

  ! Above ground grassy litter 
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%AGLit_g(off,1), displs(bidx), ierr)
  blocks(bidx) = r1len * ncp

  ! ------------- 1D vectors -------------

  ! Integer days since last rainfall
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%DSLR(off), displs(bidx), ierr)
  blocks(bidx) = i1len

  ! Real(sp) Last rainfall
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%LR(off), displs(bidx), ierr)
  blocks(bidx) = r1len
  
  ! current KBDI
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%KBDI(off), displs(bidx), ierr)
  blocks(bidx) = r1len
  
  ! DEADWOOD
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%DEADWOOD(off), displs(bidx), ierr)
  blocks(bidx) = r1len
  
  ! ------------- Wrap up -------------

  types = MPI_BYTE
 
  ! MPI: sanity check
  IF (bidx /= ntyp) THEN
     WRITE (*,*) 'invalid nrestart constant, fix it!'
     CALL MPI_Abort (comm, 1, ierr)
  END IF

  CALL MPI_Type_create_struct (bidx, blocks, displs, types, blaze_restart_t, ierr)
  CALL MPI_Type_commit (blaze_restart_t, ierr)
  
  CALL MPI_Type_size (blaze_restart_t, tsize, ierr)
  CALL MPI_Type_get_extent (blaze_restart_t, tmplb, text, ierr)
  
  WRITE (*,*) 'restart_blaze struct blocks, size, extent and lb: ',rank,bidx,tsize,text,tmplb
  
  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
  
  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  !=============================================================================
  ! Standard desired Output (blaze_out_t) 
  !=============================================================================

  ntyp = n_blaze_output_std
  IF( TRIM(BLAZE%OUTMODE) == "full") &
       ntyp = ntyp + n_blaze_output_extra
  
  ALLOCATE (blocks(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))
  
  istride  = cnt * extid ! short integer 
  r1stride = cnt * extr1 ! single precision
  r2stride = cnt * extr2 ! double precision

  bidx = 0

  ! ------------- 2D arrays -------------

  IF ( TRIM(BLAZE%OUTMODE) == "full" ) THEN 
  
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%TO(off,1), displs(bidx), ierr)
     blocks(bidx) = r1len * NTO

!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (BLAZE%AGLB_w(off,1), displs(bidx), ierr)
!CLN     blocks(bidx) = r1len * ncp
!CLN
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (BLAZE%AGLB_g(off,1), displs(bidx), ierr)
!CLN     blocks(bidx) = r1len * ncp

     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AGLit_w(off,1), displs(bidx), ierr)
     blocks(bidx) = r1len * ncp

     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AGLit_g(off,1), displs(bidx), ierr)
     blocks(bidx) = r1len * ncp

  ENDIF
     
  ! ------------- 1D vectors -------------
  
  IF ( TRIM(BLAZE%OUTMODE) == "full" ) THEN 
     
     ! Integer days since last rainfall
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%DSLR(off), displs(bidx), ierr)
     blocks(bidx) = i1len
     
     ! Real(sp) Last rainfall
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%LR(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     
     ! RAINF
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%RAINF(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     
     ! U10
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%U10(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     
     ! RH
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%RH(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     
     ! TMAX
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%TMAX(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     
     ! TMIN
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%TMIN(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     
     ! DEADWOOD
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%DEADWOOD(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     
  END IF
  
  ! current KBDI
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%KBDI(off), displs(bidx), ierr)
  blocks(bidx) = r1len
  
  ! FLIx
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%FLIx(off), displs(bidx), ierr)
  blocks(bidx) = i1len
  
  ! MacArthur Forest Fire Dangeer index 
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%FFDI(off), displs(bidx), ierr)
  blocks(bidx) = r1len
  
  ! FLI at burntime
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%FLI(off), displs(bidx), ierr)
  blocks(bidx) = r1len
  
  ! daily FLI
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%DFLI(off), displs(bidx), ierr)
  blocks(bidx) = r1len
  
  ! AB
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%AB(off), displs(bidx), ierr)
  blocks(bidx) = r1len
  
  ! available fuel w (CLN check unuit!!!)
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%w(off), displs(bidx), ierr)
  blocks(bidx) = r1len
  
  ! MacArthur drought-factor D
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%D(off), displs(bidx), ierr)
  blocks(bidx) = r1len
  
  ! Flame height Z
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%Z(off), displs(bidx), ierr)
  blocks(bidx) = r1len
  
  ! Rate-of-Spread ROS
  bidx = bidx + 1
  CALL MPI_Get_address (BLAZE%ROS(off), displs(bidx), ierr)
  blocks(bidx) = r1len
  
  ! ------------- Wrap up -------------
  
  types = MPI_BYTE
  
  CALL MPI_Type_create_struct (bidx, blocks, displs, types, blaze_out_t, ierr)
  CALL MPI_Type_commit (blaze_out_t, ierr)
  
  CALL MPI_Type_size (blaze_out_t, tsize, ierr)
  CALL MPI_Type_get_extent (blaze_out_t, tmplb, text, ierr)
  
  WRITE (*,*) 'restart struct blocks, size, extent and lb: ',rank,bidx,tsize,text,tmplb
  
  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
     
  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)
     
END SUBROUTINE worker_blaze_types

SUBROUTINE master_simfire_types(comm, wland, mp, SF, simfire_restart_ts, simfire_inp_ts,simfire_out_ts)

  ! Send simfire restart data to workers  
  
  IMPLICIT NONE
  
  INTEGER              , INTENT(IN)  :: comm ! MPI communicator to talk to the workers
  TYPE(lpdecomp_t), DIMENSION(:)     , INTENT(IN)  :: wland
  INTEGER              , INTENT(IN)  :: mp   ! Number of gridcells
  TYPE(TYPE_SIMFIRE)   , INTENT(IN)  :: SF
  INTEGER, DIMENSION(:), INTENT(OUT) :: simfire_restart_ts,simfire_inp_ts,simfire_out_ts 

  ! MPI: temp arrays for marshalling all types into a struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types
  INTEGER :: ntyp ! number of worker's types

  INTEGER :: last2d, i

  ! MPI: block lenghts for hindexed representing all vectors
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blen

  ! MPI: block lengths and strides for hvector representing matrices
  INTEGER :: r1len, r2len, i1len
  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride, istride

  INTEGER :: tsize, totalrecv, totalsend
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

  INTEGER :: rank, off, cnt, wnp
  INTEGER :: bidx, midx, vidx, ierr

  ! Restart value handles (simfire_restart_ts) 
  
 !CLN   INTEGER, DIMENSION(:), ALLOCATABLE  :: IGBP, BIOME, REGION, NDAY
 !CLN  REAL,    DIMENSION(:), ALLOCATABLE  :: POPD, MAX_NESTEROV, CNEST, FAPAR, LAT, LON, FLI
 !CLN  INTEGER   :: SYEAR, EYEAR, NCELLS
 !CLN  REAL      :: RES, RESF
 !CLN  CHARACTER :: IGBPFILE*120, HYDEPATH*100
 
  ntyp = n_simfire_restart

  ALLOCATE (blocks(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))

  istride  = mp * extid ! short integer 
  r1stride = mp * extr1 ! single precision
  r2stride = mp * extr2 ! double precision

  DO rank = 1, wnp
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     i1len = cnt * extid  
     r1len = cnt * extr1
     r2len = cnt * extr2

     bidx = 0

     ! ------------- 2D arrays -------------

     ! Last years max_Nesterov Index 
     bidx = bidx + 1
     CALL MPI_Get_address (SF%SAV_NESTEROV(off,1), displs(bidx), ierr) ! 1
     CALL MPI_Type_create_hvector (12, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1

     ! FAPAR
     bidx = bidx + 1
     CALL MPI_Get_address (SF%SAV_FAPAR(off,1), displs(bidx), ierr) ! 1
     CALL MPI_Type_create_hvector (FAPAR_AVG_INT, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1

     last2d = bidx

     ! ------------- 1D arrays -------------

     ! Integer IGBP classification
     bidx = bidx + 1
     CALL MPI_Get_address (SF%IGBP(off), displs(bidx), ierr)
     blocks(bidx) = i1len

     ! BIOME
     bidx = bidx + 1
     CALL MPI_Get_address (SF%BIOME(off), displs(bidx), ierr) 
     blocks(bidx) = i1len

     ! ------------- Wrap up -------------

     types(last2d+1:bidx) = MPI_BYTE

     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE (*,*) 'invalid blaze ntyp constant, fix it!'
        CALL MPI_Abort (comm, 1, ierr)
     END IF

     CALL MPI_Type_create_struct (bidx, blocks, displs, types, simfire_restart_ts(rank), ierr)
     CALL MPI_Type_commit (simfire_restart_ts(rank), ierr)

     CALL MPI_Type_size (simfire_restart_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent (simfire_restart_ts(rank), tmplb, text, ierr)

     WRITE (*,*) 'restart results recv from worker, size, extent, lb: ', &
   &       rank,tsize,text,tmplb

     totalrecv = totalrecv + tsize

     ! free the partial types used for matrices
     DO i = 1, last2d
        CALL MPI_Type_free (types(i), ierr)
     END DO

  END DO
    
  WRITE (*,*) 'total size of simfire restart fields received from all workers: ', totalrecv
   
  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  CALL MPI_Reduce (MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, &
    &     0, comm, ierr)

  WRITE (*,*) 'total size of simfire restart fields sent by all workers: ', totalsend

  IF (totalrecv /= totalsend) THEN
          WRITE (*,*) 'error: simfire restart fields totalsend and totalrecv differ'
          CALL MPI_Abort (comm, 0, ierr)
  END IF

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  !=============================================================================
  ! Standard input (simfire_inp_ts)
  !=============================================================================

  ntyp = n_simfire_input
  
  ALLOCATE (blocks(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))
  
  istride  = mp * extid ! short integer 
  r1stride = mp * extr1 ! single precision
  r2stride = mp * extr2 ! double precision

  DO rank = 1, wnp
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     i1len = cnt * extid  
     r1len = cnt * extr1
     r2len = cnt * extr2

     bidx   = 0
     last2d = 0
     
     ! Population density
     bidx = bidx + 1
     CALL MPI_Get_address (SF%PopDens(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! ------------- Wrap up -------------
   
     types(last2d+1:bidx) = MPI_BYTE
     
     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE (*,*) 'invalid blaze ntyp constant, fix it!'
        CALL MPI_Abort (comm, 1, ierr)
     END IF
     
     CALL MPI_Type_create_struct (bidx, blocks, displs, types, simfire_inp_ts(rank), ierr)
     CALL MPI_Type_commit (simfire_inp_ts(rank), ierr)
     
     CALL MPI_Type_size (simfire_inp_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent (simfire_inp_ts(rank), tmplb, text, ierr)
     
     WRITE (*,*) 'restart results recv from worker, size, extent, lb: ', &
          &       rank,tsize,text,tmplb
     
     totalrecv = totalrecv + tsize
     
     ! free the partial types used for matrices
     DO i = 1, last2d
        CALL MPI_Type_free (types(i), ierr)
     END DO

  END DO

  WRITE (*,*) 'total size of simfire input fields received from all workers: ', totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  CALL MPI_Reduce (MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, &
    &     0, comm, ierr)

  WRITE (*,*) 'total size of restart fields sent by all workers: ', totalsend

  IF (totalrecv /= totalsend) THEN
          WRITE (*,*) 'error: restart fields totalsend and totalrecv differ'
          CALL MPI_Abort (comm, 0, ierr)
  END IF

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)
   
  !=============================================================================
  ! Desired output (simfire_out_ts)
  !=============================================================================

!  INTEGER, DIMENSION(:), ALLOCATABLE  :: IGBP, BIOME, REGION, NDAY
!  REAL,    DIMENSION(:), ALLOCATABLE  :: POPD, MAX_NESTEROV, CNEST, FAPAR, LAT, LON, FLI


   IF ( TRIM(SF%OUTMODE) == "full" ) THEN 
  
     ntyp = n_simfire_input
     
     ALLOCATE (blocks(ntyp))
     ALLOCATE (displs(ntyp))
     ALLOCATE (types(ntyp))
     
     istride  = mp * extid ! short integer 
     r1stride = mp * extr1 ! single precision
     r2stride = mp * extr2 ! double precision
     
     DO rank = 1, wnp
        off = wland(rank)%patch0
        cnt = wland(rank)%npatch
        
        i1len = cnt * extid  
        r1len = cnt * extr1
        r2len = cnt * extr2
        
        bidx   = 0
        last2d = 0
     
        ! BIOME
        bidx = bidx + 1
        CALL MPI_Get_address (SF%BIOME(off), displs(bidx), ierr)
        blocks(bidx) = i1len
       
        ! Pop. density
        bidx = bidx + 1
        CALL MPI_Get_address (SF%POPD(off), displs(bidx), ierr)
        blocks(bidx) = r1len
       
        ! Current MAX NESTEROV
        bidx = bidx + 1
        CALL MPI_Get_address (SF%MAX_NESTEROV(off), displs(bidx), ierr)
        blocks(bidx) = r1len
       
        ! FAPAR
        bidx = bidx + 1
        CALL MPI_Get_address (SF%FAPAR(off), displs(bidx), ierr)
        blocks(bidx) = r1len
       
        ! ------------- Wrap up -------------

        types(last2d+1:bidx) = MPI_BYTE
        
        ! MPI: sanity check
        IF (bidx /= ntyp) THEN
           WRITE (*,*) 'invalid blaze ntyp constant, fix it!'
           CALL MPI_Abort (comm, 1, ierr)
        END IF
        
        CALL MPI_Type_create_struct (bidx, blocks, displs, types, simfire_out_ts(rank), ierr)
        CALL MPI_Type_commit (simfire_out_ts(rank), ierr)
        
        CALL MPI_Type_size (simfire_out_ts(rank), tsize, ierr)
        CALL MPI_Type_get_extent (simfire_out_ts(rank), tmplb, text, ierr)
        
        WRITE (*,*) 'restart results recv from worker, size, extent, lb: ', &
             &       rank,tsize,text,tmplb
        
        totalrecv = totalrecv + tsize
        
        ! free the partial types used for matrices
        DO i = 1, last2d
           CALL MPI_Type_free (types(i), ierr)
        END DO
        
     END DO
     
     WRITE (*,*) 'total size of restart fields received from all workers: ', totalrecv
     
     ! MPI: check whether total size of received data equals total
     ! data sent by all the workers
     totalsend = 0
     CALL MPI_Reduce (MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, &
          &     0, comm, ierr)
     
     WRITE (*,*) 'total size of restart fields sent by all workers: ', totalsend
     
     IF (totalrecv /= totalsend) THEN
        WRITE (*,*) 'error: restart fields totalsend and totalrecv differ'
        CALL MPI_Abort (comm, 0, ierr)
     END IF
     
     DEALLOCATE(types)
     DEALLOCATE(displs)
     DEALLOCATE(blocks)
  ELSE
     simfire_out_ts(:) = -999 ! Only dummy setting
  END IF

END SUBROUTINE master_simfire_types

SUBROUTINE worker_simfire_types(comm, mp, SF, simfire_restart_t, simfire_inp_t, simfire_out_t)

  IMPLICIT NONE
  
  INTEGER              , INTENT(IN)  :: comm ! MPI communicator to talk to the workers
  INTEGER              , INTENT(IN)  :: mp   ! Number of gridcells
  TYPE(TYPE_SIMFIRE)   , INTENT(IN)  :: SF
  INTEGER              , INTENT(OUT) :: simfire_restart_t,simfire_inp_t,simfire_out_t

  ! MPI: temp arrays for marshalling all types into a struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types
  INTEGER :: ntyp ! number of worker's types

  INTEGER :: last2d, i

  ! MPI: block lenghts for hindexed representing all vectors
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blen

  ! MPI: block lengths and strides for hvector representing matrices
  INTEGER :: r1len, r2len, i1len, llen
  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride, istride

  INTEGER :: tsize, totalrecv, totalsend
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

  INTEGER :: rank, off, cnt
  INTEGER :: bidx, midx, vidx, ierr

  CALL MPI_Comm_rank (comm, rank, ierr)

  !=============================================================================
  ! RESTART comm simfire_restart_t) 
  !=============================================================================
  
  ntyp = n_simfire_restart
  
  ALLOCATE (blocks(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))
  
  cnt  = mp
  
  r1len = cnt * extr1
  r2len = cnt * extr2
  i1len = cnt * extid
  llen  = cnt * extl

  off  = 1
  bidx = 0

  ! ------------- 2D arrays -------------

  bidx = bidx + 1
  CALL MPI_Get_address (SF%SAV_NESTEROV(off,1), displs(bidx), ierr)
  blocks(bidx) = r1len * 12
  
  bidx = bidx + 1
  CALL MPI_Get_address (SF%SAV_FAPAR(off,1), displs(bidx), ierr)
  blocks(bidx) = r1len * FAPAR_AVG_INT

  ! ------------- 1D arrays -------------
  
  ! Integer IGBP classification
  bidx = bidx + 1
  CALL MPI_Get_address (SF%IGBP(off), displs(bidx), ierr)
  blocks(bidx) = i1len
  
  ! BIOME
  bidx = bidx + 1
  CALL MPI_Get_address (SF%BIOME(off), displs(bidx), ierr) 
  blocks(bidx) = i1len
  
  ! ------------- Wrap up -------------

  types(:) = MPI_BYTE

  CALL MPI_Type_create_struct (bidx, blocks, displs, types, simfire_restart_t, ierr)
  CALL MPI_Type_commit (simfire_restart_t, ierr)
  
  CALL MPI_Type_size (simfire_restart_t, tsize, ierr)
  CALL MPI_Type_get_extent (simfire_restart_t, tmplb, text, ierr)
  
  WRITE (*,*) 'restart struct blocks, size, extent and lb: ',rank,bidx,tsize,text,tmplb
  
  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
     
  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  !=============================================================================
  ! Standard input (simfire_inp_t)
  !=============================================================================

  ntyp = n_simfire_input
  
  ALLOCATE (blocks(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))

  cnt  = mp
  
  r1len = cnt * extr1
  r2len = cnt * extr2
  i1len = cnt * extid
  llen  = cnt * extl

  off  = 1
  bidx = 0
 
  ! Population density
  bidx = bidx + 1
  CALL MPI_Get_address (SF%PopDens(off), displs(bidx), ierr)
  blocks(bidx) = r1len

  ! ------------- Wrap up -------------
   
  types(:) = MPI_BYTE
     
  CALL MPI_Type_create_struct (bidx, blocks, displs, types, simfire_inp_t, ierr)
  CALL MPI_Type_commit (simfire_inp_t, ierr)
  
  CALL MPI_Type_size (simfire_inp_t, tsize, ierr)
  CALL MPI_Type_get_extent (simfire_inp_t, tmplb, text, ierr)
  
  WRITE (*,*) 'restart struct blocks, size, extent and lb: ',rank,bidx,tsize,text,tmplb
  
  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
     
  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  !=============================================================================
  ! Desired output (simfire_out_t)
  !=============================================================================

  IF ( TRIM(SF%OUTMODE) == "full" ) THEN 
     
     ntyp = n_simfire_input
     
     ALLOCATE (blocks(ntyp))
     ALLOCATE (displs(ntyp))
     ALLOCATE (types(ntyp))
     
     istride  = mp * extid ! short integer 
     r1stride = mp * extr1 ! single precision
     r2stride = mp * extr2 ! double precision

     off  = 1
     bidx = 0

     ! BIOME
     bidx = bidx + 1
     CALL MPI_Get_address (SF%BIOME(off), displs(bidx), ierr)
     blocks(bidx) = i1len
     
     ! Pop. density
     bidx = bidx + 1
     CALL MPI_Get_address (SF%POPD(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     
     ! Current MAX NESTEROV
     bidx = bidx + 1
     CALL MPI_Get_address (SF%MAX_NESTEROV(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     
     ! FAPAR
     bidx = bidx + 1
     CALL MPI_Get_address (SF%FAPAR(off), displs(bidx), ierr)
     blocks(bidx) = r1len
       
     ! ------------- Wrap up -------------

     types(:) = MPI_BYTE

     CALL MPI_Type_create_struct (bidx, blocks, displs, types, simfire_out_t, ierr)
     CALL MPI_Type_commit (simfire_out_t, ierr)
     
     CALL MPI_Type_size (simfire_out_t, tsize, ierr)
     CALL MPI_Type_get_extent (simfire_out_t, tmplb, text, ierr)
  
     WRITE (*,*) 'restart struct blocks, size, extent and lb: ',rank,bidx,tsize,text,tmplb
     
     ! MPI: check whether total size of received data equals total
     ! data sent by all the workers
     CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
     
     DEALLOCATE(types)
     DEALLOCATE(displs)
     DEALLOCATE(blocks)
  ELSE
     simfire_out_t = -999 ! Only dummy setting
  END IF
 
END SUBROUTINE worker_simfire_types
  
END MODULE BLAZE_MPI
