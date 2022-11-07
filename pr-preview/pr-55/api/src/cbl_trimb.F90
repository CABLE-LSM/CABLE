MODULE trimb_mod

USE cbl_ssnow_data_mod, ONLY : r_2

PUBLIC  trimb

CONTAINS
! this routine solves the system:
!	   a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)=rhs(k)    for k=2,kmax-1
!	   with  b(k)*u(k)+c(k)*u(k+1)=rhs(k)	       for k=1
!	   and	 a(k)*u(k-1)+b(k)*u(k)=rhs(k)	       for k=kmax
!
!	 the Thomas algorithm is used for solving sets of linear equation
!	 rhs initially contains rhs; leaves with answer (jlm)
!	 n.b. this one does not assume b = 1-a-c

SUBROUTINE trimb (a, b, c, rhs, kmax)

IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: kmax ! no. of discrete layers

  REAL(r_2), DIMENSION(:,:), INTENT(IN) ::                                    &
       a,    & ! coef "A" in finite diff eq
       b,    & ! coef "B" in finite diff eq
       c       ! coef "C" in finite diff eq

  REAL(r_2), DIMENSION(:,:), INTENT(INOUT)  :: rhs ! right hand side of eq

  REAL(r_2), DIMENSION(SIZE(a,1),SIZE(a,2)) ::                                &
       e, temp, g

  INTEGER :: k   ! do lloop counter

  e(:,1) = c(:,1) / b(:,1)
  DO k = 2, kmax - 1
     temp(:,k) = 1. / ( b(:,k) - a(:,k) * e(:,k-1) )
     e(:,k) = c(:,k) * temp(:,k)
  END DO

  g(:,1) = rhs(:,1) / b(:,1)
  DO k = 2, kmax - 1
     g(:,k) = ( rhs(:,k) - a(:,k) * g(:,k-1) ) * temp(:,k)
  END DO

  ! do back substitution to give answer now
  rhs(:,kmax) = ( rhs(:,kmax) - a(:,kmax) * g(:,kmax-1) )                     &
       / ( b(:,kmax) - a(:,kmax) * e(:,kmax-1) )

  DO k = kmax - 1, 1, - 1
     rhs(:,k) = g(:,k) - e(:,k) * rhs(:,k + 1)
  END DO

END SUBROUTINE trimb

END MODULE trimb_mod
