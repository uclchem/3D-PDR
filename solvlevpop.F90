SUBROUTINE solvlevpop(NLEV,TRANSITION,density,SOLUTION)!,coolant)

  use definitions
  use healpix_types
!  use maincode_module, only : p,iteration,pdr,gridpoint,nrays

  implicit none  
  integer(kind=i4b), intent(in) :: NLEV
!  integer(kind=i4b), intent(in) :: coolant
  real(kind=dp), intent(in) :: density
  real(kind=dp), intent(in) :: transition(1:NLEV,1:NLEV)
  real(kind=dp), intent(out) :: SOLUTION(1:NLEV)
  integer(kind=i4b) :: i,j!,ii
  real(kind=dp) :: out1
  real(kind=dp) :: A(1:NLEV,1:NLEV)
  logical::call_writes
  !real(kind=dp) :: temp_a(1:nlev,1:nlev),temp_solution(1:nlev)

         A=0.0D0
!        Fill the matrix
         DO I=1,NLEV
            OUT1=0.0D0
            DO J=1,NLEV
               OUT1=OUT1+TRANSITION(I,J)
               A(I,J)=TRANSITION(J,I)
            ENDDO
            A(I,I)=-OUT1
         ENDDO
!        Initialize the solution array before calling the solver routine
         DO I=1,NLEV
            SOLUTION(I)=0.0D0
            A(NLEV,I)=1.0D-8 !non-zero starting parameter to avoid division by zero.
         ENDDO

         SOLUTION(NLEV)=DENSITY*1.0D-8

         !CALL GAUSS_JORDAN(A,NLEV,NLEV,SOLUTION,coolant,call_writes)
         CALL GAUSS_JORDAN(A,NLEV,NLEV,SOLUTION,call_writes)

!        Replace negative level populations due to numerical noise around 0
         DO I=1,NLEV
            if (solution(i).lt.0.0D0) solution(i)=0.0D0!1.0D-99!then !stop 'found negative solution!'
!              write(6,*) '';write(6,*) 'found negative solution in p=',p;write(6,*) 'coolant=',coolant;write(6,*)''
!              call gauss_jordan_writes(temp_a,nlev,nlev,temp_solution,coolant,i)
!              stop
!            endif
         ENDDO
         
     return
end subroutine

!C-----------------------------------------------------------------------
!C Standard Gauss-Jordon linear equation solver from Numerical Recipes
!C A(N,N) is an input matrix stored in an array of dimensions NPxNP
!C B(N,M) is an input matrix containing the M right-hand side vectors
!C stored in an array of dimensions NPxMPP
!C
!C On output, A(N,N) is replaced by its matrix inverse and B(N,M)
!C is replaced by the corresponding set of solution vectors
!C
!C Note: set NMAX to the maximum possible dimension (NLEVEL)
!C-----------------------------------------------------------------------
!   SUBROUTINE GAUSS_JORDAN(A,N,NP,B,M,MPP,coolant)
!   SUBROUTINE GAUSS_JORDAN(A,N,NP,B,coolant,call_writes)
   SUBROUTINE GAUSS_JORDAN(A,N,NP,B,call_writes)

      use definitions
      use healpix_types
      use maincode_module, only : gastemperature,p,iteration

      IMPLICIT NONE
!      integer(kind=i4b), intent(in) :: coolant
      INTEGER(kind=i4b):: I,J,K,L,LL,IROW,ICOL
      INTEGER(kind=i4b), intent(in):: N,NP!,M,MPP
      integer(kind=i4b), PARAMETER :: NMAX=100
      INTEGER(kind=i4b):: IPIV(1:NMAX),INDXR(1:NMAX),INDXC(1:NMAX)
      real(kind=dp), intent(inout) :: A(1:NP,1:NP)
      real(kind=dp), intent(inout) :: B(1:NP)!,1:MPP)
      real(kind=dp) :: BIG,DUM,PIVINV
      logical,intent(out)::call_writes

      ICOL=0
      IROW=0
      IPIV=0
      DO I=1,N
         BIG=0.0D0
         DO J=1,N
            IF(IPIV(J).NE.1) THEN
               DO K=1,N
                  IF(IPIV(K).EQ.0) THEN
                     IF(ABS(A(J,K)).GE.BIG) THEN
                        BIG=ABS(A(J,K))
                        IROW=J
                        ICOL=K
                     ENDIF
                  ELSE IF(IPIV(K).GT.1) THEN
                     PRINT *,'ERROR! Singular matrix in GAUSS_JORDAN'
call_writes=.true.
return
!                     write(6,*) 'Crashed in first loop'
!                     write(6,*) 'grid point = ',p, ' coolant = ',coolant
!                     write(6,*) 'gastemperature = ',gastemperature(p)
!                     STOP
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         IPIV(ICOL)=IPIV(ICOL)+1
         IF(IROW.NE.ICOL) THEN
            DO L=1,N
               DUM=A(IROW,L)
               A(IROW,L)=A(ICOL,L)
               A(ICOL,L)=DUM
            ENDDO
!            DO L=1,M
!               DUM=B(IROW,L)
!               B(IROW,L)=B(ICOL,L)
!               B(ICOL,L)=DUM
!            ENDDO
!================================
            DUM=B(IROW)
            B(IROW)=B(ICOL)
            B(ICOL)=DUM
!================================
         ENDIF
         INDXR(I)=IROW
         INDXC(I)=ICOL
         IF(A(ICOL,ICOL).EQ.0.0D0) THEN
            PRINT *,'ERROR! Singular matrix found by GAUSS_JORDAN'
call_writes=.true.
return
!            write(6,*) 'Crashed in second loop'
!            write(6,*) 'grid point = ',p, ' coolant = ',coolant
!            write(6,*) 'gastemperature = ',gastemperature(p)
!            STOP
         ENDIF
         PIVINV=1.0D0/A(ICOL,ICOL)
         A(ICOL,ICOL)=1.0D0
         DO L=1,N
            A(ICOL,L)=A(ICOL,L)*PIVINV
         ENDDO
!         DO L=1,M
!            B(ICOL,L)=B(ICOL,L)*PIVINV
!         ENDDO
!=======================================
         B(ICOL)=B(ICOL)*PIVINV
!=======================================
         DO LL=1,N
            IF(LL.NE.ICOL) THEN
               DUM=A(LL,ICOL)
               A(LL,ICOL)=0.0D0
               DO L=1,N
                  A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
               ENDDO
!               DO L=1,M
!                  B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
!               ENDDO
!=============================================
               B(LL)=B(LL)-B(ICOL)*DUM
!=============================================
            ENDIF
         ENDDO
      ENDDO
      DO L=N,1,-1
         IF(INDXR(L).NE.INDXC(L)) THEN
            DO K=1,N
               DUM=A(K,INDXR(L))
               A(K,INDXR(L))=A(K,INDXC(L))
               A(K,INDXC(L))=DUM
            ENDDO
         ENDIF
      ENDDO
      RETURN
      END subroutine


!C-----------------------------------------------------------------------
!C Standard Gauss-Jordon linear equation solver from Numerical Recipes
!C A(N,N) is an input matrix stored in an array of dimensions NPxNP
!C B(N,M) is an input matrix containing the M right-hand side vectors
!C stored in an array of dimensions NPxMPP
!C
!C On output, A(N,N) is replaced by its matrix inverse and B(N,M)
!C is replaced by the corresponding set of solution vectors
!C
!C Note: set NMAX to the maximum possible dimension (NLEVEL)
!C-----------------------------------------------------------------------
!   SUBROUTINE GAUSS_JORDAN(A,N,NP,B,M,MPP,coolant)
!   SUBROUTINE GAUSS_JORDAN_writes(A,N,NP,B,coolant,ill)
   SUBROUTINE GAUSS_JORDAN_writes(A,N,NP,B,ill)

      use definitions
      use healpix_types
      use maincode_module, only : gastemperature,p,iteration

      IMPLICIT NONE
      integer(kind=i4b), intent(in) :: ill!,coolant
      INTEGER(kind=i4b):: I,J,K,L,LL,IROW,ICOL
      INTEGER(kind=i4b), intent(in):: N,NP!,M,MPP
      integer(kind=i4b), PARAMETER :: NMAX=100
      INTEGER(kind=i4b):: IPIV(1:NMAX),INDXR(1:NMAX),INDXC(1:NMAX)
      real(kind=dp), intent(inout) :: A(1:NP,1:NP)
      real(kind=dp), intent(inout) :: B(1:NP)!,1:MPP)
      real(kind=dp) :: BIG,DUM,PIVINV

write(6,*) 'b'
do i=1,np
write(6,*) b(i),i
enddo

write(6,*) 'a'
do i=1,np
  do j=1,np
    write(6,*) a(i,j)
  enddo
enddo


      ICOL=0
      IROW=0
      IPIV=0
      DO I=1,N
         BIG=0.0D0
         DO J=1,N
            IF(IPIV(J).NE.1) THEN
               DO K=1,N
                  IF(IPIV(K).EQ.0) THEN
                     IF(ABS(A(J,K)).GE.BIG) THEN
                        BIG=ABS(A(J,K))
                        IROW=J
                        ICOL=K
                     ENDIF !ABS(A
                  ELSE IF(IPIV(K).GT.1) THEN
                     PRINT *,'ERROR! Singular matrix in GAUSS_JORDAN'
                     write(6,*) 'Crashed in first loop'
!                     write(6,*) 'grid point = ',p, ' coolant = ',coolant
                     write(6,*) 'gastemperature = ',gastemperature(p)
                     STOP
                  ENDIF !IPIV(K).EQ.0
               ENDDO !K=1,N
            ENDIF !IPIV(J).NE.1
         ENDDO !J=1,N
         IPIV(ICOL)=IPIV(ICOL)+1
         IF(IROW.NE.ICOL) THEN
            DO L=1,N
               DUM=A(IROW,L)
               A(IROW,L)=A(ICOL,L)
               A(ICOL,L)=DUM
            ENDDO !L=1,N
!            DO L=1,M
!               DUM=B(IROW,L)
!               B(IROW,L)=B(ICOL,L)
!               B(ICOL,L)=DUM
!            ENDDO
!================================
            DUM=B(IROW)
if (i.eq.ill) write(6,*) 'dum=',dum,'A'
!write(6,*) 'DUM=',DUM
            B(IROW)=B(ICOL)
if (i.eq.ill) write(6,*) 'b(',irow,')=',b(irow),'B'
!write(6,*) 'irow=',irow
!write(6,*) 'B(irow)=',b(irow)
            B(ICOL)=DUM
if (i.eq.ill) write(6,*) 'b(',icol,')=',b(icol),'C'
!write(6,*) 'icol=',icol
!write(6,*) 'b(icol)=',b(icol)
!================================
         ENDIF !IROW.NE.ICOL
         INDXR(I)=IROW
         INDXC(I)=ICOL
         IF(A(ICOL,ICOL).EQ.0.0D0) THEN
            PRINT *,'ERROR! Singular matrix found by GAUSS_JORDAN'
            write(6,*) 'Crashed in second loop'
!            write(6,*) 'grid point = ',p, ' coolant = ',coolant
            write(6,*) 'gastemperature = ',gastemperature(p)
            STOP
         ENDIF
         PIVINV=1.0D0/A(ICOL,ICOL)
         A(ICOL,ICOL)=1.0D0
         DO L=1,N
            A(ICOL,L)=A(ICOL,L)*PIVINV
         ENDDO
!         DO L=1,M
!            B(ICOL,L)=B(ICOL,L)*PIVINV
!         ENDDO
!=======================================
if (i.eq.ill) write(6,*) b(icol),pivinv,'D'
         B(ICOL)=B(ICOL)*PIVINV
if (i.eq.ill) write(6,*) b(icol),'E'
!write(6,*) 'pivinv=',pivinv
!write(6,*) 'b(icol)=',b(icol)
!=======================================
         DO LL=1,N
            IF(LL.NE.ICOL) THEN
               DUM=A(LL,ICOL)
               A(LL,ICOL)=0.0D0
               DO L=1,N
                  A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
               ENDDO
!               DO L=1,M
!                  B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
!               ENDDO
!=============================================
if (i.eq.ill) then
write(6,*) 'll=',ll,'F'
write(6,*) 'b(ll)=',b(ll),'G'
write(6,*) 'b(icol)=',b(icol),'H'
write(6,*) 'dum=',dum,'I'
endif
               B(LL)=B(LL)-B(ICOL)*DUM
if (i.eq.ill) write(6,*) 'b(ll) after=',b(ll),'J'
!=============================================
            ENDIF !LL.NE.ICOL
         ENDDO !LL=1,N
      ENDDO ! I=1,N
      DO L=N,1,-1
         IF(INDXR(L).NE.INDXC(L)) THEN
            DO K=1,N
               DUM=A(K,INDXR(L))
               A(K,INDXR(L))=A(K,INDXC(L))
               A(K,INDXC(L))=DUM
            ENDDO !K=1,N
         ENDIF !INDXR(L).NE.INDXC(L)
      ENDDO !L=N,1,-1
do i=1,n
write(6,*) b(i),i
enddo
      RETURN
      END subroutine
