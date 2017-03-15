!=======================================================================
!
!     Copied from Numerical Recipes
!
! Given a tabulated function YA (of size MxN) and tabulated independent
! variables X1A (M values) and X2A (N values), this routine constructs
! one-dimensional natural cubic splines of the rows of YA and returns
! the second derivatives in the array Y2A.
!
!-----------------------------------------------------------------------
      SUBROUTINE SPLIE2(X1A,X2A,YA,M,N,Y2A)

      USE DEFINITIONS
      USE HEALPIX_TYPES

      IMPLICIT NONE

      INTEGER(KIND=I4B), INTENT(IN) :: M,N
      REAL(KIND=DP), INTENT(IN)     :: X1A(1:M),X2A(1:N),YA(1:M,1:N)
      REAL(KIND=DP), INTENT(OUT)    :: Y2A(1:M,1:N)

      INTEGER(KIND=I4B) :: I,J
      REAL(KIND=DP) :: YTMP(1:N),Y2TMP(1:N)

      DO I=1,M
         DO J=1,N
            YTMP(J)=YA(I,J)
         ENDDO
!        Values of 1.0D30 indicate a natural spline
         CALL SPLINE(X2A,YTMP,N,1.0D30,1.0D30,Y2TMP)
         DO J=1,N
            Y2A(I,J)=Y2TMP(J)
         ENDDO
      ENDDO

      RETURN
      END
!=======================================================================

!=======================================================================
!
!     Calculate the cubic spline for a set of points (X,Y)
!     (c.f. Numerical Recipes, Chapter 3.3: Spline Routine)
!
!     Given the arrays X and Y (size N) containing a tabulated
!     function, i.e., Y(I)=f(X(I)), with X(1) < X(2) < ... < X(N),
!     and given values YP1 and YPN for the first derivative of the
!     interpolating function at points 1 and N, respectively, this
!     routine returns an array Y2 of length N, which contains the
!     second derivatives of the interpolating function at the
!     tabulated points X(I). If YP1 and/or YPN are equal to 1.0E+30
!     or larger, the routine is signalled to set the corresponding
!     boundary condition for a natural spline, with zero second
!     derivative at that boundary.
!
!     I/O parameters:
!     Input   X   = vector for independent variable; dimension X(1:N)
!     Input   Y   = vector for x-dependent variable; dimension Y(1:N)
!     Input   N   = dimension of vectors containing tabulated function
!     Input   YP1 = 1st derivative of the function at point 1
!     Input   YPN = 1st derivative of the function at point N
!     Output  Y2  = 2nd derivative of the function; dimension Y2(1:N)
!
!-----------------------------------------------------------------------
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)

      USE DEFINITIONS
      USE HEALPIX_TYPES

      IMPLICIT NONE

      INTEGER(KIND=I4B), INTENT(IN) :: N
      REAL(KIND=DP), INTENT(IN)     :: X(1:N),Y(1:N)
      REAL(KIND=DP), INTENT(IN)     :: YP1,YPN
      REAL(KIND=DP), INTENT(OUT)    :: Y2(1:N)

      INTEGER(KIND=I4B) :: I
      REAL(KIND=DP) :: P,QN,SIG,U(1:N),UN

      IF(YP1.GE.1.0D30) THEN
!        The lower boundary condition is either set to be "natural"
         Y2(1)=0.0D0
         U(1)=0.0D0
      ELSE
!        or to have a specified first derivative
         Y2(1)=-0.5D0
         U(1)=(3.0D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF

!     This is the decomposition loop of the tridiagonal algorithm
!     Y2 and U are used for temporary storage of the decomposed factors
      DO I=2,N-1
         SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
         P=SIG*Y2(I-1)+2.0D0
         Y2(I)=(SIG-1.0D0)/P
         U(I)=(6.0D0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))&
     &              /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      ENDDO

      IF(YPN.GE.1.0D30) THEN
!        The upper boundary condition is either set to be "natural"
         QN=0.0D0
         UN=0.0D0
      ELSE
!        or to have a specified first derivative
         QN=0.5D0
         UN=(3.0D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF

      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.0D0)
!     This is the back-substitution loop of the tridiagonal algorithm
      DO I=N-1,1,-1
         Y2(I)=Y2(I)*Y2(I+1)+U(I)
      ENDDO

      RETURN
      END
!=======================================================================

!=======================================================================
!
!     Given X1A, X2A, YA, M, N (as described in SPLIE2) and Y2A (as
!     produced by that routine), and given a desired interpolating
!     point (X1,X2), this routine returns an interpolated function
!     value Y by performing a bicubic spline interpolation.
!
!-----------------------------------------------------------------------
      SUBROUTINE SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)

      USE DEFINITIONS
      USE HEALPIX_TYPES

      IMPLICIT NONE

      INTEGER(KIND=I4B), INTENT(IN) :: M,N
      REAL(KIND=DP), INTENT(IN)     :: X1A(1:M),X2A(1:N),YA(1:M,1:N),Y2A(1:M,1:N)
      REAL(KIND=DP), INTENT(IN)     :: X1,X2
      REAL(KIND=DP), INTENT(OUT)    :: Y

      INTEGER(KIND=I4B) :: I,J
      REAL(KIND=DP) :: YTMP(1:N),Y2TMP(1:N),YYTMP(1:M),YY2TMP(1:M)

!     Perform M evaluations of the row splines constructed by
!     SPLIE2 using the one-dimensional spline evaluator SPLINT
      DO I=1,M
         DO J=1,N
            YTMP(J)=YA(I,J)
            Y2TMP(J)=Y2A(I,J)
         ENDDO
         CALL SPLINT(X2A,YTMP,Y2TMP,N,X2,YYTMP(I))
      ENDDO

!     Construct the one-dimensional column spline and evaluate it
!     Values of 1.0D30 indicate a natural spline
      CALL SPLINE(X1A,YYTMP,M,1.0D30,1.0D30,YY2TMP)
      CALL SPLINT(X1A,YYTMP,YY2TMP,M,X1,Y)

      RETURN
      END
!=======================================================================

!=======================================================================
!
!     Perform a cubic spline interpolation evaluated at the point X
!     (c.f. Numerical Recipes, Chapter 3.3: Splint Routine, 
!           Numerical Recipes, Chapter 3.4: Hunt Routine)
!
!     Given the arrays XA and YA (size N) containing a tabulated
!     function, i.e., YA(I) = f(XA(I)), with the XA(I)'s in order,
!     and given the array Y2A produced by the SPLINE routine, this
!     routine returns a cubic spline interpolated value Y.
!
!     I/O parameters:
!     Input   XA  = vector for independent variable; dimension XA(1:N)
!     Input   YA  = vector for x-dependent variable; dimension YA(1:N)
!     Input   Y2A = 2nd derivative of the function; dimension Y2A(1:N)
!     Input   N   = dimension of input vectors
!     Input   X   = x-value at which Y is to be interpolated
!     Output  Y   = result of interpolation
!
!-----------------------------------------------------------------------
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)

      USE DEFINITIONS
      USE HEALPIX_TYPES

      IMPLICIT NONE

      INTEGER(KIND=I4B), INTENT(IN) :: N
      REAL(KIND=DP), INTENT(IN)     :: XA(1:N),YA(1:N),Y2A(1:N)
      REAL(KIND=DP), INTENT(IN)     :: X
      REAL(KIND=DP), INTENT(OUT)    :: Y

      LOGICAL :: ASCND
      INTEGER(KIND=I4B) :: JLO,JHI,JMID,INC
      REAL(KIND=DP) :: A,B

      JLO=0
      JHI=0

!     ASCND is TRUE if the table values are in ascending order, FALSE otherwise
      ASCND=XA(N).GT.XA(1)

!     Find the interval XA(JLO) <= X <= XA(JLO+1) = XA(JHI)
      IF(JLO.LE.0 .OR. JLO.GT.N) THEN
!        Input guess not useful, go immediately to bisection
         JLO=0
         JHI=N+1
         GOTO 300
      ENDIF

!     Set the hunting increment
      INC=1

      IF(X.GE.XA(JLO) .EQV. ASCND) THEN
!        Hunt up:
 100     JHI=JLO+INC
         IF(JHI.GT.N) THEN
!           Done hunting, since off the end of the table
            JHI=N+1
         ELSE IF(X.GE.XA(JHI) .EQV. ASCND) THEN
!           Not done hunting...
            JLO=JHI
!           ...so double the increment...
            INC=INC+INC
!           ...and try again
            GOTO 100
         ENDIF
!     Done hunting, value bracketed
      ELSE
         JHI=JLO
!        Hunt down:
 200     JLO=JHI-INC
         IF(JLO.LT.1) THEN
!           Done hunting, since off the end of the table
            JLO=0
         ELSE IF(X.LT.XA(JLO) .EQV. ASCND) THEN
!           Not done hunting...
            JHI=JLO
!           ...so double the increment...
            INC=INC+INC
!           ...and try again
            GOTO 200
         ENDIF
!     Done hunting, value bracketed
      ENDIF

 300  IF((JHI-JLO).NE.1) THEN
!        Hunt is done, so begin the final bisection phase
         JMID=(JHI+JLO)/2
         IF(X.GT.XA(JMID) .EQV. ASCND) THEN
            JLO=JMID
         ELSE
            JHI=JMID
         ENDIF
         GOTO 300
      ENDIF

      IF(JLO.EQ.0) THEN
         JLO=1
         JHI=2
      ENDIF
      IF(JLO.EQ.N) THEN
         JLO=N-1
         JHI=N
      ENDIF

!     JLO and JHI now bracket the input value X
!     The cubic spline polynomial is now evaluated
      A=(XA(JHI)-X)/(XA(JHI)-XA(JLO))
      B=(X-XA(JLO))/(XA(JHI)-XA(JLO))
      Y=A*YA(JLO)+B*YA(JHI)+((A**3-A)*Y2A(JLO)+(B**3-B)*Y2A(JHI))&
     &  *((XA(JHI)-XA(JLO))**2)/6.0D0

      RETURN
      END
!=======================================================================
