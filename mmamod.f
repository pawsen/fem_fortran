C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE ASYMPT(ITER,M,N,XVAL,XMIN,XMAX,XOLD1,XOLD2,
     1                  XLOW,XUPP,ALFA,BETA)
C
C       Version "November 2000".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     ASYMPT calculates the asymptotes XLOW and XUPP,
C     and the bounds ALFA and BETA, for the current subproblem.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XVAL(*),XMIN(*),XMAX(*),XOLD1(*),XOLD2(*),
     1          XLOW(*),XUPP(*),ALFA(*),BETA(*)
C
      ALBEFA=0.1
      GHINIT=0.5
!     GHMOVE=0.7

      ghmove = 0.65  ! teb
      ghmove2 = 1.07 ! teb

      IF(ITER.GE.3) GOTO 350
C
C***  Here ITER = 1 or 2 .
      DO 200 J=1,N
      XMMJ=XMAX(J)-XMIN(J)
      IF(XMMJ.LT.0.00001) XMMJ=0.00001
      XLOW(J)=XVAL(J)-GHINIT*XMMJ
      XUPP(J)=XVAL(J)+GHINIT*XMMJ
  200 CONTINUE
      GOTO 500
C
C***  Here ITER is greater than 2.
  350 CONTINUE
C
      DO 400 J=1,N
      XTEST=(XVAL(J)-XOLD1(J))*(XOLD1(J)-XOLD2(J))
      FAK=1.0
      IF(XTEST.LT.0.) FAK=GHMOVE
!     IF(XTEST.GT.0.) FAK=2.0/(GHMOVE+1.0)

      if (xtest .gt. 0.) fak = ghmove2  ! teb

      XLOW(J)=XVAL(J)-FAK*(XOLD1(J)-XLOW(J))
      XUPP(J)=XVAL(J)+FAK*(XUPP(J)-XOLD1(J))
      XMMJ=XMAX(J)-XMIN(J)
      IF(XMMJ.LT.0.00001) XMMJ=0.00001
      GMINJ = XVAL(J)-100.0*XMMJ
      GMAXJ = XVAL(J)-0.00001*XMMJ
      HMINJ = XVAL(J)+0.00001*XMMJ
      HMAXJ = XVAL(J)+100.0*XMMJ
      IF(XLOW(J).LT.GMINJ) XLOW(J)=GMINJ
      IF(XLOW(J).GT.GMAXJ) XLOW(J)=GMAXJ
      IF(XUPP(J).LT.HMINJ) XUPP(J)=HMINJ
      IF(XUPP(J).GT.HMAXJ) XUPP(J)=HMAXJ
  400 CONTINUE
C
  500 CONTINUE
C
      DO 600 J=1,N ! (skal vaere n og ikke nvar)
      XMIJ=XMIN(J)-0.000001
      XMAJ=XMAX(J)+0.000001
      IF(XVAL(J).GE.XMIJ) GOTO 550
      XLOW(J)=XVAL(J)-(XMAJ-XVAL(J))/0.9
      XUPP(J)=XVAL(J)+(XMAJ-XVAL(J))/0.9
      GOTO 600
  550 CONTINUE
      IF(XVAL(J).LE.XMAJ) GOTO 600
      XLOW(J)=XVAL(J)-(XVAL(J)-XMIJ)/0.9
      XUPP(J)=XVAL(J)+(XVAL(J)-XMIJ)/0.9
  600 CONTINUE
C
      DO 700 J=1,N
      ALFA(J)=XLOW(J)+ALBEFA*(XVAL(J)-XLOW(J))
      BETA(J)=XUPP(J)-ALBEFA*(XUPP(J)-XVAL(J))
      IF(ALFA(J).LT.XMIN(J)) ALFA(J)=XMIN(J)
      IF(BETA(J).GT.XMAX(J)) BETA(J)=XMAX(J)
  700 CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE GENSUB(M,N,XVAL,XMIN,XMAX,F0VAL,DF0DX,FMAX,FVAL,
     1                  DFDX,P,Q,B,P0,Q0,R0,XLOW,XUPP)
C
C       Version "November 2000".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     GENSUB calculates P( ),Q( ),B( ),P0( ),Q0( ) and R0
C     for the current subproblem.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XVAL(*),XMIN(*),XMAX(*),DF0DX(*),FMAX(*),FVAL(*),
     1          DFDX(*),P(*),Q(*),B(*),P0(*),Q0(*),XLOW(*),XUPP(*)
C
      RAA0=0.000001
      R0=F0VAL
      DO 20 I=1,M
      B(I)=FMAX(I)-FVAL(I)
   20 CONTINUE
C
      DO 50 J=1,N
      MJ=M*(J-1)
      UJLJ=XUPP(J)-XLOW(J)
      UJXJ=XUPP(J)-XVAL(J)
      XJLJ=XVAL(J)-XLOW(J)
      UJXJ2=UJXJ*UJXJ
      XJLJ2=XJLJ*XJLJ
      P0J=0.5*RAA0/UJLJ
      Q0J=0.5*RAA0/UJLJ
      IF(DF0DX(J).GT.0.) P0J=P0J+1.001*DF0DX(J)
      IF(DF0DX(J).GT.0.) Q0J=Q0J+0.001*DF0DX(J)
      IF(DF0DX(J).LT.0.) Q0J=Q0J-1.001*DF0DX(J)
      IF(DF0DX(J).LT.0.) P0J=P0J-0.001*DF0DX(J)
      P0J=P0J*UJXJ2
      Q0J=Q0J*XJLJ2
      P0(J)=P0J
      Q0(J)=Q0J
      R0=R0-P0J/UJXJ-Q0J/XJLJ
C
      DO 40 I=1,M
      IJ=MJ+I
      PIJ=0.5*RAA0/UJLJ
      QIJ=0.5*RAA0/UJLJ
      DFIJ=DFDX(IJ)
      IF(DFIJ.GT.0.) PIJ=PIJ+1.001*DFIJ
      IF(DFIJ.GT.0.) QIJ=QIJ+0.001*DFIJ
      IF(DFIJ.LT.0.) QIJ=QIJ-1.001*DFIJ
      IF(DFIJ.LT.0.) PIJ=PIJ-0.001*DFIJ
      PIJ=PIJ*UJXJ2
      QIJ=QIJ*XJLJ2
      P(IJ)=PIJ
      Q(IJ)=QIJ
      B(I)=B(I)+PIJ/UJXJ+QIJ/XJLJ
C
   40 CONTINUE
   50 CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE MAXIM(M,N,GEPS,IYFREE,GRADF,DSRCH,HESSF,X,Y,Z,ULAM,
     1                 UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     MAXIM solves the dual MMA subproblem.
C     The dual variables are ULAM(I), I=1,..,M.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GRADF(*),DSRCH(*),HESSF(*),X(*),Y(*),ULAM(*),
     1          UU(*),XLOW(*),XUPP(*),ALFA(*),BETA(*),
     2          A(*),B(*),C(*),P(*),Q(*),P0(*),Q0(*)
      INTEGER IYFREE(*)
C
      ITR=0
      M3=3*M+10
C
      DO 10 I=1,M
      ULAM(I)=0.
      IYFREE(I)=1
 10   CONTINUE
C
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA,
     1           A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
C
      GMX=0.
      DO 20 I=1,M
      IYFREE(I)=0
      IF(GRADF(I).GE.GEPS) IYFREE(I)=1
      IF(GRADF(I).GE.GMX) GMX=GRADF(I)
 20   CONTINUE
C
      IF(GMX.LE.GEPS) GOTO 100
C
 30   CONTINUE
      ITR=ITR+1
      IF(ITR.GT.M3) GOTO 100
C
      CALL SUBSPA(ITR,M,N,GEPS,F,IYFREE,GRADF,DSRCH,HESSF,
     1            X,Y,ULAM,UU,XLOW,XUPP,ALFA,BETA,A,B,C,
     2            P,Q,P0,Q0,IHITY)
C
      IF(IHITY.EQ.0) GOTO 40
      IYFREE(IHITY)=0
      ULAM(IHITY)=0.
      GOTO 30
C
 40   CONTINUE
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA,
     1           A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
C
      GMX=0.
      IGMX=0
      DO 50 I=1,M
      IF(IYFREE(I).EQ.1) GOTO 50
      IF(GRADF(I).LE.GMX) GOTO 50
      GMX=GRADF(I)
      IGMX=I
 50   CONTINUE
C
      IF(GMX.LE.GEPS) GOTO 100
      IYFREE(IGMX)=1
      GOTO 30
C
 100  CONTINUE
      CALL XYZLAM(M,N,X,Y,Z,ULAM,XLOW,XUPP,ALFA,BETA,
     1            A,B,C,P,Q,P0,Q0,IYFREE)
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE SUBSPA(ITR,M,N,GEPS,F,IYFREE,GRADF,DSRCH,HESSF,
     1                  X,Y,ULAM,UU,XLOW,XUPP,ALFA,BETA,
     2                  A,B,C,P,Q,P0,Q0,IHITY)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C    SUBSPA maximizes the dual objective function on the subspace
C    defined by ULAM(I) = 0 for every I such that IYFREE(I) = 0.
C    The first iteration is a steepest ascent step,
C    the second and third iterations are conjugate gradient steps,
C    and after that a Newton method is used.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GRADF(*),DSRCH(*),HESSF(*),X(*),Y(*),
     1          ULAM(*),UU(*),XLOW(*),XUPP(*),ALFA(*),BETA(*),
     2          A(*),B(*),C(*),P(*),Q(*),P0(*),Q0(*)
      INTEGER IYFREE(*)
C
      IHITY=0
      ITESUB=0
      NYDIM=0
      DSRTOL=-0.0000001*GEPS
      GNORM2=0.
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA,
     1           A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
C
      DO 10 I=1,M
      DSRCH(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 10
      NYDIM=NYDIM+1
      GNORM2=GNORM2+GRADF(I)**2
      DSRCH(I)=GRADF(I)
 10   CONTINUE
C
      IF(NYDIM.EQ.0) GOTO 180
      ITEMAX=50+5*NYDIM
C
 15   ITESUB=ITESUB+1
      TMAX=1.0D20
      IHITY=0
C
      GTD=0.
      DO 20 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 20
      GTD=GTD+GRADF(I)*DSRCH(I)
      IF(DSRCH(I).GT.DSRTOL) GOTO 20
      T=ULAM(I)/(-DSRCH(I))
      IF(T.GE.TMAX) GOTO 20
      TMAX=T
      IHITY=I
 20   CONTINUE
      IF(GTD.LE.0.) GOTO 180
C
      CALL LINESE(M,N,TMAX,TOPT,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,
     2            IYFREE,IHITY,IHITMX,ITESUB)
C
      DO 30 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 30
      ULAM(I)=ULAM(I)+TOPT*DSRCH(I)
 30   CONTINUE
C
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA,
     1           A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
C
      IF(IHITMX.EQ.1) GOTO 180
      IHITY=0
      GNOLD2=GNORM2
      GNORM2=0.
      IOPT=1
C
      DO 40 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 40
      IF(DABS(GRADF(I)).GT.GEPS) IOPT=0
      GNORM2=GNORM2+GRADF(I)**2
 40   CONTINUE
C
      IF(IOPT.EQ.1) GOTO 180
      IF(ITESUB.GT.ITEMAX) GOTO 175
      GKVOT=GNORM2/GNOLD2
C
      DO 50 I=1,M
      DSRCH(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 50
      DSRCH(I)=GRADF(I)+GKVOT*DSRCH(I)
   50 CONTINUE
C
      IF(ITESUB.LE.2) GOTO 15
C
      CALL HESSI(M,N,ULAM,HESSF,X,Y,ALFA,BETA,
     1           A,B,C,P,Q,P0,Q0,XLOW,XUPP,IYFREE)
C
      IK=0
      IKRED=0
      DO 70 K=1,M
      DO 65 I=K,M
      IK=IK+1
      IF(IYFREE(K).EQ.0) GOTO 65
      IF(IYFREE(I).EQ.0) GOTO 65
      IKRED=IKRED+1
      HESSF(IKRED)=HESSF(IK)
 65   CONTINUE
 70   CONTINUE
C
      HTRACE=0.
      IKRED=0
      ZZZZ=0.
      DO 73 K=1,NYDIM
      DO 72 I=K,NYDIM
      IKRED=IKRED+1
      IF(I.EQ.K) HTRACE=HTRACE+HESSF(IKRED)
      IF(I.EQ.K) ZZZZ=ZZZZ+1.
 72   CONTINUE
 73   CONTINUE
C
      HESMOD=0.0001*HTRACE/ZZZZ
      IF(HESMOD.LT.GEPS) HESMOD=GEPS
      IKRED=0
      DO 77 K=1,NYDIM
      DO 76 I=K,NYDIM
      IKRED=IKRED+1
      IF(I.EQ.K) HESSF(IKRED)=HESSF(IKRED)+HESMOD
 76   CONTINUE
 77   CONTINUE
C
      CALL LDLFAC(NYDIM,GEPS,HESSF,UU)
C
      IRED=0
      DO 79 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 79
      IRED=IRED+1
      UU(IRED)=GRADF(I)
 79   CONTINUE
C
      CALL LDLSOL(NYDIM,UU,HESSF,DSRCH)
C
      DO 80 I=1,M
      UU(I)=DSRCH(I)
 80   CONTINUE
C
      IRED=0
      DO 85 I=1,M
      DSRCH(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 85
      IRED=IRED+1
      DSRCH(I)=UU(IRED)
 85   CONTINUE
C
      GOTO 15
C
 175  CONTINUE
      WRITE(*,911)
C
 180  CONTINUE
C
      DO 90 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 90
      IF(ULAM(I).GE.0.) GOTO 90
      ULAM(I)=0.
   90 CONTINUE
C
 911  FORMAT(' WARNING IN SUBROUTINE SUBSPA')
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE HESSI(M,N,ULAM,HESSF,X,Y,ALFA,BETA,
     1                 A,B,C,P,Q,P0,Q0,XLOW,XUPP,IYFREE)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C   HESSI calculates HESSF = minus the reduced Hessian matrix of the
C   dual objective function, given the vector ULAM of dual variables.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ULAM(*),HESSF(*),X(*),Y(*),ALFA(*),BETA(*),A(*),
     1          B(*),C(*),P(*),Q(*),P0(*),Q0(*),XLOW(*),XUPP(*)
      INTEGER IYFREE(*)
C
      IK=0
      DO 12 K=1,M
      DO 11 I=K,M
      IK=IK+1
      HESSF(IK)=0.
 11   CONTINUE
 12   CONTINUE
C
      ULAMTA=0.
      II=1
      DO 15 I=1,M
      IF(I.GT.1) II=II+M+2-I
      IF(IYFREE(I).EQ.0) GOTO 15
      IF(ULAM(I).GT.C(I)) HESSF(II)=1.
      ULAMTA=ULAMTA+ULAM(I)*A(I)
 15   CONTINUE
C
      IF(ULAMTA.LE.1.) GOTO 40
C
      KK=1
      DO 30 K=1,M
      IF(K.GT.1) KK=KK+M+2-K
      IF(IYFREE(K).EQ.0) GOTO 30
      DO 20 I=K,M
      IF(IYFREE(I).EQ.0) GOTO 20
      IK=KK+I-K
      HESSF(IK)=HESSF(IK)+10.*A(I)*A(K)
 20   CONTINUE
 30   CONTINUE
C
 40   CONTINUE
C
      DO 100 J=1,N
      PJ=P0(J)
      QJ=Q0(J)
      MJ1=M*(J-1)
      DO 50 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 50
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      PJ=PJ+ULAM(I)*PIJ
      QJ=QJ+ULAM(I)*QIJ
 50   CONTINUE
C
      SRPJ=DSQRT(PJ)
      SRQJ=DSQRT(QJ)
      XJ=(SRPJ*XLOW(J)+SRQJ*XUPP(J))/(SRPJ+SRQJ)
      IF(XJ.GE.BETA(J)) GOTO 100
      IF(XJ.LE.ALFA(J)) GOTO 100
C
      UJXJ=XUPP(J)-XJ
      XJLJ=XJ-XLOW(J)
      UJXJ2=UJXJ**2
      XJLJ2=XJLJ**2
      RR=2.*PJ/UJXJ**3+2.*QJ/XJLJ**3
C
      KK=1
      DO 80 K=1,M
      IF(K.GT.1) KK=KK+M+2-K
      IF(IYFREE(K).EQ.0) GOTO 80
      PKJ=P(MJ1+K)
      QKJ=Q(MJ1+K)
      TTK=PKJ/UJXJ2-QKJ/XJLJ2
      DO 70 I=K,M
      IF(IYFREE(I).EQ.0) GOTO 70
      IK=KK+I-K
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      TTI=PIJ/UJXJ2-QIJ/XJLJ2
      HESSF(IK)=HESSF(IK)+TTI*TTK/RR
 70   CONTINUE
 80   CONTINUE
C
 100  CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE LDLFAC(N,EPS,ADL,E)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C    LDLFAC makes a factorization of a given symmetric matrix A.
C    If A is positive definite, then A = L*D*LT.
C    If A is not positive definite, then A + E = L*D*LT,
C    where E is a positive semidefinite diagonal matrix such that
C    A + E is positive definite.
C    On entry, ADL defines the given matrix A.
C    On leave, ADL defines the calculated matrices D and L.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ADL(*),E(*)
C
      JJ=1
C
      DO 100 J=1,N
      E(J)=0.
      IF(J.GT.1) JJ=JJ+N+2-J
      IF(J.EQ.1) GOTO 25
      KK=JJ
      JK=JJ
      DO 20 L=1,J-1
      K=J-L
      JK=JK-N+K
      KK=KK-N+K-1
      ADL(JJ)=ADL(JJ)-ADL(KK)*ADL(JK)*ADL(JK)
 20   CONTINUE
 25   IF(ADL(JJ).GE.EPS) GOTO 35
      E(J)=EPS-ADL(JJ)
      ADL(JJ)=EPS
 35   IF(J.EQ.N) GOTO 100
      IJ=JJ
      DO 50 I=J+1,N
      IJ=IJ+1
      ADLIJ=ADL(IJ)
      IF(J.EQ.1) GOTO 45
      IK=IJ
      JK=JJ
      KK=JJ
      DO 40 L=1,J-1
      K=J-L
      IK=IK-N+K
      JK=JK-N+K
      KK=KK-N+K-1
      ADLIJ=ADLIJ-ADL(KK)*ADL(IK)*ADL(JK)
 40   CONTINUE
 45   ADL(IJ)=ADLIJ/ADL(JJ)
 50   CONTINUE
 100  CONTINUE
C
      RETURN
      END
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE LDLSOL(N,B,DL,X)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     LDLSOL solves a system of linear equations: A*X = B,
C     where A has already been factorized as L*D*Ltranspose.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(*),DL(*),X(*)
C
      JJ=1
C
      DO 30 J=1,N
      X(J)=B(J)
      IF(J.EQ.1) GOTO 30
      JJ=JJ+N+2-J
      JK=JJ
      DO 20 L=1,J-1
      K=J-L
      JK=JK-N+K
      X(J)=X(J)-DL(JK)*X(K)
 20   CONTINUE
 30   CONTINUE
C
      JJ=1
      DO 40 J=1,N
      IF(J.GT.1) JJ=JJ+N+2-J
      X(J)=X(J)/DL(JJ)
 40   CONTINUE
C
      DO 60 L=1,N-1
      J=N-L
      JJ=JJ-N+J-1
      KJ=JJ
      DO 50 K=J+1,N
      KJ=KJ+1
      X(J)=X(J)-DL(KJ)*X(K)
 50   CONTINUE
 60   CONTINUE
C
      RETURN
      END
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE XYZLAM(M,N,X,Y,Z,ULAM,XLOW,XUPP,ALFA,BETA,
     1                  A,B,C,P,Q,P0,Q0,IYFREE)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     XYZLAM calculates the X,Y,Z that minimize the Lagrange
C     function, given the vector ULAM of Lagrange multipliers.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),Y(*),ULAM(*),XLOW(*),XUPP(*),ALFA(*),BETA(*),
     1          A(*),B(*),C(*),P(*),Q(*),P0(*),Q0(*)
      INTEGER IYFREE(*)
C
      DO 30 J=1,N
      PJ=P0(J)
      QJ=Q0(J)
      MJ1=M*(J-1)
C
      DO 20 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 20
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      PJ=PJ+ULAM(I)*PIJ
      QJ=QJ+ULAM(I)*QIJ
 20   CONTINUE
C
      SRPJ=DSQRT(PJ)
      SRQJ=DSQRT(QJ)
      XJ=(SRPJ*XLOW(J)+SRQJ*XUPP(J))/(SRPJ+SRQJ)
      IF(XJ.LT.ALFA(J)) XJ=ALFA(J)
      IF(XJ.GT.BETA(J)) XJ=BETA(J)
      X(J)=XJ
 30   CONTINUE
C
      UA=0.
      DO 40 I=1,M
      Y(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 40
      UA=UA+ULAM(I)*A(I)
      YI=ULAM(I)-C(I)
      IF(YI.GT.0.) Y(I)=YI
 40   CONTINUE
C
      Z=0.
      UA1=UA-1.
      IF(UA1.GT.0.) Z=10.*UA1
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA,
     1                 A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     GRADI calculates the gradient GRADF of the dual
C     objective function, given the vector ULAM of dual variables.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),Y(*),ULAM(*),XLOW(*),XUPP(*),ALFA(*),BETA(*),
     1          A(*),B(*),C(*),P(*),Q(*),P0(*),Q0(*),GRADF(*)
      INTEGER IYFREE(*)
C
      CALL XYZLAM(M,N,X,Y,Z,ULAM,XLOW,XUPP,ALFA,BETA,
     1            A,B,C,P,Q,P0,Q0,IYFREE)
C
      DO 10 I=1,M
      GRADF(I)=-B(I)-Y(I)-A(I)*Z
 10   CONTINUE
C
      DO 30 J=1,N
      MJ1=M*(J-1)
      UJXJ=XUPP(J)-X(J)
      XJLJ=X(J)-XLOW(J)
C
      DO 20 I=1,M
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      GRADF(I)=GRADF(I)+PIJ/UJXJ+QIJ/XJLJ
 20   CONTINUE
 30   CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE LINDER(M,N,T,DFDT,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1                  ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     LINDER calculates the scalar product DFDT of GRADF and DSRCH
C     (= the directional derivative) at the point ULAM + T*DSRCH.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ULAM(*),DSRCH(*),X(*),Y(*),UU(*),
     1          XLOW(*),XUPP(*),ALFA(*),BETA(*),
     2          A(*),B(*),C(*),P(*),Q(*),P0(*),Q0(*)
      INTEGER IYFREE(*)
C
      DO 10 I=1,M
      UU(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 10
      UU(I)=ULAM(I)+T*DSRCH(I)
      IF(UU(I).LT.0.) UU(I)=0.
 10   CONTINUE
C
      CALL XYZLAM(M,N,X,Y,Z,UU,XLOW,XUPP,ALFA,BETA,
     1            A,B,C,P,Q,P0,Q0,IYFREE)
C
      DO 20 I=1,M
      UU(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 20
      UU(I)=-B(I)-Y(I)-A(I)*Z
   20 CONTINUE
C
      DO 40 J=1,N
      MJ1=M*(J-1)
      UJXJ=XUPP(J)-X(J)
      XJLJ=X(J)-XLOW(J)
      DO 30 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 30
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      UU(I)=UU(I)+PIJ/UJXJ+QIJ/XJLJ
   30 CONTINUE
   40 CONTINUE
C
      DFDT=0.
      DO 50 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 50
      DFDT=DFDT+UU(I)*DSRCH(I)
   50 CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE LINESE(M,N,TMAX,TOPT,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1                  ALFA,BETA,A,B,C,P,Q,P0,Q0,
     2                  IYFREE,IHITY,IHITMX,ITESUB)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     LINESE makes an approximate line search (maximization) in the
C     direction DSRCH from the point ULAM.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ULAM(*),DSRCH(*),X(*),Y(*),UU(*),XLOW(*),XUPP(*),
     1          ALFA(*),BETA(*),A(*),B(*),C(*),P(*),Q(*),P0(*),Q0(*)
      INTEGER IYFREE(*)
C
      ITT1=0
      ITT2=0
      ITT3=0
C
      IF(IHITY.EQ.0) GOTO 40
      CALL LINDER(M,N,TMAX,DFDTMX,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDTMX.GT.0.) GOTO 80
      IF(TMAX.GT.1.) GOTO 40
      T2=TMAX
C
 30   ITT1=ITT1+1
      IF(ITT1.GT.12) GOTO 90
      T1=T2/2.
      IF(ITESUB.LE.3) T1=T2/256.
      CALL LINDER(M,N,T1,DFDT1,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDT1.GT.0.) GOTO 60
      T2=T1
      GOTO 30
C
 40   T1=1.
      T2=T1
      CALL LINDER(M,N,T1,DFDT1,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(ITESUB.LE.5.AND.DFDT1.GT.0.) GOTO 50
      IF(ITESUB.GE.6.AND.DFDT1.GE.0.) GOTO 90
      GOTO 30
C
 50   ITT2=ITT2+1
      IF(ITT2.GT.10) GOTO 90
      T2=2.*T1
      IF(ITESUB.LE.3) T2=256.*T1
      IF(IHITY.EQ.0) GOTO 55
      IF(T2.LT.TMAX) GOTO 55
      T2=TMAX
      GOTO 60
 55   CALL LINDER(M,N,T2,DFDT2,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDT2.LE.0.) GOTO 60
      T1=T2
      GOTO 50
C
 60   CONTINUE
      SQT1=DSQRT(T1)
      SQT2=DSQRT(T2)
 62   ITT3=ITT3+1
      IF(ITT3.GT.9) GOTO 90
      TM=SQT1*SQT2
      CALL LINDER(M,N,TM,DFDTM,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDTM.GT.0.) GOTO 65
      T2=TM
      TKVOT=T1/T2
      IF(TKVOT.GT.0.97) GOTO 90
      SQT2=DSQRT(T2)
      GOTO 62
 65   T1=TM
      TKVOT=T1/T2
      IF(TKVOT.GT.0.97) GOTO 90
      SQT1=DSQRT(T1)
      GOTO 62
C
 80   TOPT=TMAX
      IHITMX=1
      GOTO 100
 90   TOPT=T1
      IHITMX=0
 100  CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE MMASUB(ITER,M,N,GEPS,IYFREE,XVAL,XMMA,
     1                  XMIN,XMAX,XOLD1,XOLD2,XLOW,XUPP,
     2                  ALFA,BETA,A,B,C,Y,Z,ULAM,
     3                  F0VAL,FVAL,FMAX,DF0DX,DFDX,
     4                  P,Q,P0,Q0,UU,GRADF,DSRCH,HESSF)
C
C       Version "November 2000".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C    Use of this code is for academic purposes only,
C    regulated by an agreement with Krister Svanberg.
C    The code is not to be redistributed.
C
C    MMASUB generates and solves the MMA subproblem,
C    which is of the following form in the variables
C    x_1,...,x_N, y_1,...,y_M, and z.
C
C   minimize h_0(x) + r_0 + z + 0.05*z^2 + sum{c_i*y_i + 0.5*(y_i)^2}
C
C subject to h_i(x) - a_i*z - y_i <= b_i ,     i=1,..,M
C                   alfa_j <= x_j <= beta_j ,  j=1,..,N
C                             y_i >= 0 ,       i=1,..,M
C                               z >= 0 .
C
C    with h_i(x) = sum{p_ij/(xupp_j-x_j) + q_ij/(x_j-xlow_j)}.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION XVAL(*),XMMA(*),XMIN(*),XMAX(*),
     1          XOLD1(*),XOLD2(*),XLOW(*),XUPP(*),
     2          ALFA(*),BETA(*),A(*),B(*),C(*),Y(*),ULAM(*),
     3          FVAL(*),FMAX(*),DF0DX(*),DFDX(*),
     4          P(*),Q(*),P0(*),Q0(*),UU(*),
     5          GRADF(*),DSRCH(*),HESSF(*)
      INTEGER IYFREE(*)
C
C********+*********+*********+*********+*********+*********+*********+
C  The sizes of the above areas must be at least as follows:
C
C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
C    IYFREE(M)
C********+*********+*********+*********+*********+*********+*********+
C
C***  Input to the subroutine MMASUB:
C
C  ITER  = Current iteration number ( =1 the first iteration).
C     N  = Number of variables x_j in the problem.
C     M  = Number of constraints in the problem (not including
C          the simple upper and lower bounds on the variables).
C  GEPS  = Tolerance parameter for the constraints.
C          (Used in the termination criteria for the subproblem.)
C   XVAL(j) = Current value of the variable x_j.
C   XMIN(j) = Original lower bound for the variable x_j.
C   XMAX(j) = Original upper bound for the variable x_j.
C  XOLD1(j) = Value of the variable x_j one iteration ago.
C  XOLD2(j) = Value of the variable x_j two iterations ago.
C   XLOW(j) = Current value of the lower asymptot l_j.
C   XUPP(j) = Current value of the upper asymptot u_j.
C      A(i) = Coefficient a_i for the minimax variable z.
C      C(i) = Coefficient c_i for the artificial variable y_i.
C    F0VAL  = Value of the objective function f_0(x)
C   FVAL(i) = Value of the i:th constraint function f_i(x).
C   FMAX(i) = Right hand side of the i:th constraint.
C  DF0DX(j) = Derivative of f_0(x) with respect to x_j.
C   DFDX(k) = Derivative of f_i(x) with respect to x_j,
C             where k = (j-1)*M + i.
C
C*** Output from the subroutine MMASUB:
C
C   XMMA(j) = Optimal value of x_j in the MMA subproblem.
C      Y(i) = Optimal value of the "artificial" variable y_i.
C      Z    = Optimal value of the "minimax" variable z.
C   ULAM(i) = Optimal value of the dual variable lambda_i.
C   XLOW(j) = New values of the lower asymptot l_j.
C   XUPP(j) = New values of the upper asymptot u_j.
C
C*** Working areas and their usage in MMASUB:
C
C   ALFA(j) = Lower bound for x_j in the MMA subproblem.
C   BETA(j) = Upper bound for x_j in the MMA subproblem.
C      P(k) = Coefficient p_ij in the MMA subproblem,
C             where k = (j-1)*M + i.
C      Q(k) = Coefficient q_ij in the MMA subproblem,
C             where k = (j-1)*M + i.
C     P0(j) = Coefficient p_0j in the MMA subproblem.
C     Q0(j) = Coefficient q_0j in the MMA subproblem.
C      B(i) = Right hand side b_i in the MMA subproblem.
C  GRADF(i) = Gradient component of the dual objective function.
C  DSRCH(i) = Search direction component in the dual subproblem.
C  HESSF(k) = Hessian matrix component of the dual function.
C     UU(i) = Component in a working area.
C IYFREE(i) = 0 for dual variables which are fixed to zero in
C               the current subspace of the dual subproblem,
C           = 1 for dual variables which are "free" in
C               the current subspace of the dual subproblem.
C
C********+*********+*********+*********+*********+*********+*********+
C
      CALL ASYMPT(ITER,M,N,XVAL,XMIN,XMAX,XOLD1,XOLD2,
     1            XLOW,XUPP,ALFA,BETA)
C
C****  ASYMPT calculates the asymptotes XLOW(j) and XUPP(j),
C****  and the bounds ALFA(j) and BETA(j).
C
      CALL GENSUB(M,N,XVAL,XMIN,XMAX,F0VAL,DF0DX,FMAX,FVAL,
     1            DFDX,P,Q,B,P0,Q0,R0,XLOW,XUPP)
C
C***** GENSUB generates the MMA subproblem by calculating the
C***** coefficients P(i,j),Q(i,j),B(i),P0(j),Q0(j) and R0.
C     
      CALL MAXIM(M,N,GEPS,IYFREE,GRADF,DSRCH,HESSF,XMMA,Y,Z,
     1           ULAM,UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
C
C***** MAXIM solves the MMA subproblem by maximizing the dual
C***** objective function subject to non-negativity constraints
C***** on the dual variables.
C***** ULAM = optimal dual solution of the MMA subproblem.
C***** XMMA,Y,Z = optimal (primal) solution of the MMA subproblem.
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
