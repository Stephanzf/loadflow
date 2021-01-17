	PROGRAM LOADFLOW
	PARAMETER (NBU=65,NG=50,IPQLFI=15, ILDFLO=16,
     *                           IPQLFO=17, IZLEQO=18 )
	COMMON/S15/Y(NBU,NBU),BUSTYPE(NBU),   PG(NBU),
     *                QG(NBU), PL(NBU),QL(NBU),VSPEC(NBU),
     *                QMIN(NBU),QMAX(NBU) ,VMAG(NBU), VANG(NBU),
     *                DF(NBU)  ,DX(NBU)   ,A(NBU,NBU)
	COMMON/S16/B1P(NBU,NBU), B11P(NBU,NBU)
	COMMON/S17/GENBUS( NG ),V( NBU )

	INTEGER        BUSTYPE   , NBUS    , NGEN  , PUSW
	INTEGER        GENBUS    , ITERMAX , NWRDSW
	REAL           ERRMAXPQ
	COMPLEX        Y         , V
	COMPLEX        CMPLX     , CONJG
	COMPLEX        YIJ       , SUM

	OPEN( IPQLFI, FILE= ' PQLF.INP '  , STATUS = 'OLD' )
	OPEN( IPQLFO, FILE= ' PQLF.OUT '  )
	OPEN( ILDFLO, FILE= ' LDFL.OUT '  )
	OPEN( IZLEQO, FILE= ' ZLEQ.OUT '  )

C    -------------------------
	CALL STRING (IPQLFI, ILDFLO )
100     READ(IPQLFI, *, ERR =110, END =110 ) NBUS,
     *                            ITERMAX,  ERRMAXPQ , NWRDSW ,PUSW
500     FORMAT(T8,I3 ,T19,I3 , T33,F9.6 , T49, I2, T61, I2/)
	
110     NWRDSW = 1 
	write(*,*)
	write(*, '(A45)') '  Enter 0 -->              LF study , '
	write(*,*)
	write(*, '(A45)') '  Any Key --> LF et NTWK Equivalent ! '
	read (*,*, ERR = 8) NWRDSW
8       continue 
	WRITE(ILDFLO, 500 ) NBUS, ITERMAX,  ERRMAXPQ , NWRDSW ,PUSW  

	DO  I = 1 , NBUS
	DO  J = 1 , NBUS
	B1P(I,J)  = 0.0
	B11P(I,J) = 0.0
	Y(I,J)    = CMPLX(0.0,0.0)
	END DO
	END DO
C
C   loop here read each line data & assembly Y matrix:
C   read each line data & assembly  Y, B1P, B11P matrices
C
	CALL STRING (IPQLFI, ILDFLO )
165     READ(IPQLFI, *, ERR=180, END=180 ) I, J, R, R2, X, X2, BK
180     IF (I.EQ.0) GOTO 5
	WRITE(ILDFLO, 181 ) I, J, R, R2, X, X2, BK
181     FORMAT(T5, I4, T10, I4, T17, F8.4, T26, F8.4, T35, F8.4,
     *                                     T44, F8.4, T55, F8.4  )
C
C   line branch assembly:
C        
	IF (J.GT.0) THEN
	  YIJ=(1.0,0.0)/CMPLX(R,X)
	  NUL=0.0
	  Y(I,I)=Y(I,I)+YIJ+CMPLX(NUL,BK/2.0)
	  Y(J,J)=Y(J,J)+YIJ+CMPLX(NUL,BK/2.0)
	  Y(I,J)=Y(I,J)-YIJ
	  Y(J,I)=Y(J,I)-YIJ
	  B1P(I,J) = B1P(I,J) - AIMAG(YIJ)
	  B1P(J,I) = B1P(I,J)
	  B1P(I,I) = B1P(I,I) + AIMAG(YIJ)
	  B1P(J,J) = B1P(J,J) + AIMAG(YIJ)
	END IF
C
C   transformer line assembly:
C        
	IF (J.LT.0) THEN
	  J        = -J
	  YIJ      = (1.0,0.0)/CMPLX(R,X)
	  Y(I,I)   = Y(I,I)+YIJ
	  Y(J,J)   = Y(J,J)+YIJ/BK/BK
	  Y(I,J)   = Y(I,J)-YIJ/BK
	  Y(J,I)   = Y(J,I)-YIJ/BK
	  B1P(I,J) = B1P(I,J) - AIMAG(YIJ)/BK
	  B1P(J,I) = B1P(I,J)
	  B1P(I,I) = B1P(I,I) + AIMAG(YIJ)/BK
	  B1P(J,J) = B1P(J,J) + AIMAG(YIJ)/BK
	END IF
C
C   grounded line assembly:
C        
	IF (J.EQ.0) THEN
	  IF (BK.EQ.1.0) THEN
	  YIJ      = (1.0,0.0)/CMPLX(R,X)
	  Y(I,I)   = Y(I,I)+YIJ
	  B1P(I,I) = B1P(I,I) + AIMAG(YIJ)
	  END IF
	  IF (BK.EQ.-1.0) THEN
	  YIJ      = CMPLX(R,X)
	  Y(I,I)   = Y(I,I)+YIJ
	  B1P(I,I) = B1P(I,I) + AIMAG(YIJ)
	  END IF
	END IF
C
C    go back to read next line data
C        
	GOTO 165
5       CONTINUE
C
C    initialize bus type & voltage unkown:
C        
	DO I = 1 , NBUS
	BUSTYPE(I) = 1
	VMAG(I) = 1.0
	VANG(I) = 0.0
	END DO

	NPVBUS  = 0
C
C   read in bus data
C
	CALL STRING (IPQLFI, ILDFLO )
2       READ(IPQLFI, *, END=280 ) I, I1, P1, Q1, P2, Q2, Vspc
280     IF ( I.EQ.0 ) GOTO 6
	BUSTYPE(I) = I1
	PG(I)    = P1
	QG(I)    = Q1
	PL(I)    = P2
	QL(I)    = Q2
	VSPEC(I) = Vspc
	WRITE(ILDFLO, 281) I,BUSTYPE(I),PG(I),QG(I), PL(I), 
     *                                       QL(I), VSPEC(I) 
281     FORMAT(T5, I4, T13, I2,  T22, F8.3,
     *         T31, F8.3, T40, F8.3, T49, F8.3, T60, F6.3  )
	IF (BUSTYPE(I).EQ.2) NPVBUS = NPVBUS + 1
	IF((BUSTYPE(I).EQ.2).OR.(BUSTYPE(I).EQ.3)) VMAG(I) = VSPEC(I)
	
	GOTO 2
6       CONTINUE
	CALL STRING (IPQLFI, ILDFLO )
C       
C  form B' matrix and facterize it:
C        
	DO I = 1, NBUS
	  IF (BUSTYPE(I).EQ.3)  THEN
	    DO 20 K = I, NBUS-1
	    DO 20 J = 1, NBUS
20          B1P(K,J) = B1P(K+1,J)
	    DO 30 K = I, NBUS-1
	    DO 30 L = 1, NBUS
30          B1P(L,K) = B1P(L,K+1)
	  END IF
	  GOTO 35
	END DO
C
C   form B" matrix and facterize it:
C
35      L = 1
	DO I = 1, NBUS
	  IF (BUSTYPE(I).EQ.1) THEN
	    K = 1
	    DO J = 1, NBUS
	    IF (BUSTYPE(J).EQ.1) THEN
	    B11P(L,K) = AIMAG(Y(I,J))
	    K = K + 1
	    END IF
	    END DO
	    L = L +1
	  END IF
	END DO
	CALL FACTERIZE( B1P , NBUS - 1 )   
	CALL FACTERIZE( B11P, NBUS-1-NPVBUS )
C                          
C   iteration begins here
C
	ITER = 0
	KP   = 1
	KQ   = 1
	WRITE(*,678)
	WRITE(ildflo,678)
678     format(/1x,'iter',1x, T7 ,  'bus#', T16, 'DQmax' ,
     +                       T24,  'bus#', T33, 'DPmax' /   
     +  ' -------------------------------------')
C
C   convert generation and demand to per unit on 100MVA base
C
	DO I = 1, NBUS
	IF ( PUSW.NE.0 ) THEN
	   PG( I ) = PG( I ) / 100.0
	   QG( I ) = QG( I ) / 100.0
	   PL( I ) = PL( I ) / 100.0
	   QL( I ) = QL( I ) / 100.0
	END IF
	END DO

7777    CONTINUE
	IF ( ITER.GT.ITERMAX ) THEN
	  WRITE(*,*) 'P-Q Decoupled LF diverged !'
	  STOP
	END IF
	K    = 1        
	DO I = 1, NBUS   ! I loop begins
	IF (BUSTYPE(I).NE.3) THEN
	  X1 = 0.0
	  DO J = 1, NBUS
	     XIJ = VANG(I) - VANG(J)
	     X1  = X1 + VMAG(J)*( REAL(Y(I,J))*COS(XIJ) 
     *               + AIMAG(Y(I,J))*SIN(XIJ) )
	  END DO
	  DF(K) = (PG(K)-PL(K))/VMAG(K) - X1
	  K = K + 1
	END IF 
	END DO            
	DPMAX = 0.0
	DO I = 1, NBUS - 1
	IF ( DABS(DF(I)).GT.ERRMAXPQ ) THEN
	     DPMAX  = DABS(DF(I))
	     IDPMAX = I
	END IF
	END DO

	IF ( DPMAX.LE.ERRMAXPQ ) THEN
	  KP = 0
	  IF ( KQ.EQ.0 ) GOTO       2121
	  GOTO                      1212
	END IF
	CALL AXEQB( B1P, DX, DF, NBUS-1 )
	DO I = 1, NBUS - 1
	VANG(I) = VANG(I) -DX(I)/VMAG(I)
	END DO
	ITER = ITER + 1
	KQ   = 1
1212    K    =1
	DO I = 1, NBUS   ! I loop begins
	IF (BUSTYPE(I).EQ.1) THEN
	  X1 = 0.0
	  DO J = 1, NBUS
	     XIJ = VANG(I) - VANG(J)
	     X1  = X1 + VMAG(J)*(  REAL(Y(I,J))*SIN(XIJ) 
     *                             - AIMAG(Y(I,J))*COS(XIJ) )
	  END DO
	  DF(K) = (QG(I)-QL(I))/VMAG(I) - X1
	  K = K + 1
	END IF 
	END DO           
	DQMAX = 0.0

	DO I = 1, NBUS - 1 -NPVBUS
	IF ( DABS(DF(I)).GT.ERRMAXPQ ) THEN
	     DQMAX  = DABS(DF(I))
	     IDQMAX = I
	END IF
	END DO
	it = int( iter/2 )
	write(ILDFLO,677) it, idqmax, dqmax , idpmax, dpmax
	write(*,677) it, idqmax, dqmax , idpmax, dpmax
677     format( 2x, i3,1x, T7,i4, T12, F10.6, T24, i4,T29,f10.6 )

	IF ( DQMAX.LE.ERRMAXPQ ) THEN
	  KQ = 0
	  IF ( KP.EQ.0 ) GOTO       2121
	  GOTO                      7777
	END IF

	CALL AXEQB( B11P, DX, DF, NBUS-1-NPVBUS )
	L1 = 1
	DO I = 1, NBUS  
	IF (BUSTYPE(I).EQ.1) THEN
	    VMAG(I) = VMAG(I) -DX(L1)
	    L1 = L1 + 1
	END IF
	END DO
	ITER = ITER + 1
	KP   = 1
	GOTO                         7777
2121    CONTINUE
C   calculate reactive power for pv buses
	DO I = 1, NBUS
	IF ( BUSTYPE(I).EQ.2) THEN
	   SUM  = ( 0.0, 0.0 )
	   DO J = 1, NBUS
	   YIJ  = CMPLX( VMAG(J)*COS(VANG(J)), VMAG(J)*SIN(VANG(J)) )
	   SUM  = SUM +CONJG(Y(I,J))*CONJG(YIJ)
	   END DO
	   SUM=SUM*CMPLX( VMAG(I)*COS(VANG(I)),VMAG(I)*SIN(VANG(I)))  
	   if( (PG(I).NE.0.0).OR.(QG(I).NE.0.0) ) then
	       QG(I) = AIMAG( SUM ) + QL(I)
	   else 
	       QL(I) = - AIMAG( SUM )
	   end if
	END IF
	END DO
C   calculate real and reactive power for the slack bus        
	DO I = 1, NBUS 
	IF ( BUSTYPE(I).EQ.3) THEN
	   SUM  = ( 0.0, 0.0 )
	   DO J = 1, NBUS
	   YIJ  = CMPLX( VMAG(J)*COS(VANG(J)), VMAG(J)*SIN(VANG(J)) )
	   SUM  = SUM +CONJG(Y(I,J))*CONJG(YIJ)
	   END DO
	   SUM   =  SUM*VMAG(I)
	   PG(I) =  REAL( SUM ) + PL(I)
	   QG(I) = AIMAG( SUM ) + QL(I)
	   GOTO                      8888
	END IF
	END DO
8888    CONTINUE    
	
	WRITE(ILDFLO,540) IT
	WRITE(*,540)      IT
540     FORMAT(1X,'P-Q Dcoupled LF converged at',I3, 'th iteration!'/)
	IF ( NWRDSW.EQ.1 ) GOTO 555
	DO I = 1, NBUS
	  IF ( (PL(I).NE.0.0 ).OR.( QL(I).NE.0.0 ) ) THEN
	  SUM = CMPLX( PL(I), -QL(I) ) / VMAG(I)**2
	  WRITE( IZLEQO, 545 ) I, SUM 
	  END IF
	END DO
545     FORMAT( I5, 2F15.7 )

555     WRITE(ILDFLO,550)
550     FORMAT(3X,'I',1X,'TYP',5X,'PG',5X,5X,'QG',5X,5X,'PL',5X,
     *  5X,'QL',5X,3X,'VMAG', 2X,4X,'ANGL',4X/)

560     FORMAT(1X, I3, I3, 4F12.5, F9.5, F12.5)

C
C   convert generation and demand to nominal unit
C
	IF ( PUSW.NE.0 ) THEN
	DO I = 1, NBUS
	   PG( I ) = PG( I ) * 100.0
	   QG( I ) = QG( I ) * 100.0
	   PL( I ) = PL( I ) * 100.0
	   QL( I ) = QL( I ) * 100.0
	END DO
	END IF

	DO I = 1, NBUS
	VANG(I) = VANG(I)*180.0/3.141592
	WRITE(ILDFLO,560) I,BUSTYPE(I),PG(I),QG(I),PL(I),QL(I),
     *                   VMAG(I), VANG(I) 
	END DO
C
C   P, Q are printed out in MW and MVA on their bases, VMAG in p.u.
C
	IF ( NWRDSW.EQ.0 ) GOTO 600
	DO I = 1, NBUS
	BUSTYPE ( I ) = 0
	END DO
	IG = 1
	DO I = 1,  NBUS
	IF ( (PG(I).NE.0.0).OR.(QG(I).NE.0.0) ) THEN
	   IF ( PUSW.EQ.0 ) THEN
	      WRITE( IPQLFO, 575 ) IG, ( PG(I)-PL(I))*100.0,
     *                ( QG(I)-QL(I) )*100.0 ,  VMAG(I) , VANG(I)
	   ELSE
	      WRITE( IPQLFO, 575 ) IG, PG(I)-PL(I), QG(I)-QL(I),
     *                              VMAG(I)    , VANG(I)
	   END IF
	   BUSTYPE (I)  = 1
	   GENBUS  (IG) = I
	   IG = IG + 1
	END IF
	END DO
575     FORMAT( 1X, I4, 2X, 4F15.6 )
	NGEN = IG - 1
C
C   call network reduction subroutine to form various Ys
C
	IF ( PUSW.EQ.0 ) THEN
	DO I = 1, NBUS
	   PG( I ) = PG( I ) * 100.0
	   QG( I ) = QG( I ) * 100.0
	   PL( I ) = PL( I ) * 100.0
	   QL( I ) = QL( I ) * 100.0
	END DO
	END IF

	CALL NWRD( NBUS, NGEN , IPQLFO )
	WRITE (*, '(A34)' ) 'NTWK Equivalent is also executed.'
	write (*, *)
600     CONTINUE
	CLOSE ( IPQLFI )
	CLOSE ( IPQLFO )
	CLOSE ( ILDFLO )
	CLOSE ( IZLEQO )
	STOP
	END
C   This subroutine facterize a matrix: A = [ L D U ]
	SUBROUTINE FACTERIZE( A, NN)
	DIMENSION  A(65,65)
	REAL       A
	INTEGER    NN
	DO  I = 1, NN  
	   DO 1 K = I + 1, NN
	   DO 1 J = I + 1, NN
1          A(K,J) = A(K,J) - A(K,I)*A(I,J)/A(I,I)
	END DO
	DO 2 I = 1, NN
	   A(I,I) = 1.0/A(I,I)
	   DO 2 J  = I + 1, NN
2          A(I,J) = A(I,J)*A(I,I)
	RETURN
	END
C
C   This subroutine solves the eqn [ A ] [ X ] = [ B ]
C       
	SUBROUTINE  AXEQB(A, DX, DF, NN)
	DIMENSION   A(65, 65), DX(65),  DF(65)
	REAL        A, DX,  DF , X
	INTEGER     NN
C
C   forward process
	DF(1) = DF(1)*A(1,1)
	DO I = 2, NN
	X = 0.0
	  DO J = 1, I-1
	  X = X + A(I,J)*DF(J)
	  END DO
	DF(I) = ( DF(I) - X )*A(I,I)
	END DO
C
C   backward process
	DX(NN) = DF(NN)
	DO I  = NN - 1, 1, -1
	  X = 0.0
	  DO J = I + 1, NN
	  X = X + A(I,J)*DX(J)
	  END DO
	DX(I) = DF(I) - X
	END DO
	RETURN
	END

C   This subroutine reads character strings from datafile

	SUBROUTINE STRING( IPQLFI, ILDFLO )
	CHARACTER*1      STR
	CHARACTER*75     STRR

10      READ ( IPQLFI, 20 , END = 30)  STR
	IF ( STR.EQ. '!' ) THEN
	   BACKSPACE IPQLFI
	   READ ( IPQLFI, 22 )  STRR
	   WRITE( ILDFLO, 25 )  STRR
	   GOTO  10
	END IF 
20      FORMAT( A1  )
22      FORMAT( A75 )
25      FORMAT( 2X, A75 )
	BACKSPACE IPQLFI
30      RETURN
	END 
C
C   This subroutine is for network reduction according to
C         fault location and type concerned
C
	SUBROUTINE NWRD ( NBUS, NGEN , IPQLFO )
	PARAMETER ( NBU=65, NG=50,  IZGNDO=18 )
	COMMON/S15/Y(NBU,NBU),BUSTYPE(NBU),   PG(NBU),
     *                QG(NBU), PL(NBU),QL(NBU),VSPEC(NBU),
     *                QMIN(NBU),QMAX(NBU) ,VMAG(NBU), VANG(NBU),
     *                DF(NBU)  ,DX(NBU)   ,A(NBU,NBU)
	COMMON/S17/GENBUS( NG ),V( NBU )
	dimension   Y0( NBU, NBU )
	dimension   BUSTYPE0(NBU)

	INTEGER        BUSTYPE , NBUS ,  NGEN , GENBUS  ,  FLTBUS
	COMPLEX        CMPLX   , CONJG,  Y    , Y0      ,  V 
	COMPLEX        CURR    , Zeq
	REAL           VMAG    , VANG ,  PL   , QL , PG ,  QG
	INTEGER        BUSTYPE0 
	common /sss /Zeq
C
C   open file created by the STCC program        
C
	OPEN (IZGNDO, FILE = ' ZGND.OUT ' , STATUS = ' OLD ' )

	DO I=1,NBUS
	VANG(I) = VANG( I ) * 3.141592/180.0
	V(I)    = VMAG( I ) * CMPLX( COS(VANG(I)), SIN(VANG(I)) )    
	Y(I,I)  = Y(I,I) + CMPLX( PL(I),-QL(I) )*0.01/VMAG(I)**2
	END DO
	
	DO 10 I = 1, NBUS
	BUSTYPE0(I) = BUSTYPE(I)
	DO 10 J = 1, NBUS
10      Y0( I, J ) = Y( I, J )  
C
C   read fault information output from STCC program
C
	LL = 0
5       READ (IZGNDO, 280, END = 7 ) XI, XJ, XL, ST, CT, RT, Zeq
	LL = LL +  1
	GOTO    5
280     FORMAT( 8F10.4 )
7       NFAULT = LL 
	FLTBUS = INT( XI )

	IF ( FLTBUS.EQ.0 ) THEN
	   CALL   KRON ( NBUS , 0  )   
	   CALL PYTIME ( NGEN , 50.0, IPQLFO )
C
C   check if reduction is correct
840        WRITE( *, 850 )
850        FORMAT( 2X, 'COMPUTED GENERATION'/)
	   DO I = 1 , NGEN
	   CURR = ( 0.0 , 0.0 )
	     DO J = 1 , NGEN
	     CURR = CURR + Y( GENBUS(I),GENBUS(J) )*V( GENBUS(J) )
	     ENDDO
	   CURR = V( GENBUS(I) ) * CONJG(CURR)
	   WRITE( * , 900 ) GENBUS(I) , CURR
	   END DO
900        FORMAT(2X,I5,2P2F14.4)
	   RETURN 
	END IF
	IF ( ST.NE.0.0 ) THEN
	   CALL  KRON  ( NBUS , 0 )   
	   CALL PYTIME ( NGEN , ST, IPQLFO )
	END IF
C        
C   start the fault and modify the healthy Y ( Y0 ) matrix
C
	   DO 20  I = 1, NBUS
	   BUSTYPE(I) = BUSTYPE0(I)
	   DO 20  J = 1, NBUS
20         Y( I,J ) = Y0( I, J )  
	IF ( CABS(Zeq).EQ.0.0 )  THEN
	   DO I=1,NBUS
	   Y( I , FLTBUS ) = ( 0.0, 0.0 )
	   Y( FLTBUS, I  ) = ( 0.0, 0.0 )
	   ENDDO
	ELSE
	   Y( FLTBUS, FLTBUS ) = Y( FLTBUS, FLTBUS ) + 1.0 / Zeq
	END IF
	CALL  KRON  ( NBUS , FLTBUS ) 
	CALL PYTIME ( NGEN , CT, IPQLFO )
C        
C   restore to the healthy Y0 matrix used upto end of simulation
C
	DO 30  I = 1, NBUS
	BUSTYPE(I) = BUSTYPE0(I)
	DO 30  J = 1, NBUS
30      Y( I,J ) = Y0( I, J )  
	CALL  KRON  ( NBUS , 0  ) 
	CALL PYTIME ( NGEN , 50.0, IPQLFO )
C
	CLOSE(IZGNDO)
	RETURN
	END
C
C   This subroutine is a Kron elimilation process  
C
	SUBROUTINE KRON ( NBUS, FLTBUS )
	PARAMETER ( NBU=65 )
	COMMON/S15/Y(NBU,NBU),BUSTYPE(NBU),   PG(NBU),
     *                QG(NBU), PL(NBU),QL(NBU),VSPEC(NBU),
     *                QMIN(NBU),QMAX(NBU) ,VMAG(NBU), VANG(NBU),
     *                DF(NBU)  ,DX(NBU)   ,A(NBU,NBU)

	INTEGER       BUSTYPE  , FLTBUS
	COMPLEX       Y        , Zeq
	common /sss /Zeq

	DO K = 1, NBUS
	IF( (  (K.EQ.FLTBUS).AND.(CABS(Zeq).EQ.0.0)  ).OR.
     *      (   BUSTYPE(K).EQ.1  ).OR.
     *      ( CABS(Y(K,K)).EQ.0.0)    )  GOTO   200
	BUSTYPE(K) = -1
	DO I = 1 , NBUS
	 IF (BUSTYPE(I).NE.-1)  THEN
	   DO J = 1 , NBUS
	   IF(BUSTYPE(J).NE.-1)   THEN
	   Y( I,J ) = Y( I,J ) - Y( I,K ) * Y( K,J ) / Y( K,K )
	   ENDIF
	   END DO
	 END IF
	END DO
200     END DO
	RETURN
	END
C           
C   This subroutine prints various Ys and TIME to a file
C
	SUBROUTINE PYTIME ( NGEN , XT, IPQLFO )
	PARAMETER ( NBU=65 , NG=50 )
	COMMON/S15/Y(NBU,NBU),BUSTYPE(NBU),   PG(NBU),
     *                QG(NBU), PL(NBU),QL(NBU),VSPEC(NBU),
     *                QMIN(NBU),QMAX(NBU) ,VMAG(NBU), VANG(NBU),
     *                DF(NBU)  ,DX(NBU)   ,A(NBU,NBU)
	COMMON/S17/GENBUS( NG ),V( NBU )

	INTEGER       BUSTYPE  , GENBUS
	COMPLEX       Y        , V
	REAL          XT   

	DO  I=1,NGEN
	WRITE(IPQLFO,950)( Y( GENBUS(I),GENBUS(J) ),J=1,NGEN)
	ENDDO
	WRITE(IPQLFO,960) XT
950     FORMAT(1X, 15(2F14.6, 2X) )
960     FORMAT(1X, F14.6 ) 
	RETURN
	END
