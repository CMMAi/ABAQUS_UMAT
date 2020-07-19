C========================== User MATerial  =========================
C
C===========================UMAT================================C
C                                                               C
C     User MATerial for power law creep                         C
C                                                               C
C     Programmer: Chi-Hua Yu PhD                                C
C                 National Taiwan University                    C
C     Version   : 2.0                                           C
C     Est. date : 01/20/2011                                    C
C     Last upd. : 07/10/2013                                    C
C     Documents : UMAT programming documents.docx               C
C     Descrp.   : This UMAT subroutine for ABAQUS is based on   C
C                 Dr. Gao's UMAT for BMG.                       C
C                                                               C
C===============================================================C

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 MATERL
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),
     4 DFGRD0(3,3),DFGRD1(3,3)
C -----------------------------------------------------------
C	Please refer to ABAQUS manual for the meaning of these arrays
C -----------------------------------------------------------
C
C      DIMENSION EELAS(6),EPLAS(6),FLOW(6),FINCR(2)
	DIMENSION SSDEVN(NTENS),SSDEVNP(NTENS),SSDEVNT(NTENS),DSNDEV(NTENS)
C     Varialbes table
C 
C     SSDEVN(NTENS)     : deviatoric stress in previous time step t
C     SSDEVNP(NTENS)    : deviatoric stress in present time step t+dt
C     SSDEVNT(NTENS)    : deviatoric stress preidctor
C     DSNDEV(NTENS)     : deviatoric strain increment
C     
C     DSMISES           : Mises stress increment
C     SMISESN           : Mises stress in previous time step t
C     SMISEST           : Mises stress predictor
C     SMISESNP          : Mises stress in present time step t+dt
C     SSTRACEN          : stress invariant I1
C     DSNMIS            : J2 invariant of DSNDEV
C     DECREQ            : effective strian increment
C     DECRNEW           : effective strain increment after NR iteration
C     FINCR             : -F(x)/F'(x) after NR iteration 
      PARAMETER (ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,SIX=6.0D0)
      DATA NEWTON,TOLER/1000000,1.D-10/
C
C     PROPS(1) - E
C	PROPS(2) - NU
C     PROPS(3) - SIMGA0
C     PROPS(4) - REFFERENCE STRAIN RATE
C     PROPS(5) - n
C
C	STATEV(1) - DECREQ
C     STATEV(2) - FLAG(0 FOR INITIAL, OTHERS IS 1)
C -----------------------------------------------------------
C
      IF (NDI.NE.3) THEN
         WRITE(6,1)
1       FORMAT(//,30X,'***ERROR - THIS UMAT MAY ONLY BE USED FOR ',
     1          'ELEMENTS WITH THREE DIRECT STRESS COMPONENTS')
C	which means that this code is applicable for CPE4, CAX4, C3D8, etc.
      ENDIF
C
C     INITIAL VALUES
C
	IF (STATEV(2).EQ.0) THEN
C        initial guess of effective creep strain rate
	   STATEV(1)=1.D-11
	   STATEV(2)=1

	ENDIF
C
	EMOD=PROPS(1)
	ENU=PROPS(2)
      IF(ENU.GT.0.4999.AND.ENU.LT.0.5001) ENU=0.499
	EBULK3=EMOD/(ONE-TWO*ENU)
      EG2=EMOD/(ONE+ENU)
      EG=EG2/TWO
      ELAM=(EBULK3-EG2)/THREE
      SIGMA0=PROPS(3)
      DE0=PROPS(4)
      RN=PROPS(5)
C
      DECREQ = STATEV(1)
C
C	STEP (i) -- compute deviatoric stress
C
	SMISESN=(STRESS(1)-STRESS(2))*(STRESS(1)-STRESS(2))+
	1	(STRESS(2)-STRESS(3))*(STRESS(2)-STRESS(3))+
	1	(STRESS(3)-STRESS(1))*(STRESS(3)-STRESS(1))
	DO 10 K1=NDI+1,NTENS
		SMISESN=SMISESN+SIX*STRESS(K1)*STRESS(K1)
 10	CONTINUE
	SMISESN=SQRT(SMISESN/TWO)
	SSTRACEN=STRESS(1)+STRESS(2)+STRESS(3)
	DO 20 K1=1,NDI
		SSDEVN(K1)=STRESS(K1)-SSTRACEN/THREE
 20	CONTINUE
	DO 30 K1=NDI+1,NTENS
		SSDEVN(K1)=STRESS(K1)
 30	CONTINUE
C	
C	STEP (ii) -- compute deviatoric strain
C
	DO 40 K1=1,NDI
		DSNDEV(K1)=DSTRAN(K1)-(DSTRAN(1)+DSTRAN(2)+DSTRAN(3))/THREE
 40	CONTINUE
      DO 50 K1=NDI+1,NTENS
		DSNDEV(K1)=DSTRAN(K1)
 50	CONTINUE
C	
C	J2 INVARIANT OF DSNDEV
C
	DSNMIS=(DSNDEV(1)-DSNDEV(2))*(DSNDEV(1)-DSNDEV(2))+
	1	(DSNDEV(2)-DSNDEV(3))*(DSNDEV(2)-DSNDEV(3))+
	1	(DSNDEV(3)-DSNDEV(1))*(DSNDEV(3)-DSNDEV(1))
	DO 55 K1=NDI+1,NTENS
		DSNMIS=DSNMIS+SIX*DSNDEV(K1)*DSNDEV(K1)
 55	CONTINUE
C
C	IF DSNMIS=0, WE CANNOT PROCEED, E.G., BETASTAR=0/0 IN STEP (vi).
C	JUMP TO THE	ELASTIC DEFORMAITON CASE.
C
	IF (DSNMIS.EQ.0.D0) GOTO 200
C
C	STEP (iii) -- compute trial stress, Eq. (6)
C
	DO 60 K1=1,NDI
		SSDEVNT(K1)=SSDEVN(K1)+EG2*DSNDEV(K1)
 60	CONTINUE
	DO 70 K1=NDI+1,NTENS
		SSDEVNT(K1)=SSDEVN(K1)+EG*DSNDEV(K1)
 70	CONTINUE
	SMISEST=(SSDEVNT(1)-SSDEVNT(2))*(SSDEVNT(1)-SSDEVNT(2))+
	1	(SSDEVNT(2)-SSDEVNT(3))*(SSDEVNT(2)-SSDEVNT(3))+
	1	(SSDEVNT(3)-SSDEVNT(1))*(SSDEVNT(3)-SSDEVNT(1))
	DO 80 K1=NDI+1,NTENS
		SMISEST=SMISEST+SIX*SSDEVNT(K1)*SSDEVNT(K1)
 80	CONTINUE
	SMISEST=SQRT(SMISEST/TWO)
C
C	STEP (iv) -- compute and Mises-stress increments, Appendix A
C
C      WRITE(*,*) "STEP IV"
	CALL DERIVATIVES(DECREQ,DTIME,RN,EMOD,ENU,
     1 SIGMA0,DE0,SMISEST,FINCR)
C
      DECRNEW = DECREQ+FINCR
C
	DO 90 KEWTON=1,NEWTON
C		NEED TO MODIFY THE FOLLOWING IF CONDITION, AT LEAST NORMALIZE 
	 IF ((ABS(DECRNEW-DECREQ).LT.TOLER))
	1    GOTO 100
	   DECREQ=DECRNEW
	   CALL DERIVATIVES(DECREQ,DTIME,RN,EMOD,ENU,
     1 SIGMA0,DE0,SMISEST,FINCR)
        DECRNEW = DECREQ+FINCR
C         DECRNEW = FINCR
C        WRITE(*,*) "DECREQ = ",DECREQ
C        WRITE(*,*) "FINCR = ",FINCR
 90	CONTINUE
	WRITE(6,2) NEWTON
 2    FORMAT(//,30X,'***WARNING - PLASTICITY ALGORITHM DID NOT ',
     1        'CONVERGE AFTER ',I10,' ITERATIONS')
	PNEWDT=0.3D-2
 100	CONTINUE
C
C	UPDATE STATE VARIABLES
C
	STATEV(1)=DECRNEW
C	WRITE(*,*) "EFFECTIVE STRAIN RATE = ", DECRNEW

C
C	STEP (v) -- compute Mises-strain increment, Eq. (11)
C

	BETASTAR=1.D0-THREE*EMOD/TWO/(1.D0+ENU)/SMISEST*DECRNEW

	SMISESNP=BETASTAR*SMISEST

C
C	STEP (vi) -- compute Eq. (12)
C	
	DO 110 K1=1,NDI
		SSDEVNP(K1)=BETASTAR*SSDEVNT(K1)
 110	CONTINUE
	DO 120 K1=NDI+1,NTENS
		SSDEVNP(K1)=BETASTAR*SSDEVNT(K1)
 120	CONTINUE
C
C	STEP (vii) -- compute the updated stress tensor, Eqs. (4) and (5)
C
	SSTRACENP=SSTRACEN+EBULK3*(DSTRAN(1)+DSTRAN(2)+DSTRAN(3))
	DO 130 K1=1,NDI
		STRESS(K1)=SSDEVNP(K1)+SSTRACENP/THREE
 130	CONTINUE
	DO 140 K1=NDI+1,NTENS
		STRESS(K1)=SSDEVNP(K1)
 140	CONTINUE
C
C	STEP (viii) -- compute psi in Eq. (B7)
C

C
C	STEP (ix) -- compute the elastic-plastic consistent tangent
C
	C1=THREE*EMOD/TWO/(1.D0+ENU)
	C2=C1/SMISESNP
	C3=C2+1./RN/DECREQ
	GAMMA=C3*BETASTAR
	TEMP1=DECREQ-1./GAMMA
	TEMP2=TEMP1/SMISESNP
	CONS=TEMP2*BETASTAR/SMISESNP/SMISESNP
	C1=C1**2*CONS
	DO 150 K1=1,NDI
		DDSDDE(K1,K1)=TWO*EMOD/THREE/(1.D0+ENU)*BETASTAR+
	1		C1*SSDEVNP(K1)*SSDEVNP(K1)+EBULK3/THREE
 150	CONTINUE
	DDSDDE(1,2)=-EG2/THREE*BETASTAR+EBULK3/THREE+
	1	C1*SSDEVNP(1)*SSDEVNP(2)
	DDSDDE(2,1)=DDSDDE(1,2)
	DDSDDE(1,3)=-EG2/THREE*BETASTAR+EBULK3/THREE+
	1	C1*SSDEVNP(1)*SSDEVNP(3)
	DDSDDE(3,1)=DDSDDE(1,3)
	DDSDDE(2,3)=-EG2/THREE*BETASTAR+EBULK3/THREE+
	1	C1*SSDEVNP(2)*SSDEVNP(3)
	DDSDDE(3,2)=DDSDDE(2,3)
	DO 160 K1=1,NDI
		DO 170 K2=NDI+1,NTENS
			DDSDDE(K1,K2)=C1*SSDEVNP(K1)*SSDEVNP(K2)
			DDSDDE(K2,K1)=DDSDDE(K1,K2)
 170		CONTINUE
 160	CONTINUE
	DO 180 K1=NDI+1,NTENS
		DDSDDE(K1,K1)=EG*BETASTAR+C1*SSDEVNP(K1)*SSDEVNP(K1)
 180	CONTINUE
	IF (NTENS.EQ.6) THEN
		DDSDDE(4,5)=C1*SSDEVNP(4)*SSDEVNP(5)
		DDSDDE(5,4)=DDSDDE(4,5)
		DDSDDE(4,6)=C1*SSDEVNP(4)*SSDEVNP(6)
		DDSDDE(6,4)=DDSDDE(4,6)
		DDSDDE(5,6)=C1*SSDEVNP(5)*SSDEVNP(6)
		DDSDDE(6,5)=DDSDDE(5,6)
	ELSEIF (NTENS.EQ.4) THEN
C		DDSDDE(4,4) HAS BEEN GIVEN
	ELSE
		WRITE(6,3) 
 3		FORMAT(//,30X,'***ERROR - THIS UMAT MAY ONLY BE USED FOR ',
     1          'ELEMENTS WITH ONE OR THREE SHEAR COMPONENTS')
	ENDIF
C
	GOTO 300
C

 200	CONTINUE
C	
C	ELASTIC STIFFNESS MATRIX
C

      DO 210 K1=1,NDI
        DO 220 K2=1,NDI
           DDSDDE(K2,K1)=ELAM
 220    CONTINUE
        DDSDDE(K1,K1)=EG2+ELAM
 210  CONTINUE
      DO 230 K1=NDI+1,NTENS
        DDSDDE(K1,K1)=EG
 230  CONTINUE

 300	CONTINUE

C
      RETURN
      END
C
C	Refer to Appendix A for the following subroutine
C
	SUBROUTINE DERIVATIVES(DECREQ,DTIME,RN,EMOD,ENU,
     1 SIGMA0,DE0,SMISEST,FINCR)
C	
	INCLUDE 'ABA_PARAM.INC'
C     
	IF(DECREQ .LT. 0.D0) DECREQ = 1.D-20 
	TEMP1=1-3*EMOD*DECREQ/2/(1+ENU)/SMISEST
	TEMP2=DECREQ/DTIME/DE0
C
C      WRITE(*,*) DECREQ
	F1=SMISEST*TEMP1/SIGMA0-TEMP2**(1./RN)
	C2=-3*EMOD/2/(1+ENU)/SIGMA0
	F2=C2-TEMP2**(1./RN)/RN/DECREQ
C     
      FINCR=-F1/F2
C
      RETURN
      END
