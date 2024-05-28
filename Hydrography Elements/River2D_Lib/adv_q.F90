!==============================================================================|
!   Calculate the Turbulent Kinetic Energy and Mixing Length Based  on         |
!   The Mellor-Yamada Level 2.5 Turbulent Closure Model                        |
!==============================================================================|

   SUBROUTINE ADV_Q(Q,QF)               

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_WD
   USE MOD_SPHERICAL
   USE MOD_NORTHPOLE
   IMPLICIT NONE
   REAL(SP), DIMENSION(0:MT,KB)     :: Q,QF,XFLUX
   REAL(SP), DIMENSION(M)           :: PUPX,PUPY,PVPX,PVPY  
   REAL(SP), DIMENSION(M)           :: PQPX,PQPY,PQPXD,PQPYD,VISCOFF
   REAL(SP), DIMENSION(3*(NT))      :: DTIJ 
   REAL(SP), DIMENSION(3*(NT),KBM1) :: UVN
   REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
   REAL(SP) :: Q1MIN, Q1MAX, Q2MIN, Q2MAX
   REAL(DP) :: TY,TXPI,TYPI
   REAL(DP) :: XTMP1,XTMP
   REAL(DP) :: X1_DP,Y1_DP,X2_DP,Y2_DP,XII,YII
   REAL(DP) :: X11_TMP,Y11_TMP,X33_TMP,Y33_TMP
   REAL(DP) :: VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP
   REAL(DP) :: TXPI_TMP,TYPI_TMP
   REAL(SP) :: QMEAN1
   REAL(SP), DIMENSION(0:NT,KB)    :: UQ,VQ
!------------------------------------------------------------------------------!

   QMEAN1 = 1.E-8
   
   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF
     
!
!--Initialize Fluxes-----------------------------------------------------------!
!
   QF    = 0.0_SP
   XFLUX = 0.0_SP
   
   UQ = 0.0_SP
   VQ = 0.0_SP
   UVN = 0.0_SP
   
   DO K=2,KBM1
     DO I=1,NT
       UQ(I,K) = (U(I,K)*DZ(K-1)+U(I,K-1)*DZ(K))/(DZ(K)+DZ(K-1))
       VQ(I,K) = (V(I,K)*DZ(K-1)+V(I,K-1)*DZ(K))/(DZ(K)+DZ(K-1))
     END DO
   END DO     

!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
   DO I=1,NCV
     I1=NTRG(I)
     DTIJ(I)=DT1(I1)
     DO K=2,KBM1
       UVN(I,K) = VQ(I1,K)*DLTXE(I) - UQ(I1,K)*DLTYE(I)
     END DO
   END DO

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!

   DO K=2,KBM1
     PQPX  = 0.0_SP 
     PQPY  = 0.0_SP 
     PQPXD = 0.0_SP 
     PQPYD = 0.0_SP
     DO I=1,M
       DO J=1,NTSN(I)-1
         I1=NBSN(I,J)
         I2=NBSN(I,J+1)
!         FFD=0.5_SP*(Q(I1,K)+Q(I2,K)-QMEAN1(I1,K)-QMEAN1(I2,K))
         FFD=0.5_SP*(Q(I1,K)+Q(I2,K)-QMEAN1-QMEAN1)
         FF1=0.5_SP*(Q(I1,K)+Q(I2,K))
        if (SPHERICAL==1)then
         XTMP  = VX(I2)*TPI-VX(I1)*TPI
	 XTMP1 = VX(I2)-VX(I1)
	 IF(XTMP1 >  180.0_SP)THEN
	   XTMP = -360.0_SP*TPI+XTMP
	 ELSE IF(XTMP1 < -180.0_SP)THEN
	   XTMP =  360.0_SP*TPI+XTMP
	 END IF  
         TXPI=XTMP*COS(DEG2RAD*VY(I))
         TYPI=(VY(I1)-VY(I2))*TPI
    if (NORTHPOLE==1)then
         IF(NODE_NORTHAREA(I) == 1)THEN
           VX1_TMP = REARTH * COS(VY(I1)*PI/180.0_SP) * COS(VX(I1)*PI/180.0_SP) &
                     * 2._SP /(1._SP+sin(VY(I1)*PI/180.0_SP))
           VY1_TMP = REARTH * COS(VY(I1)*PI/180.0_SP) * SIN(VX(I1)*PI/180.0_SP) &
                     * 2._SP /(1._SP+sin(VY(I1)*PI/180.0_SP))

           VX2_TMP = REARTH * COS(VY(I2)*PI/180.0_SP) * COS(VX(I2)*PI/180.0_SP) &
                     * 2._SP /(1._SP+sin(VY(I2)*PI/180.0_SP))
           VY2_TMP = REARTH * COS(VY(I2)*PI/180.0_SP) * SIN(VX(I2)*PI/180.0_SP) &
                     * 2._SP /(1._SP+sin(VY(I2)*PI/180.0_SP))

           TXPI = (VX2_TMP-VX1_TMP)/(2._SP /(1._SP+sin(VY(I)*PI/180.0_SP)))
           TYPI = (VY1_TMP-VY2_TMP)/(2._SP /(1._SP+sin(VY(I)*PI/180.0_SP)))
  	   IF(I /= NODE_NORTHPOLE)THEN
	     TXPI_TMP = TYPI*COS(VX(I)*PI/180.)-TXPI*SIN(VX(I)*PI/180.)
	     TYPI_TMP = TXPI*COS(VX(I)*PI/180.)+TYPI*SIN(VX(I)*PI/180.)
	     TYPI_TMP = -TYPI_TMP
	    
	     TXPI = TXPI_TMP
	     TYPI = TYPI_TMP
	   END IF  
	 END IF 
    endif	 
         PQPX(I)=PQPX(I)+FF1*TYPI
         PQPY(I)=PQPY(I)+FF1*TXPI
         PQPXD(I)=PQPXD(I)+FFD*TYPI
         PQPYD(I)=PQPYD(I)+FFD*TXPI
        else
         PQPX(I)=PQPX(I)+FF1*(VY(I1)-VY(I2))
         PQPY(I)=PQPY(I)+FF1*(VX(I2)-VX(I1))
         PQPXD(I)=PQPXD(I)+FFD*(VY(I1)-VY(I2))
         PQPYD(I)=PQPYD(I)+FFD*(VX(I2)-VX(I1))
        endif
       END DO
       PQPX(I)=PQPX(I)/ART2(I)
       PQPY(I)=PQPY(I)/ART2(I)
       PQPXD(I)=PQPXD(I)/ART2(I)
       PQPYD(I)=PQPYD(I)/ART2(I)
     END DO
          
     DO I=1,M
!       PUPX(I)=0.0_SP
!       PUPY(I)=0.0_SP
!       PVPX(I)=0.0_SP
!       PVPY(I)=0.0_SP
!       J=1
!       I1=NBVE(I,J)
!       JTMP=NBVT(I,J)
!       J1=JTMP+1-(JTMP+1)/4*3
!       J2=JTMP+2-(JTMP+2)/4*3
!       X11=0.5_SP*(VX(I)+VX(NV(I1,J1)))
!       Y11=0.5_SP*(VY(I)+VY(NV(I1,J1)))
!       X22=XC(I1)
!       Y22=YC(I1)
!       X33=0.5_SP*(VX(I)+VX(NV(I1,J2)))
!       Y33=0.5_SP*(VY(I)+VY(NV(I1,J2)))

!#      if defined (SPHERICAL)
!       X1_DP=VX(I)
!       Y1_DP=VY(I)
!       X2_DP=VX(NV(I1,J1))
!       Y2_DP=VY(NV(I1,J1))
!       CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X11_TMP,Y11_TMP)
!       X11=X11_TMP
!       Y11=Y11_TMP
!      X2_DP=VX(NV(I1,J2))
!       Y2_DP=VY(NV(I1,J2))
!       CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X33_TMP,Y33_TMP)
!       X33=X33_TMP
!       Y33=Y33_TMP

!       XTMP  = X33*TPI-X11*TPI
!       XTMP1 = X33-X11
!       IF(XTMP1 >  180.0_SP)THEN
!        XTMP = -360.0_SP*TPI+XTMP
!       ELSE IF(XTMP1 < -180.0_SP)THEN
!         XTMP =  360.0_SP*TPI+XTMP
!       END IF	 
!       TXPI=XTMP*COS(DEG2RAD*VY(I))
!       TYPI=(Y11-Y33)*TPI
       
!       PUPX(I)=PUPX(I)+UQ(I1,K)*TYPI
!       PUPY(I)=PUPY(I)+UQ(I1,K)*TXPI
!      PVPX(I)=PVPX(I)+VQ(I1,K)*TYPI
!      PVPY(I)=PVPY(I)+VQ(I1,K)*TXPI
!#      else
!       PUPX(I)=PUPX(I)+UQ(I1,K)*(Y11-Y33)
!       PUPY(I)=PUPY(I)+UQ(I1,K)*(X33-X11)
!       PVPX(I)=PVPX(I)+VQ(I1,K)*(Y11-Y33)
!       PVPY(I)=PVPY(I)+VQ(I1,K)*(X33-X11)
!#      endif

!       IF(ISONB(I) /= 0) THEN
!#        if defined (SPHERICAL)
!         XTMP  = X11*TPI-VX(I)*TPI
!         XTMP1 = X11-VX(I)
!         IF(XTMP1 >  180.0_SP)THEN
!	   XTMP = -360.0_SP*TPI+XTMP
!         ELSE IF(XTMP1 < -180.0_SP)THEN
!	   XTMP =  360.0_SP*TPI+XTMP
!	 END IF  
!         TXPI=XTMP*COS(DEG2RAD*VY(I))
!         TYPI=(VY(I)-Y11)*TPI
	 
!         PUPX(I)=PUPX(I)+UQ(I1,K)*TYPI
!        PUPY(I)=PUPY(I)+UQ(I1,K)*TXPI
!         PVPX(I)=PVPX(I)+VQ(I1,K)*TYPI
!         PVPY(I)=PVPY(I)+VQ(I1,K)*TXPI
!#        else
!         PUPX(I)=PUPX(I)+UQ(I1,K)*(VY(I)-Y11)
!         PUPY(I)=PUPY(I)+UQ(I1,K)*(X11-VX(I))
!         PVPX(I)=PVPX(I)+VQ(I1,K)*(VY(I)-Y11)
!         PVPY(I)=PVPY(I)+VQ(I1,K)*(X11-VX(I))
!#        endif
!       END IF

!       DO J=2,NTVE(I)-1
!         I1=NBVE(I,J)
!         JTMP=NBVT(I,J)
!         J1=JTMP+1-(JTMP+1)/4*3
!         J2=JTMP+2-(JTMP+2)/4*3
!         X11=0.5_SP*(VX(I)+VX(NV(I1,J1)))
!         Y11=0.5_SP*(VY(I)+VY(NV(I1,J1)))
!         X22=XC(I1)
!         Y22=YC(I1)
!         X33=0.5_SP*(VX(I)+VX(NV(I1,J2)))
!         Y33=0.5_SP*(VY(I)+VY(NV(I1,J2)))

!#        if defined (SPHERICAL)
!         X1_DP=VX(I)
!         Y1_DP=VY(I)
!         X2_DP=VX(NV(I1,J1))
!         Y2_DP=VY(NV(I1,J1))
!         CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X11_TMP,Y11_TMP)
!	 X11=X11_TMP
!	 Y11=Y11_TMP
!         X2_DP=VX(NV(I1,J2))
!         Y2_DP=VY(NV(I1,J2))
!         CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X33_TMP,Y33_TMP)
!	 X33=X33_TMP
!	 Y33=Y33_TMP

!         XTMP  = X33*TPI-X11*TPI
!         XTMP1 = X33-X11
!         IF(XTMP1 >  180.0_SP)THEN
!	   XTMP = -360.0_SP*TPI+XTMP
!         ELSE IF(XTMP1 < -180.0_SP)THEN
!	   XTMP =  360.0_SP*TPI+XTMP
!	 END IF  
!         TXPI=XTMP*COS(DEG2RAD*VY(I))
!         TYPI=(Y11-Y33)*TPI

!         PUPX(I)=PUPX(I)+UQ(I1,K)*TYPI
!         PUPY(I)=PUPY(I)+UQ(I1,K)*TXPI
!         PVPX(I)=PVPX(I)+VQ(I1,K)*TYPI
!         PVPY(I)=PVPY(I)+VQ(I1,K)*TXPI
!#        else
!         PUPX(I)=PUPX(I)+UQ(I1,K)*(Y11-Y33)
!         PUPY(I)=PUPY(I)+UQ(I1,K)*(X33-X11)
!         PVPX(I)=PVPX(I)+VQ(I1,K)*(Y11-Y33)
!         PVPY(I)=PVPY(I)+VQ(I1,K)*(X33-X11)
!#        endif
!       END DO
!       J=NTVE(I)
!       I1=NBVE(I,J)
!       JTMP=NBVT(I,J)
!       J1=JTMP+1-(JTMP+1)/4*3
!       J2=JTMP+2-(JTMP+2)/4*3
!       X11=0.5_SP*(VX(I)+VX(NV(I1,J1)))
!       Y11=0.5_SP*(VY(I)+VY(NV(I1,J1)))
!       X22=XC(I1)
!       Y22=YC(I1)
!       X33=0.5_SP*(VX(I)+VX(NV(I1,J2)))
!       Y33=0.5_SP*(VY(I)+VY(NV(I1,J2)))

!#      if defined (SPHERICAL)
!       X1_DP=VX(I)
!       Y1_DP=VY(I)
!       X2_DP=VX(NV(I1,J1))
!       Y2_DP=VY(NV(I1,J1))
!       CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X11_TMP,Y11_TMP)
!       X11=X11_TMP
!       Y11=Y11_TMP
!       X2_DP=VX(NV(I1,J2))
!       Y2_DP=VY(NV(I1,J2))
!       CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X33_TMP,Y33_TMP)
!       X33=X33_TMP
!       Y33=Y33_TMP

!       XTMP  = X33*TPI-X11*TPI
!       XTMP1 = X33-X11
!       IF(XTMP1 >  180.0_SP)THEN
!         XTMP = -360.0_SP*TPI+XTMP
!       ELSE IF(XTMP1 < -180.0_SP)THEN
!         XTMP =  360.0_SP*TPI+XTMP
!       END IF	 
!       TXPI=XTMP*COS(DEG2RAD*VY(I))
!       TYPI=(Y11-Y33)*TPI

!       PUPX(I)=PUPX(I)+UQ(I1,K)*TYPI
!       PUPY(I)=PUPY(I)+UQ(I1,K)*TXPI
!       PVPX(I)=PVPX(I)+VQ(I1,K)*TYPI
!       PVPY(I)=PVPY(I)+VQ(I1,K)*TXPI
!#      else
!       PUPX(I)=PUPX(I)+UQ(I1,K)*(Y11-Y33)
!       PUPY(I)=PUPY(I)+UQ(I1,K)*(X33-X11)
!       PVPX(I)=PVPX(I)+VQ(I1,K)*(Y11-Y33)
!       PVPY(I)=PVPY(I)+VQ(I1,K)*(X33-X11)
!#      endif

!       IF(ISONB(I) /= 0) THEN
!#      if defined (SPHERICAL)
!         XTMP  = VX(I)*TPI-X11*TPI
!         XTMP1 = VX(I)-X11
!         IF(XTMP1 >  180.0_SP)THEN
!	   XTMP = -360.0_SP*TPI+XTMP
!         ELSE IF(XTMP1 < -180.0_SP)THEN
!	   XTMP =  360.0_SP*TPI+XTMP
!	 END IF  
!         TXPI=XTMP*COS(DEG2RAD*VY(I))
!         TYPI=(Y11-VY(I))*TPI

!         PUPX(I)=PUPX(I)+UQ(I1,K)*TYPI
!         PUPY(I)=PUPY(I)+UQ(I1,K)*TXPI
!         PVPX(I)=PVPX(I)+VQ(I1,K)*TYPI
!         PVPY(I)=PVPY(I)+VQ(I1,K)*TXPI
!#        else
!         PUPX(I)=PUPX(I)+UQ(I1,K)*(Y11-VY(I))
!         PUPY(I)=PUPY(I)+UQ(I1,K)*(VX(I)-X11)
!         PVPX(I)=PVPX(I)+VQ(I1,K)*(Y11-VY(I))
!         PVPY(I)=PVPY(I)+VQ(I1,K)*(VX(I)-X11)
!#        endif
!       END IF
!       PUPX(I)=PUPX(I)/ART1(I)
!       PUPY(I)=PUPY(I)/ART1(I)
!       PVPX(I)=PVPX(I)/ART1(I)
!       PVPY(I)=PVPY(I)/ART1(I)
!       TMP1=PUPX(I)**2+PVPY(I)**2
!       TMP2=0.5_SP*(PUPY(I)+PVPX(I))**2
!       VISCOFF(I)=SQRT(TMP1+TMP2)*ART1(I)
 
       VISCOFF(I) = (VISCOFH(I,K)*DZ(K-1)+VISCOFH(I,K-1)*DZ(K))/(DZ(K)+DZ(K-1))
       
     END DO
!     IF(K == KBM1) THEN
!       AH_BOTTOM(1:M) = HORCON*(FACT*VISCOFF(1:M) + FM1)
!     END IF


     DO I=1,NCV_I
       IA=NIEC(I,1)
       IB=NIEC(I,2)
       XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
       YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
      if (SPHERICAL==1)then
       X1_DP=XIJE(I,1)
       Y1_DP=YIJE(I,1)
       X2_DP=XIJE(I,2)
       Y2_DP=YIJE(I,2)
       CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,XII,YII)
       XI=XII		
       XTMP  = XI*TPI-VX(IA)*TPI
       XTMP1 = XI-VX(IA)
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP
       END IF	 

       DXA=XTMP*COS(DEG2RAD*VY(IA))    
       DYA=(YI-VY(IA))*TPI
       XTMP  = XI*TPI-VX(IB)*TPI
       XTMP1 = XI-VX(IB)
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP
       END IF	 

       DXB=XTMP*COS(DEG2RAD*VY(IB)) 
       DYB=(YI-VY(IB))*TPI
      else
       DXA=XI-VX(IA)
       DYA=YI-VY(IA)
       DXB=XI-VX(IB)
       DYB=YI-VY(IB)
      endif

       FIJ1=Q(IA,K)+DXA*PQPX(IA)+DYA*PQPY(IA)
       FIJ2=Q(IB,K)+DXB*PQPX(IB)+DYB*PQPY(IB)

       Q1MIN=MINVAL(Q(NBSN(IA,1:NTSN(IA)-1),K))
       Q1MIN=MIN(Q1MIN, Q(IA,K))
       Q1MAX=MAXVAL(Q(NBSN(IA,1:NTSN(IA)-1),K))
       Q1MAX=MAX(Q1MAX, Q(IA,K))
       Q2MIN=MINVAL(Q(NBSN(IB,1:NTSN(IB)-1),K))
       Q2MIN=MIN(Q2MIN, Q(IB,K))
       Q2MAX=MAXVAL(Q(NBSN(IB,1:NTSN(IB)-1),K))
       Q2MAX=MAX(Q2MAX, Q(IB,K))
       IF(FIJ1 < Q1MIN) FIJ1=Q1MIN
       IF(FIJ1 > Q1MAX) FIJ1=Q1MAX
       IF(FIJ2 < Q2MIN) FIJ2=Q2MIN
       IF(FIJ2 > Q2MAX) FIJ2=Q2MAX
    
       UN=UVN(I,K)

!       VISCOF=HORCON*(FACT*0.5_SP*(VISCOFF(IA)+VISCOFF(IB))/HPRNU + FM1)
       VISCOF=HORCON*(FACT*0.5_SP*(VISCOFF(IA)+VISCOFF(IB)) + FM1)/HPRNU

       TXX=0.5_SP*(PQPXD(IA)+PQPXD(IB))*VISCOF
       TYY=0.5_SP*(PQPYD(IA)+PQPYD(IB))*VISCOF

       FXX=-DTIJ(I)*TXX*DLTYE(I)
       FYY= DTIJ(I)*TYY*DLTXE(I)

       EXFLUX=-UN*DTIJ(I)* &
          ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP+FXX+FYY

       XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX
       XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX

     END DO

    if  (SPHERICAL==1)then
    if  (NORTHPOLE==1)then
     CALL ADV_Q_XY(XFLUX,PQPX,PQPY,PQPXD,PQPYD,VISCOFF,Q,UQ,VQ,K)
    endif
    endif  

   END DO !!SIGMA LOOP

!
!-Accumulate Fluxes at Boundary Nodes
!

 
!--------------------------------------------------------------------
!   The central difference scheme in vertical advection
!--------------------------------------------------------------------
   DO K=2,KBM1
     DO I=1,M
       IF(ISWETN(I)*ISWETNT(I) == 1) THEN
         TEMP=WTS(I,K-1)*Q(I,K-1)-WTS(I,K+1)*Q(I,K+1)
         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)/(DZ(K-1)+DZ(K))
       END IF
     END DO
   END DO  !! SIGMA LOOP

!
!--Update Q or QL-------------------------------------------------------------!
!

   DO I=1,M
     IF(ISWETN(I)*ISWETNT(I) == 1 )THEN
     DO K=2,KBM1
       QF(I,K)=(Q(I,K)-XFLUX(I,K)/ART1(I)*(DTI/DT(I)))*(DT(I)/D(I))
     END DO
     ELSE
     DO K=2,KBM1
       QF(I,K)=Q(I,K)
     END DO
     END IF
   END DO

   RETURN
   END SUBROUTINE ADV_Q
!==============================================================================|
