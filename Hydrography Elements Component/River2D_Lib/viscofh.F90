!==============================================================================|
!   Calculate Advection and Horizontal Diffusion Terms for Temperature         |
!==============================================================================|

   SUBROUTINE VISCOF_H               

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
   USE MOD_WD 
   USE MOD_SPHERICAL
   IMPLICIT NONE
   REAL(SP), DIMENSION(M)            :: PUPX,PUPY,PVPX,PVPY  
   REAL(SP), DIMENSION(M)            :: VISCOFF
   REAL(SP) :: X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2
   INTEGER  :: I,I1,IA,IB,J,J1,J2,K,JTMP
   REAL(DP) :: TXPI,TYPI
   REAL(DP) :: XTMP1,XTMP
   REAL(DP) :: X1_DP,Y1_DP,X2_DP,Y2_DP
   REAL(DP) :: X11_TMP,Y11_TMP,X33_TMP,Y33_TMP

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!

   DO K=1,KBM1
     DO I=1,M
       PUPX(I)=0.0_SP
       PUPY(I)=0.0_SP
       PVPX(I)=0.0_SP
       PVPY(I)=0.0_SP
       J=1
       I1=NBVE(I,J)
       JTMP=NBVT(I,J)
       J1=JTMP+1-(JTMP+1)/4*3
       J2=JTMP+2-(JTMP+2)/4*3
       X11=0.5_SP*(VX(I)+VX(NV(I1,J1)))
       Y11=0.5_SP*(VY(I)+VY(NV(I1,J1)))
       X22=XC(I1)
       Y22=YC(I1)
       X33=0.5_SP*(VX(I)+VX(NV(I1,J2)))
       Y33=0.5_SP*(VY(I)+VY(NV(I1,J2)))

      if (SPHERICAL==1)then
       X1_DP=VX(I)
       Y1_DP=VY(I)
       X2_DP=VX(NV(I1,J1))
       Y2_DP=VY(NV(I1,J1))
       CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X11_TMP,Y11_TMP)
       X11=X11_TMP
       Y11=Y11_TMP
       X2_DP=VX(NV(I1,J2))
       Y2_DP=VY(NV(I1,J2))
       CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X33_TMP,Y33_TMP)
       X33=X33_TMP
       Y33=Y33_TMP

       XTMP  = X33*TPI-X11*TPI
       XTMP1 = X33-X11
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP
       END IF	 
       TXPI=XTMP*COS(DEG2RAD*VY(I))
       TYPI=(Y11-Y33)*TPI

       PUPX(I)=PUPX(I)+U(I1,K)*TYPI
       PUPY(I)=PUPY(I)+U(I1,K)*TXPI
       PVPX(I)=PVPX(I)+V(I1,K)*TYPI
       PVPY(I)=PVPY(I)+V(I1,K)*TXPI
      else
       PUPX(I)=PUPX(I)+U(I1,K)*(Y11-Y33)
       PUPY(I)=PUPY(I)+U(I1,K)*(X33-X11)
       PVPX(I)=PVPX(I)+V(I1,K)*(Y11-Y33)
       PVPY(I)=PVPY(I)+V(I1,K)*(X33-X11)
      endif

       IF(ISONB(I) /= 0) THEN
        if  (SPHERICAL==1) then
         XTMP  = X11*TPI-VX(I)*TPI
         XTMP1 = X11-VX(I)
         IF(XTMP1 >  180.0_SP)THEN
	   XTMP = -360.0_SP*TPI+XTMP
         ELSE IF(XTMP1 < -180.0_SP)THEN
	   XTMP =  360.0_SP*TPI+XTMP
         END IF  
         TXPI=XTMP*COS(DEG2RAD*VY(I))
         TYPI=(VY(I)-Y11)*TPI

         PUPX(I)=PUPX(I)+U(I1,K)*TYPI
         PUPY(I)=PUPY(I)+U(I1,K)*TXPI
         PVPX(I)=PVPX(I)+V(I1,K)*TYPI
         PVPY(I)=PVPY(I)+V(I1,K)*TXPI
        else
         PUPX(I)=PUPX(I)+U(I1,K)*(VY(I)-Y11)
         PUPY(I)=PUPY(I)+U(I1,K)*(X11-VX(I))
         PVPX(I)=PVPX(I)+V(I1,K)*(VY(I)-Y11)
         PVPY(I)=PVPY(I)+V(I1,K)*(X11-VX(I))
        endif
       END IF

       DO J=2,NTVE(I)-1
         I1=NBVE(I,J)
         JTMP=NBVT(I,J)
         J1=JTMP+1-(JTMP+1)/4*3
         J2=JTMP+2-(JTMP+2)/4*3
         X11=0.5_SP*(VX(I)+VX(NV(I1,J1)))
         Y11=0.5_SP*(VY(I)+VY(NV(I1,J1)))
         X22=XC(I1)
         Y22=YC(I1)
         X33=0.5_SP*(VX(I)+VX(NV(I1,J2)))
         Y33=0.5_SP*(VY(I)+VY(NV(I1,J2)))

        if  (SPHERICAL==1)then
         X1_DP=VX(I)
         Y1_DP=VY(I)
         X2_DP=VX(NV(I1,J1))
         Y2_DP=VY(NV(I1,J1))
         CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X11_TMP,Y11_TMP)
	 X11=X11_TMP
	 Y11=Y11_TMP
         X2_DP=VX(NV(I1,J2))
         Y2_DP=VY(NV(I1,J2))
         CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X33_TMP,Y33_TMP)
	 X33=X33_TMP
	 Y33=Y33_TMP

         XTMP  = X33*TPI-X11*TPI
         XTMP1 = X33-X11
         IF(XTMP1 >  180.0_SP)THEN
	   XTMP = -360.0_SP*TPI+XTMP
         ELSE IF(XTMP1 < -180.0_SP)THEN
	   XTMP =  360.0_SP*TPI+XTMP
	 END IF  
         TXPI=XTMP*COS(DEG2RAD*VY(I))
         TYPI=(Y11-Y33)*TPI
         PUPX(I)=PUPX(I)+U(I1,K)*TYPI
         PUPY(I)=PUPY(I)+U(I1,K)*TXPI
         PVPX(I)=PVPX(I)+V(I1,K)*TYPI
         PVPY(I)=PVPY(I)+V(I1,K)*TXPI
        else
         PUPX(I)=PUPX(I)+U(I1,K)*(Y11-Y33)
         PUPY(I)=PUPY(I)+U(I1,K)*(X33-X11)
         PVPX(I)=PVPX(I)+V(I1,K)*(Y11-Y33)
         PVPY(I)=PVPY(I)+V(I1,K)*(X33-X11)
        endif
       END DO
       J=NTVE(I)
       I1=NBVE(I,J)
       JTMP=NBVT(I,J)
       J1=JTMP+1-(JTMP+1)/4*3
       J2=JTMP+2-(JTMP+2)/4*3
       X11=0.5_SP*(VX(I)+VX(NV(I1,J1)))
       Y11=0.5_SP*(VY(I)+VY(NV(I1,J1)))
       X22=XC(I1)
       Y22=YC(I1)
       X33=0.5_SP*(VX(I)+VX(NV(I1,J2)))
       Y33=0.5_SP*(VY(I)+VY(NV(I1,J2)))

      if  (SPHERICAL==1) then
       X1_DP=VX(I)
       Y1_DP=VY(I)
       X2_DP=VX(NV(I1,J1))
       Y2_DP=VY(NV(I1,J1))
       CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X11_TMP,Y11_TMP)
       X11=X11_TMP
       Y11=Y11_TMP
       X2_DP=VX(NV(I1,J2))
       Y2_DP=VY(NV(I1,J2))
       CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,X33_TMP,Y33_TMP)
       X33=X33_TMP
       Y33=Y33_TMP

       XTMP  = X33*TPI-X11*TPI
       XTMP1 = X33-X11
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP
       END IF
       TXPI=XTMP*COS(DEG2RAD*VY(I))
       TYPI=(Y11-Y33)*TPI

       PUPX(I)=PUPX(I)+U(I1,K)*TYPI
       PUPY(I)=PUPY(I)+U(I1,K)*TXPI
       PVPX(I)=PVPX(I)+V(I1,K)*TYPI
       PVPY(I)=PVPY(I)+V(I1,K)*TXPI
      else
       PUPX(I)=PUPX(I)+U(I1,K)*(Y11-Y33)
       PUPY(I)=PUPY(I)+U(I1,K)*(X33-X11)
       PVPX(I)=PVPX(I)+V(I1,K)*(Y11-Y33)
       PVPY(I)=PVPY(I)+V(I1,K)*(X33-X11)
      endif

       IF(ISONB(I) /= 0) THEN
       if (SPHERICAL==1)then
         XTMP  = VX(I)*TPI-X11*TPI
         XTMP1 = VX(I)-X11
         IF(XTMP1 >  180.0_SP)THEN
	   XTMP = -360.0_SP*TPI+XTMP
         ELSE IF(XTMP1 < -180.0_SP)THEN
	   XTMP =  360.0_SP*TPI+XTMP
	 END IF  
         TXPI=XTMP*COS(DEG2RAD*VY(I))
         TYPI=(Y11-VY(I))*TPI
	 
         PUPX(I)=PUPX(I)+U(I1,K)*TYPI
         PUPY(I)=PUPY(I)+U(I1,K)*TXPI
         PVPX(I)=PVPX(I)+V(I1,K)*TYPI
         PVPY(I)=PVPY(I)+V(I1,K)*TXPI
        else
         PUPX(I)=PUPX(I)+U(I1,K)*(Y11-VY(I))
         PUPY(I)=PUPY(I)+U(I1,K)*(VX(I)-X11)
         PVPX(I)=PVPX(I)+V(I1,K)*(Y11-VY(I))
         PVPY(I)=PVPY(I)+V(I1,K)*(VX(I)-X11)
        endif
       END IF
       PUPX(I)=PUPX(I)/ART1(I)
       PUPY(I)=PUPY(I)/ART1(I)
       PVPX(I)=PVPX(I)/ART1(I)
       PVPY(I)=PVPY(I)/ART1(I)
       TMP1=PUPX(I)**2+PVPY(I)**2
       TMP2=0.5_SP*(PUPY(I)+PVPX(I))**2
       VISCOFF(I)=SQRT(TMP1+TMP2)*ART1(I)
       
       VISCOFH(I,K) = VISCOFF(I)

     END DO
   END DO  
    
   RETURN
   END SUBROUTINE VISCOF_H
!==============================================================================|
