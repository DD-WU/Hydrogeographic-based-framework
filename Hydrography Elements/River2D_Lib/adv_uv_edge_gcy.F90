!==============================================================================!

   SUBROUTINE ADV_UV_EDGE_GCY

!==============================================================================!
! this subroutine calculate advective, coriolis, pressure gradient, etc in     !
! x and y momentum equations except vertical diffusion terms for internal mode ! 
!==============================================================================!

   USE ALL_VARS
   USE MOD_SPHERICAL
   USE MOD_NORTHPOLE
   USE BCS
   USE MOD_WD
   USE MOD_EQUITIDE

   IMPLICIT NONE
   REAL(SP) :: XFLUX(0:NT,KB),YFLUX(0:NT,KB)
   REAL(SP) :: PSTX_TM(0:NT,KB),PSTY_TM(0:NT,KB)
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: FXX,FYY,XADV,YADV,TXXIJ,TYYIJ,TXYIJ,UN
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP,TPA,TPB
   REAL(SP) :: XIJA,YIJA,XIJB,YIJB,UIJ,VIJ,UIJ1,VIJ1,UIJ2,VIJ2
   REAL(SP) :: DIJ,ELIJ,TMPA,TMPB,TMP,XFLUXV,YFLUXV
   REAL(SP) :: FACT,FM1,EXFLUX,ISWETTMP
   INTEGER  :: I,IA,IB,J1,J2,K1,K2,K3,K4,K5,K6,K,II,J,I1,I2

!Q&C< used for new BC 06/22/04
   REAL(SP) :: UK1(KB),UK2(KB),UK3(KB),UK4(KB),UK5(KB),UK6(KB), &
               VK1(KB),VK2(KB),VK3(KB),VK4(KB),VK5(KB),VK6(KB)
!>Q&C 06/22/04
   REAL(DP) :: XTMP,XTMP1
   REAL(SP) :: U_TMP,V_TMP,UF_TMP,VF_TMP


!------------------------------------------------------------------------------!

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

!
!-----Initialize Flux Variables------------------------------------------------!
!
   XFLUX  = 0.0_SP
   YFLUX  = 0.0_SP
   PSTX_TM = 0.0_SP
   PSTY_TM = 0.0_SP

!
!-----Loop Over Edges and Accumulate Flux--------------------------------------!
!
   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)
     DIJ=0.5_SP*(DT(J1)+DT(J2))
     ELIJ=0.5_SP*(EGF(J1)+EGF(J2))
    if  (EQUI_TIDE==1) then
     ELIJ=ELIJ-0.5_SP*(EGF_EQI(J1)+EGF_EQI(J2))
    endif


     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
     K4=NBE(IB,1)
     K5=NBE(IB,2)
     K6=NBE(IB,3)
    if  (SPHERICAL==1) then
     XIJA=DLTXNE(I,1)
     YIJA=DLTYNE(I,1)
     XIJB=DLTXNE(I,2)
     YIJB=DLTYNE(I,2)
    else
     XIJA=XIJC(I)-XC(IA)
     YIJA=YIJC(I)-YC(IA)
     XIJB=XIJC(I)-XC(IB)
     YIJB=YIJC(I)-YC(IB)
    endif

!Q&C< used for new BC 07/26/04
     UK1 = U(K1,:)
     UK2 = U(K2,:)
     UK3 = U(K3,:)
     UK4 = U(K4,:)
     UK5 = U(K5,:)
     UK6 = U(K6,:)
     VK1 = V(K1,:)
     VK2 = V(K2,:)
     VK3 = V(K3,:)
     VK4 = V(K4,:)
     VK5 = V(K5,:)
     VK6 = V(K6,:)
	 if(IA*IB/=0) then
     IF(K1 == 0) CALL GHOSTUV3(IA,1,UK1,VK1)
     IF(K2 == 0) CALL GHOSTUV3(IA,2,UK2,VK2)
     IF(K3 == 0) CALL GHOSTUV3(IA,3,UK3,VK3)
     IF(K4 == 0) CALL GHOSTUV3(IB,1,UK4,VK4)
     IF(K5 == 0) CALL GHOSTUV3(IB,2,UK5,VK5)
     IF(K6 == 0) CALL GHOSTUV3(IB,3,UK6,VK6)
	 endif
!>Q&C 07/26/04

     DO K=1,KBM1
      IF(ISWETCT(IA)*ISWETC(IA) == 1 .OR. ISWETCT(IB)*ISWETC(IB) == 1)THEN

!Q&C< used for new BC 06/22/04
       COFA1=A1U(IA,1)*U(IA,K)+A1U(IA,2)*UK1(K)+A1U(IA,3)*UK2(K)+A1U(IA,4)*UK3(K)
       COFA2=A2U(IA,1)*U(IA,K)+A2U(IA,2)*UK1(K)+A2U(IA,3)*UK2(K)+A2U(IA,4)*UK3(K)
       COFA5=A1U(IA,1)*V(IA,K)+A1U(IA,2)*VK1(K)+A1U(IA,3)*VK2(K)+A1U(IA,4)*VK3(K)
       COFA6=A2U(IA,1)*V(IA,K)+A2U(IA,2)*VK1(K)+A2U(IA,3)*VK2(K)+A2U(IA,4)*VK3(K)
!       COFA1=A1U(IA,1)*U(IA,K)+A1U(IA,2)*U(K1,K)+A1U(IA,3)*U(K2,K)+A1U(IA,4)*U(K3,K)
!       COFA2=A2U(IA,1)*U(IA,K)+A2U(IA,2)*U(K1,K)+A2U(IA,3)*U(K2,K)+A2U(IA,4)*U(K3,K)
!       COFA5=A1U(IA,1)*V(IA,K)+A1U(IA,2)*V(K1,K)+A1U(IA,3)*V(K2,K)+A1U(IA,4)*V(K3,K)
!       COFA6=A2U(IA,1)*V(IA,K)+A2U(IA,2)*V(K1,K)+A2U(IA,3)*V(K2,K)+A2U(IA,4)*V(K3,K)
!Q&C> 06/22/04

       UIJ1=U(IA,K)+COFA1*XIJA+COFA2*YIJA
       VIJ1=V(IA,K)+COFA5*XIJA+COFA6*YIJA

!Q&C< used for new BC 06/22/04
       COFA3=A1U(IB,1)*U(IB,K)+A1U(IB,2)*UK4(K)+A1U(IB,3)*UK5(K)+A1U(IB,4)*UK6(K)
       COFA4=A2U(IB,1)*U(IB,K)+A2U(IB,2)*UK4(K)+A2U(IB,3)*UK5(K)+A2U(IB,4)*UK6(K)
       COFA7=A1U(IB,1)*V(IB,K)+A1U(IB,2)*VK4(K)+A1U(IB,3)*VK5(K)+A1U(IB,4)*VK6(K)
       COFA8=A2U(IB,1)*V(IB,K)+A2U(IB,2)*VK4(K)+A2U(IB,3)*VK5(K)+A2U(IB,4)*VK6(K)
!       COFA3=A1U(IB,1)*U(IB,K)+A1U(IB,2)*U(K4,K)+A1U(IB,3)*U(K5,K)+A1U(IB,4)*U(K6,K)
!       COFA4=A2U(IB,1)*U(IB,K)+A2U(IB,2)*U(K4,K)+A2U(IB,3)*U(K5,K)+A2U(IB,4)*U(K6,K)
!       COFA7=A1U(IB,1)*V(IB,K)+A1U(IB,2)*V(K4,K)+A1U(IB,3)*V(K5,K)+A1U(IB,4)*V(K6,K)
!       COFA8=A2U(IB,1)*V(IB,K)+A2U(IB,2)*V(K4,K)+A2U(IB,3)*V(K5,K)+A2U(IB,4)*V(K6,K)
!!>Q&C 06/22/04

       UIJ2=U(IB,K)+COFA3*XIJB+COFA4*YIJB
       VIJ2=V(IB,K)+COFA7*XIJB+COFA8*YIJB

       UIJ=0.5_SP*(UIJ1+UIJ2)
       VIJ=0.5_SP*(VIJ1+VIJ2)
       EXFLUX = DIJ*(-UIJ*DLTYC(I) + VIJ*DLTXC(I))

!
!-------ADD THE VISCOUS TERM & ADVECTION TERM---------------------------------!
!

       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)

!       VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON
       VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON/HPRNU

       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
       FYY=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))

       XADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*UIJ2+(1.0_SP+SIGN(1.0_SP,EXFLUX))*UIJ1)*0.5_SP
       YADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*VIJ2+(1.0_SP+SIGN(1.0_SP,EXFLUX))*VIJ1)*0.5_SP

       !!CALCULATE BOUNDARY FLUX AUGMENTERS
       TPA = FLOAT(1-ISBC(I))*EPOR(IA)
       TPB = FLOAT(1-ISBC(I))*EPOR(IB)


       !!ACCUMULATE ADVECTIVE + DIFFUSIVE + BAROTROPIC PRESSURE GRADIENT TERMS
       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+FXX*TPA
       YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+FYY*TPA
       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-FXX*TPB
       YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-FYY*TPB

    END IF
!JQI<for spherical coordinator and domain across 360^o latitude
    if  (SPHERICAL==1)then
        XTMP  = VX(J2)*TPI-VX(J1)*TPI
        XTMP1 = VX(J2)-VX(J1)
        IF(XTMP1 >  180.0_SP)THEN
	  XTMP = -360.0_SP*TPI+XTMP
        ELSE IF(XTMP1 < -180.0_SP)THEN
	  XTMP =  360.0_SP*TPI+XTMP
	END IF
	  
        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV*DT1(IA)*ELIJ*DLTYC(I)
        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV*DT1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA)) 
        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV*DT1(IB)*ELIJ*DLTYC(I)
        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV*DT1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB)) 
    else
!JQI>
        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV*DT1(IA)*ELIJ*DLTYC(I)
        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV*DT1(IA)*ELIJ*DLTXC(I)
        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV*DT1(IB)*ELIJ*DLTYC(I)
        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV*DT1(IB)*ELIJ*DLTXC(I)
!JQI<
    endif     
!JQI>
     END DO
   END DO

      DO I=1,N
       ISWETTMP = ISWETCT(I)*ISWETC(I)
        DO K=1,KBM1
	 XFLUX(I,K)  = XFLUX(I,K)*ISWETTMP
	 YFLUX(I,K)  = YFLUX(I,K)*ISWETTMP
         PSTX_TM(I,K)= PSTX_TM(I,K)*ISWETTMP
         PSTY_TM(I,K)= PSTY_TM(I,K)*ISWETTMP
        END DO
       DO K=1,KBM1
        XFLUX(I,K)=XFLUX(I,K)+PSTX_TM(I,K)
        YFLUX(I,K)=YFLUX(I,K)+PSTY_TM(I,K)
       END DO
      END DO

!
!-------ADD VERTICAL CONVECTIVE FLUX, CORIOLIS TERM AND BAROCLINIC PG TERM----!
!
   DO I=1,N
     IF(ISWETCT(I)*ISWETC(I) == 1)THEN
     DO K=1,KBM1
       IF(K == 1) THEN
         XFLUXV=-W(I,K+1)*(U(I,K)*DZ(K+1)+U(I,K+1)*DZ(K))/&
                 (DZ(K)+DZ(K+1))
         YFLUXV=-W(I,K+1)*(V(I,K)*DZ(K+1)+V(I,K+1)*DZ(K))/&
                 (DZ(K)+DZ(K+1))
       ELSE IF(K == KBM1) THEN
         XFLUXV= W(I,K)*(U(I,K)*DZ(K-1)+U(I,K-1)*DZ(K))/&
                 (DZ(K)+DZ(K-1))
         YFLUXV= W(I,K)*(V(I,K)*DZ(K-1)+V(I,K-1)*DZ(K))/&
                 (DZ(K)+DZ(K-1))
       ELSE
         XFLUXV= W(I,K)*(U(I,K)*DZ(K-1)+U(I,K-1)*DZ(K))/&
                 (DZ(K)+DZ(K-1))-&
                 W(I,K+1)*(U(I,K)*DZ(K+1)+U(I,K+1)*DZ(K))/&
                 (DZ(K)+DZ(K+1))
         YFLUXV= W(I,K)*(V(I,K)*DZ(K-1)+V(I,K-1)*DZ(K))/&
                 (DZ(K)+DZ(K-1))-&
                 W(I,K+1)*(V(I,K)*DZ(K+1)+V(I,K+1)*DZ(K))/&
                 (DZ(K)+DZ(K+1))
       END IF
      if (SPHERICAL==1) then
       XFLUX(I,K)=XFLUX(I,K)+XFLUXV/DZ(K)*ART(I)&
                 +DRHOX(I,K)-COR(I)*V(I,K)*DT1(I)*ART(I)&
                 -U(I,K)*V(I,K)/REARTH*TAN(YC(I)*PI/180.0_SP)*DT1(I)*ART(I)&
                 +0.5_SP*U(I,K)*(W(I,K+1)+W(I,K))/REARTH*DT1(I)*ART(I)
       YFLUX(I,K)=YFLUX(I,K)+YFLUXV/DZ(K)*ART(I)&
                 +DRHOY(I,K)+COR(I)*U(I,K)*DT1(I)*ART(I)&
                 +U(I,K)*U(I,K)/REARTH*TAN(YC(I)*PI/180.0_SP)*DT1(I)*ART(I)&
                 +0.5_SP*V(I,K)*(W(I,K+1)+W(I,K))/REARTH*DT1(I)*ART(I)
      else
       XFLUX(I,K)=XFLUX(I,K)+XFLUXV/DZ(K)*ART(I)&
                 +DRHOX(I,K)-COR(I)*V(I,K)*DT1(I)*ART(I)
       YFLUX(I,K)=YFLUX(I,K)+YFLUXV/DZ(K)*ART(I)&
                 +DRHOY(I,K)+COR(I)*U(I,K)*DT1(I)*ART(I)
      endif

     END DO
    END IF
   END DO

  if  (SPHERICAL==1)then
  if  (NORTHPOLE==1)then
   CALL ADV_UV_EDGE_XY(XFLUX,YFLUX)
  endif
  endif
   
   DO I=1,N
     IF(ISBCE(I) == 2) THEN
       DO K=1,KBM1
         XFLUX(I,K)=0.0_SP
         YFLUX(I,K)=0.0_SP
       END DO
     END IF
   END DO


   !ADJUST FLUX AT RIVER INFLOWS
   IF(NUMQBC >= 1) THEN
     IF(INFLOW_TYPE == 'node') THEN
       DO II=1,NUMQBC
         J=INODEQ(II)
         I1=NBVE(J,1)
         I2=NBVE(J,NTVE(J))
         DO K=1,KBM1
           VLCTYQ(II)=QDIS(II)/QAREA(II)
           TEMP=0.5_SP*QDIS(II)*VQDIST(II,K)*VLCTYQ(II)
           XFLUX(I1,K)=XFLUX(I1,K)-TEMP/DZ(K)*COS(ANGLEQ(II))
           XFLUX(I2,K)=XFLUX(I2,K)-TEMP/DZ(K)*COS(ANGLEQ(II))
           YFLUX(I1,K)=YFLUX(I1,K)-TEMP/DZ(K)*SIN(ANGLEQ(II))
           YFLUX(I2,K)=YFLUX(I2,K)-TEMP/DZ(K)*SIN(ANGLEQ(II))
         END DO
       END DO
     ELSE IF(INFLOW_TYPE == 'edge') THEN
       DO II=1,NUMQBC
         I1=ICELLQ(II)
         DO K=1,KBM1
           VLCTYQ(II)=QDIS(II)/QAREA(II)
           TEMP=QDIS(II)*VQDIST(II,K)*VLCTYQ(II)
           XFLUX(I1,K)=XFLUX(I1,K)-TEMP/DZ(K)*COS(ANGLEQ(II))
           YFLUX(I1,K)=YFLUX(I1,K)-TEMP/DZ(K)*SIN(ANGLEQ(II))
         END DO
       END DO
     ELSE
       PRINT*,'INFLOW_TYPE NOT CORRECT'
       CALL PSTOP
     END IF
   END IF

   DO I=1,N
    IF(ISWETCT(I)*ISWETC(I) == 1)THEN
!JQI<
 if  (SPHERICAL==1)then
 if  (NORTHPOLE==1)then
   IF(CELL_NORTHAREA(I) == 1)THEN
     DO K=1,KBM1
       U_TMP = -V(I,K)*COS(XC(I)*PI/180.)-U(I,K)*SIN(XC(I)*PI/180.)
       V_TMP = -V(I,K)*SIN(XC(I)*PI/180.)+U(I,K)*COS(XC(I)*PI/180.)
       UF_TMP=U_TMP*DT1(I)/D1(I)-DTI*XFLUX(I,K)/ART(I)/D1(I)
       VF_TMP=V_TMP*DT1(I)/D1(I)-DTI*YFLUX(I,K)/ART(I)/D1(I)
			     
       UF(I,K)  = VF_TMP*COS(XC(I)*PI/180.)-UF_TMP*SIN(XC(I)*PI/180.)
       VF(I,K)  = UF_TMP*COS(XC(I)*PI/180.)+VF_TMP*SIN(XC(I)*PI/180.)
       VF(I,K)  = -VF(I,K)			    
			     
     END DO
   ELSE
    
!JQI>
     DO K=1,KBM1
       UF(I,K)=U(I,K)*DT1(I)/D1(I)-DTI*XFLUX(I,K)/ART(I)/D1(I)
       VF(I,K)=V(I,K)*DT1(I)/D1(I)-DTI*YFLUX(I,K)/ART(I)/D1(I)
     END DO
!JQI<
   END IF
endif
else
 DO K=1,KBM1
       UF(I,K)=U(I,K)*DT1(I)/D1(I)-DTI*XFLUX(I,K)/ART(I)/D1(I)
       VF(I,K)=V(I,K)*DT1(I)/D1(I)-DTI*YFLUX(I,K)/ART(I)/D1(I)
     END DO


endif       
!JQI>
    ELSE
     DO K=1,KBM1
       UF(I,K)=0.0_SP
       VF(I,K)=0.0_SP
     END DO
    END IF
   END DO

   RETURN
   END SUBROUTINE ADV_UV_EDGE_GCY
!==============================================================================!
