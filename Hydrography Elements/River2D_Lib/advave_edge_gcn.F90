!==============================================================================|
!   CALCULATE CONVECTION AND DIFFUSION FLUXES FOR EXTERNAL MODE                !
!==============================================================================|
   SUBROUTINE ADVAVE_EDGE_GCN(XFLUX,YFLUX)
!==============================================================================|

   USE ALL_VARS
   USE MOD_SPHERICAL
   USE MOD_NORTHPOLE
   USE BCS
   USE MOD_OBCS
   USE MOD_WD
   USE MOD_EQUITIDE
   use wave
!QXU{
   USE MOD_BALANCE_2D
!QXU}

   USE MOD_MEANFLOW
   USE MOD_OBCS3

   IMPLICIT NONE
   INTEGER  :: I,J,K,IA,IB,J1,J2,K1,K2,K3,I1,I2,l1,l2,ierr
   REAL(SP) :: DIJ,ELIJ,XIJ,YIJ,UIJ,VIJ,UIJ1,VIJ1,UIJ2,VIJ2
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: FXX,FYY,XADV,YADV,TXXIJ,TYYIJ,TXYIJ,UN_TMP
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP
   REAL(SP) :: XFLUX(0:NT),YFLUX(0:NT)
   REAL(SP) :: FACT,FM1,ISWETTMP, ufact, pre1,pre2,preij

   REAL(SP) :: TPA,TPB

   REAL(DP) :: XTMP,XTMP1

!------------------------------------------------------------------------------!

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

!
!-------------------------INITIALIZE FLUXES------------------------------------!
!
   XFLUX = 0.0_SP
   YFLUX = 0.0_SP
   PSTX  = 0.0_SP
   PSTY  = 0.0_SP
!QXU{
  if  (BALANCE_2D==1)then

   ADFXA =0.0_SP
   ADFYA =0.0_SP
  endif
!QXU}

!
!-------------------------ACCUMULATE FLUX OVER ELEMENT EDGES-------------------!
!
   DO I=1,NE
     IA=IEC(I,1)               !   第i条边的相邻两单元
     IB=IEC(I,2)
     J1=IENODE(I,1)            !   第i条边的相邻两结点 
     J2=IENODE(I,2)
     DIJ=0.5_SP*(D(J1)+D(J2))
     ELIJ=0.5_SP*(EL(J1)+EL(J2))

	 if(icyclone==1)then      
     CALL BRACKET(WND_TM,THOUR,L1,L2,FACT,UFACT,IERR)    ! zhangzhuo
     pre1 = UFACT*pre(J1,L1) + FACT*pre(J1,L2)
     pre2 = ufact*pre(j2,L1) + fact*pre(j2,l2)
	 preij=0.5*(pre1+pre2)
     endif

    if  (EQUI_TIDE==1)then
     ELIJ=ELIJ-0.5_SP*(EL_EQI(J1)+EL_EQI(J2))
     endif
     IF(ISWETCE(IA)*ISWETC(IA) == 1 .OR. ISWETCE(IB)*ISWETC(IB) == 1)THEN
!    FLUX FROM LEFT
     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
         
     COFA1=A1U(IA,1)*UA(IA)+A1U(IA,2)*UA(K1)+A1U(IA,3)*UA(K2)+A1U(IA,4)*UA(K3)
     COFA2=A2U(IA,1)*UA(IA)+A2U(IA,2)*UA(K1)+A2U(IA,3)*UA(K2)+A2U(IA,4)*UA(K3)
     COFA5=A1U(IA,1)*VA(IA)+A1U(IA,2)*VA(K1)+A1U(IA,3)*VA(K2)+A1U(IA,4)*VA(K3)
     COFA6=A2U(IA,1)*VA(IA)+A2U(IA,2)*VA(K1)+A2U(IA,3)*VA(K2)+A2U(IA,4)*VA(K3)
     
    if  (SPHERICAL==1)then
     UIJ1=UA(IA)+COFA1*DLTXNE(I,1)+COFA2*DLTYNE(I,1)
     VIJ1=VA(IA)+COFA5*DLTXNE(I,1)+COFA6*DLTYNE(I,1)
    else
     XIJ=XIJC(I)-XC(IA)         ! xijc 边中点
     YIJ=YIJC(I)-YC(IA)
     UIJ1=UA(IA)+COFA1*XIJ+COFA2*YIJ
     VIJ1=VA(IA)+COFA5*XIJ+COFA6*YIJ
    endif

!    FLUX FROM RIGHT
     K1=NBE(IB,1)
     K2=NBE(IB,2)
     K3=NBE(IB,3)
          
     COFA3=A1U(IB,1)*UA(IB)+A1U(IB,2)*UA(K1)+A1U(IB,3)*UA(K2)+A1U(IB,4)*UA(K3)
     COFA4=A2U(IB,1)*UA(IB)+A2U(IB,2)*UA(K1)+A2U(IB,3)*UA(K2)+A2U(IB,4)*UA(K3)
     COFA7=A1U(IB,1)*VA(IB)+A1U(IB,2)*VA(K1)+A1U(IB,3)*VA(K2)+A1U(IB,4)*VA(K3)
     COFA8=A2U(IB,1)*VA(IB)+A2U(IB,2)*VA(K1)+A2U(IB,3)*VA(K2)+A2U(IB,4)*VA(K3)
     
    if  (SPHERICAL==1) then
     UIJ2=UA(IB)+COFA3*DLTXNE(I,2)+COFA4*DLTYNE(I,2)
     VIJ2=VA(IB)+COFA7*DLTXNE(I,2)+COFA8*DLTYNE(I,2)
    else
     XIJ=XIJC(I)-XC(IB)
     YIJ=YIJC(I)-YC(IB)
     UIJ2=UA(IB)+COFA3*XIJ+COFA4*YIJ
     VIJ2=VA(IB)+COFA7*XIJ+COFA8*YIJ
    endif

!    NORMAL VELOCITY
     UIJ=0.5_SP*(UIJ1+UIJ2)
     VIJ=0.5_SP*(VIJ1+VIJ2)
     UN_TMP=-UIJ*DLTYC(I) + VIJ*DLTXC(I)

!    VISCOSITY COEFFICIENT
     VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
     VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)
!     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2) + FM1)/HPRNU

!    SHEAR STRESSES
     TXXIJ=(COFA1+COFA3)*VISCOF
     TYYIJ=(COFA6+COFA8)*VISCOF
     TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
     FXX=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
     FYY=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))

!    ADD CONVECTIVE AND VISCOUS FLUXES
     XADV=DIJ*UN_TMP*&
          ((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1)*0.5_SP
     YADV=DIJ*UN_TMP* &
          ((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1)*0.5_SP

!    ACCUMULATE FLUX
   if  (MEAN_FLOW/=1)then
     XFLUX(IA)=XFLUX(IA)+(XADV+FXX*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
     YFLUX(IA)=YFLUX(IA)+(YADV+FYY*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
     XFLUX(IB)=XFLUX(IB)-(XADV+FXX*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
     YFLUX(IB)=YFLUX(IB)-(YADV+FYY*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
  else
     XFLUX(IA)=XFLUX(IA)+(XADV+FXX)*(1.0_SP-ISBC(I))*IUCP(IA)
     YFLUX(IA)=YFLUX(IA)+(YADV+FYY)*(1.0_SP-ISBC(I))*IUCP(IA)
     XFLUX(IB)=XFLUX(IB)-(XADV+FXX)*(1.0_SP-ISBC(I))*IUCP(IB)
     YFLUX(IB)=YFLUX(IB)-(YADV+FYY)*(1.0_SP-ISBC(I))*IUCP(IB)
    endif


     END IF

!QXU{
  if  (BALANCE_2D==1)then
     ADFXA(IA)=ADFXA(IA)+FXX
     ADFYA(IA)=ADFYA(IA)+FYY
     ADFXA(IB)=ADFXA(IB)-FXX
     ADFYA(IB)=ADFYA(IB)-FYY
  endif
!QXU}


!    ACCUMULATE BAROTROPIC FLUX
!JQI for spherical coordinator and domain across 360^o latitude         
    if  (SPHERICAL==1)then
     XTMP  = VX(J2)*TPI-VX(J1)*TPI
     XTMP1 = VX(J2)-VX(J1)
     IF(XTMP1 >  180.0_SP)THEN
       XTMP = -360.0_SP*TPI+XTMP
     ELSE IF(XTMP1 < -180.0_SP)THEN
       XTMP =  360.0_SP*TPI+XTMP
     END IF  
                                                                                   
     PSTX(IA)=PSTX(IA)-GRAV*D1(IA)*ELIJ*DLTYC(I)                                   !  正压项  
     PSTY(IA)=PSTY(IA)+GRAV*D1(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))

     PSTX(IB)=PSTX(IB)+GRAV*D1(IB)*ELIJ*DLTYC(I)
     PSTY(IB)=PSTY(IB)-GRAV*D1(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
    else
!JQI>
     PSTX(IA)=PSTX(IA)-GRAV*D1(IA)*ELIJ*DLTYC(I)
     PSTY(IA)=PSTY(IA)+GRAV*D1(IA)*ELIJ*DLTXC(I)
     PSTX(IB)=PSTX(IB)+GRAV*D1(IB)*ELIJ*DLTYC(I)
     PSTY(IB)=PSTY(IB)-GRAV*D1(IB)*ELIJ*DLTXC(I)



!--------------------------------------------------------------------------------------------------------
       if(icyclone ==1) then                       !  增加气压梯度计算 zhangzhuo
    pstx(ia)=pstx(ia)-D1(IA)/1000.*preij*dltyc(i)
    psty(ia)=psty(ia)+D1(IA)/1000.*preij*dltxc(i)
    pstx(ib)=pstx(ib)+D1(ib)/1000.*preij*dltyc(i)
    psty(ib)=psty(ib)-D1(ib)/1000.*preij*dltxc(i)
        endif
		if(iwave==1)then 
     pstx(ia)=pstx(ia)+(-radxxij(i)*dltyc(i)+radxyij(i)*dltxc(i))
     psty(ia)=psty(ia)+(-radxyij(i)*dltyc(i)+radyyij(i)*dltxc(i))
     pstx(ib)=pstx(ib)+(radxxij(i)*dltyc(i)-radxyij(i)*dltxc(i))
     psty(ib)=psty(ib)+(radxyij(i)*dltyc(i)-radyyij(i)*dltxc(i))

       endif
!-----------------------------------------------------------------------------------------------------------



!JQI<
    endif     
!JQI>
   END DO

  if (SPHERICAL==1)then
  if (NORTHPOLE==1)then
   CALL ADVAVE_EDGE_XY(XFLUX,YFLUX)
  endif 
  endif  

   DO I = 1,N
     ISWETTMP = ISWETCE(I)*ISWETC(I)
     XFLUX(I) = XFLUX(I)*ISWETTMP
     YFLUX(I) = YFLUX(I)*ISWETTMP
   END DO


!
!-------------------------SET BOUNDARY VALUES----------------------------------!
!

  if  (MEAN_FLOW/=1)then
      DO I=1,N
        IF(ISBCE(I) == 2) THEN
          XFLUX(I)=(XFLUX(I)+Fluxobn(I)*UA(I))*IUCP(I)
          YFLUX(I)=(YFLUX(I)+Fluxobn(I)*VA(I))*IUCP(I)
        ENDIF
      END DO
  else
   IF (nmfcell > 0) THEN
     DO K=1,nmfcell
       I1=I_MFCELL_N(K)
       XFLUX(I1) = XFLUX(I1) + FLUXOBC2D_X(K)*IUCP(I1)
       YFLUX(I1) = YFLUX(I1) + FLUXOBC2D_Y(K)*IUCP(I1)
     END DO
   END IF
  endif 



!  ADJUST FLUX FOR RIVER INFLOW
   IF(NUMQBC > 0) THEN
     IF(INFLOW_TYPE == 'node')THEN
       DO K=1,NUMQBC
         J=INODEQ(K)
         I1=NBVE(J,1)
         I2=NBVE(J,NTVE(J))
         VLCTYQ(K)=QDIS(K)/QAREA(K)
         XFLUX(I1)=XFLUX(I1)-0.5_SP*QDIS(K)*VLCTYQ(K)*COS(ANGLEQ(K))
         YFLUX(I1)=YFLUX(I1)-0.5_SP*QDIS(K)*VLCTYQ(K)*SIN(ANGLEQ(K))
         XFLUX(I2)=XFLUX(I2)-0.5_SP*QDIS(K)*VLCTYQ(K)*COS(ANGLEQ(K))
         YFLUX(I2)=YFLUX(I2)-0.5_SP*QDIS(K)*VLCTYQ(K)*SIN(ANGLEQ(K))
       END DO
     ELSE IF(INFLOW_TYPE == 'edge') THEN
       DO K=1,NUMQBC
         I1=ICELLQ(K)
         VLCTYQ(K)=QDIS(K)/QAREA(K)
         TEMP=QDIS(K)*VLCTYQ(K)
         XFLUX(I1)=XFLUX(I1)-TEMP*COS(ANGLEQ(K))
         YFLUX(I1)=YFLUX(I1)-TEMP*SIN(ANGLEQ(K))
       END DO
     END IF
   END IF

!  ADJUST FLUX FOR OPEN BOUNDARY MEAN FLOW
  if  (MEAN_FLOW==1)then
   IF(nmfcell > 0) THEN
     DO K=1,nmfcell
       I1=I_MFCELL_N(K)
       VLCTYMF(K)=MFQDIS(K)/MFAREA(K)
       TEMP=MFQDIS(K)*VLCTYMF(K)
       XFLUX(I1)=XFLUX(I1)-TEMP*COS(ANGLEMF(K))
       YFLUX(I1)=YFLUX(I1)-TEMP*SIN(ANGLEMF(K))
     END DO
   END IF
  endif

   RETURN
   END SUBROUTINE ADVAVE_EDGE_GCN
!==============================================================================|
