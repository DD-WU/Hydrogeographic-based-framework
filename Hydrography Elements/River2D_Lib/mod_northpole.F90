MODULE MOD_NORTHPOLE
   USE ALL_VARS 
   USE MOD_SPHERICAL
   USE MOD_EQUITIDE
!  if(ATMO_TIDE==1) USE MOD_ATMOTIDE

   IMPLICIT NONE
   SAVE
   INTEGER :: NODE_NORTHPOLE            !Node index at the north pole point
   INTEGER :: MP,NP,NPE,NPCV
   INTEGER, ALLOCATABLE :: NODE_NORTHAREA(:)
   INTEGER, ALLOCATABLE :: CELL_NORTHAREA(:)   
   INTEGER, ALLOCATABLE :: NPEDGE_LST(:)   
   INTEGER, ALLOCATABLE :: NCEDGE_LST(:)   
   INTEGER, ALLOCATABLE :: MP_LST(:),NP_LST(:)   
   REAL(DP), ALLOCATABLE :: A1U_XY(:,:),A2U_XY(:,:)
   REAL(DP), ALLOCATABLE :: AW0_XY(:,:),AWX_XY(:,:),AWY_XY(:,:)

   CONTAINS
!==============================================================================|
     SUBROUTINE FIND_NORTHPOLE
     IMPLICIT NONE
     INTEGER :: I,ITMP,NODE_NORTHPOLE_GL
     INTEGER,ALLOCATABLE :: TMP(:)
     
     NODE_NORTHPOLE = 0
     NODE_NORTHPOLE_GL = 0
     
     DO I = 1,MGL
       IF(ABS(YG(I)-90.0_SP) < 1.0E-4)THEN
         NODE_NORTHPOLE_GL = I
       END IF
     END DO
     
     IF(NODE_NORTHPOLE_GL == 0)THEN
       PRINT*,"NO NODE POINT ON THE NORTH POLE."
       PRINT*,"PLEASE MOVE ONE NODE TO THE NORTH POLE."
       PRINT*,"STOP RUNNING......"
       CALL PSTOP
     END IF  
     
     DO I = 1,MT
       IF(ABS(VY(I)-90.0_SP) < 1.0E-4)THEN
         NODE_NORTHPOLE = I
       END IF
     END DO
     print*,'NORTH POLE = ',NODE_NORTHPOLE
     
     
     ALLOCATE(NODE_NORTHAREA(0:MT)); NODE_NORTHAREA = 0
     ALLOCATE(CELL_NORTHAREA(0:NT)); CELL_NORTHAREA = 0

     IF(NODE_NORTHPOLE >0) THEN
       ITMP = NODE_NORTHPOLE
       ALLOCATE(TMP(MT)); TMP = 0     
       MP = 0
       DO I=1,NTSN(ITMP)-1
	 MP = MP + 1
	 TMP(MP) = NBSN(ITMP,I)
         NODE_NORTHAREA(NBSN(ITMP,I)) = 1
       END DO  
       MP = MP + 1
       TMP(MP) = ITMP
       NODE_NORTHAREA(ITMP) = 1

       ALLOCATE(MP_LST(MP))
       MP_LST(1:MP) = TMP(1:MP)
       DEALLOCATE(TMP)
       
       ALLOCATE(TMP(NT)); TMP = 0          
       NP = 0 
       DO I=1,NTVE(ITMP)
         NP = NP + 1
	 TMP(NP) = NBVE(ITMP,I)
         CELL_NORTHAREA(NBVE(ITMP,I)) = 1
       END DO
     
       ALLOCATE(NP_LST(NP))
       NP_LST(1:NP) = TMP(1:NP)
       DEALLOCATE(TMP)
       
       print*,"MP,NP=",MP,NP
       print*,MP_LST
       print*,NP_LST
     ENDIF
     
     RETURN
     END SUBROUTINE FIND_NORTHPOLE
!==============================================================================|
     
!==============================================================================|

     SUBROUTINE FIND_CELLSIDE
     
     IMPLICIT NONE
     INTEGER  ::  I,IA,IB
     INTEGER, ALLOCATABLE :: TEMP(:)
     
     ALLOCATE(TEMP(NE));  TEMP = ZERO
     NPE = 0
     
     DO I=1,NE
       IA = IEC(I,1)
       IB = IEC(I,2)
       IF(CELL_NORTHAREA(IA) == 1 .OR. CELL_NORTHAREA(IB) == 1)THEN
         NPE = NPE + 1
	 TEMP(NPE) = I
       END IF
     END DO
     
     ALLOCATE(NPEDGE_LST(NPE))
     NPEDGE_LST(1:NPE) = TEMP(1:NPE)
     DEALLOCATE(TEMP)
     
     ALLOCATE(TEMP(NCV));  TEMP = ZERO
     NPCV = 0
     
     DO I=1,NCV
       IA = NIEC(I,1)
       IB = NIEC(I,2)
       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         NPCV = NPCV + 1
	 TEMP(NPCV) = I
       END IF
     END DO
     
     ALLOCATE(NCEDGE_LST(NPCV))
     NCEDGE_LST(1:NPCV) = TEMP(1:NPCV)
     DEALLOCATE(TEMP)
     
     RETURN
     END SUBROUTINE FIND_CELLSIDE
       	 
!==============================================================================|

!==============================================================================|
   SUBROUTINE ADVAVE_EDGE_XY(XFLUX,YFLUX)

   USE MOD_OBCS
   IMPLICIT NONE
   INTEGER  :: I,J,K,IA,IB,J1,J2,K1,K2,K3,I1,I2,II
   REAL(SP) :: DIJ,ELIJ,XIJ,YIJ,UIJ,VIJ,UIJ1,VIJ1,UIJ2,VIJ2
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: FXX,FYY,XADV,YADV,TXXIJ,TYYIJ,TXYIJ 
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP
   REAL(SP) :: XFLUX(0:NT),YFLUX(0:NT)
   REAL(SP) :: FACT,FM1,ISWETTMP

   REAL(SP) :: TPA,TPB

   REAL(SP) :: UIJ1_TMP,VIJ1_TMP,UIJ2_TMP,VIJ2_TMP,TXXIJ_TMP,TYYIJ_TMP
   REAL(SP) :: XADV_TMP,YADV_TMP,PSTX_TMP,PSTY_TMP
   REAL(SP) :: UIJ_TMP,VIJ_TMP,UN_TMP
   REAL(SP) :: DLTXC_TMP,DLTYC_TMP
   REAL(SP) :: VX1_TMP,VX2_TMP,VY1_TMP,VY2_TMP
   REAL(SP) :: UAIA,VAIA,UAIB,VAIB,UAK1,VAK1,UAK2,VAK2,UAK3,VAK3
   REAL(SP) :: XIJC_TMP,YIJC_TMP,XCIA_TMP,YCIA_TMP,XCIB_TMP,YCIB_TMP
   REAL(SP) :: XIJ_TMP,YIJ_TMP
  
!------------------------------------------------------------------------------!

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

!
!--Initialize Variables--------------------------------------------------------|
!
   DO II=1,NPE
     I=NPEDGE_LST(II)
     IA=IEC(I,1)
     IB=IEC(I,2)
     IF(CELL_NORTHAREA(IA) == 1)THEN
       XFLUX(IA) = 0.0_SP
       YFLUX(IA) = 0.0_SP
       PSTX(IA)  = 0.0_SP
       PSTY(IA)  = 0.0_SP
     END IF  
     IF(CELL_NORTHAREA(IB) == 1)THEN  
       XFLUX(IB) = 0.0_SP
       YFLUX(IB) = 0.0_SP
       PSTX(IB)  = 0.0_SP
       PSTY(IB)  = 0.0_SP
     END IF  
   END DO  
!
!-------------------------ACCUMULATE FLUX OVER ELEMENT EDGES-------------------!
!
   DO II=1,NPE
     I=NPEDGE_LST(II)
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)
     DIJ=0.5_SP*(D(J1)+D(J2))
     ELIJ=0.5_SP*(EL(J1)+EL(J2))
      ELIJ=ELIJ-0.5_SP*(EL_EQI(J1)+EL_EQI(J2))

 !   if (ATMO_TIDE==1)    ELIJ=ELIJ-0.5_SP*(EL_ATMO(J1)+EL_ATMO(J2))
 

!    FLUX FROM LEFT
     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
         
     UAIA = -VA(IA)*COS(XC(IA)*PI/180.)-UA(IA)*SIN(XC(IA)*PI/180.)
     VAIA = -VA(IA)*SIN(XC(IA)*PI/180.)+UA(IA)*COS(XC(IA)*PI/180.)
     UAK1 = -VA(K1)*COS(XC(K1)*PI/180.)-UA(K1)*SIN(XC(K1)*PI/180.)
     VAK1 = -VA(K1)*SIN(XC(K1)*PI/180.)+UA(K1)*COS(XC(K1)*PI/180.)
     UAK2 = -VA(K2)*COS(XC(K2)*PI/180.)-UA(K2)*SIN(XC(K2)*PI/180.)
     VAK2 = -VA(K2)*SIN(XC(K2)*PI/180.)+UA(K2)*COS(XC(K2)*PI/180.)
     UAK3 = -VA(K3)*COS(XC(K3)*PI/180.)-UA(K3)*SIN(XC(K3)*PI/180.)
     VAK3 = -VA(K3)*SIN(XC(K3)*PI/180.)+UA(K3)*COS(XC(K3)*PI/180.)
     
     COFA1=A1U_XY(IA,1)*UAIA+A1U_XY(IA,2)*UAK1   &
          +A1U_XY(IA,3)*UAK2+A1U_XY(IA,4)*UAK3
     COFA2=A2U_XY(IA,1)*UAIA+A2U_XY(IA,2)*UAK1   &
          +A2U_XY(IA,3)*UAK2+A2U_XY(IA,4)*UAK3
     COFA5=A1U_XY(IA,1)*VAIA+A1U_XY(IA,2)*VAK1   &
          +A1U_XY(IA,3)*VAK2+A1U_XY(IA,4)*VAK3
     COFA6=A2U_XY(IA,1)*VAIA+A2U_XY(IA,2)*VAK1   &
          +A2U_XY(IA,3)*VAK2+A2U_XY(IA,4)*VAK3
     
     XIJC_TMP = REARTH * COS(YIJC(I)*PI/180.0_SP) * COS(XIJC(I)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YIJC(I)*PI/180.0_SP))
     YIJC_TMP = REARTH * COS(YIJC(I)*PI/180.0_SP) * SIN(XIJC(I)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YIJC(I)*PI/180.0_SP))
     XCIA_TMP = REARTH * COS(YC(IA)*PI/180.0_SP) * COS(XC(IA)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IA)*PI/180.0_SP))
     YCIA_TMP = REARTH * COS(YC(IA)*PI/180.0_SP) * SIN(XC(IA)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IA)*PI/180.0_SP))
		  
     XIJ_TMP = XIJC_TMP-XCIA_TMP
     YIJ_TMP = YIJC_TMP-YCIA_TMP

     UIJ1=UAIA+COFA1*XIJ_TMP+COFA2*YIJ_TMP
     VIJ1=VAIA+COFA5*XIJ_TMP+COFA6*XIJ_TMP

!    FLUX FROM RIGHT
     K1=NBE(IB,1)
     K2=NBE(IB,2)
     K3=NBE(IB,3)
          
     UAIB = -VA(IB)*COS(XC(IB)*PI/180.)-UA(IB)*SIN(XC(IB)*PI/180.)
     VAIB = -VA(IB)*SIN(XC(IB)*PI/180.)+UA(IB)*COS(XC(IB)*PI/180.)
     UAK1 = -VA(K1)*COS(XC(K1)*PI/180.)-UA(K1)*SIN(XC(K1)*PI/180.)
     VAK1 = -VA(K1)*SIN(XC(K1)*PI/180.)+UA(K1)*COS(XC(K1)*PI/180.)
     UAK2 = -VA(K2)*COS(XC(K2)*PI/180.)-UA(K2)*SIN(XC(K2)*PI/180.)
     VAK2 = -VA(K2)*SIN(XC(K2)*PI/180.)+UA(K2)*COS(XC(K2)*PI/180.)
     UAK3 = -VA(K3)*COS(XC(K3)*PI/180.)-UA(K3)*SIN(XC(K3)*PI/180.)
     VAK3 = -VA(K3)*SIN(XC(K3)*PI/180.)+UA(K3)*COS(XC(K3)*PI/180.)
     
     COFA3=A1U_XY(IB,1)*UAIB+A1U_XY(IB,2)*UAK1   &
          +A1U_XY(IB,3)*UAK2+A1U_XY(IB,4)*UAK3
     COFA4=A2U_XY(IB,1)*UAIB+A2U_XY(IB,2)*UAK1   &
          +A2U_XY(IB,3)*UAK2+A2U_XY(IB,4)*UAK3
     COFA7=A1U_XY(IB,1)*VAIB+A1U_XY(IB,2)*VAK1   &
          +A1U_XY(IB,3)*VAK2+A1U_XY(IB,4)*VAK3
     COFA8=A2U_XY(IB,1)*VAIB+A2U_XY(IB,2)*VAK1   &
          +A2U_XY(IB,3)*VAK2+A2U_XY(IB,4)*VAK3
     
     XCIB_TMP = REARTH * COS(YC(IB)*PI/180.0_SP) * COS(XC(IB)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IB)*PI/180.0_SP))
     YCIB_TMP = REARTH * COS(YC(IB)*PI/180.0_SP) * SIN(XC(IB)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IB)*PI/180.0_SP))
		  
     XIJ_TMP = XIJC_TMP-XCIB_TMP
     YIJ_TMP = YIJC_TMP-YCIB_TMP
     UIJ2=UAIB+COFA3*XIJ_TMP+COFA4*YIJ_TMP
     VIJ2=VAIB+COFA7*XIJ_TMP+COFA8*YIJ_TMP

!    VISCOSITY COEFFICIENT
     VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
     VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)

!     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2) + FM1)/HPRNU


     VX1_TMP = REARTH * COS(VY(IENODE(I,1))*PI/180.0_SP) * COS(VX(IENODE(I,1))*PI/180.0_SP) &
               * 2._SP /(1._SP+sin(VY(IENODE(I,1))*PI/180.0_SP))
     VY1_TMP = REARTH * COS(VY(IENODE(I,1))*PI/180.0_SP) * SIN(VX(IENODE(I,1))*PI/180.0_SP) &
               * 2._SP /(1._SP+sin(VY(IENODE(I,1))*PI/180.0_SP))

     VX2_TMP = REARTH * COS(VY(IENODE(I,2))*PI/180.0_SP) * COS(VX(IENODE(I,2))*PI/180.0_SP) &
               * 2._SP /(1._SP+sin(VY(IENODE(I,2))*PI/180.0_SP))
     VY2_TMP = REARTH * COS(VY(IENODE(I,2))*PI/180.0_SP) * SIN(VX(IENODE(I,2))*PI/180.0_SP) &
               * 2._SP /(1._SP+sin(VY(IENODE(I,2))*PI/180.0_SP))

     DLTXC_TMP = VX2_TMP-VX1_TMP
     DLTYC_TMP = VY2_TMP-VY1_TMP

!    SHEAR STRESSES         
     TXXIJ=(COFA1+COFA3)*VISCOF
     TYYIJ=(COFA6+COFA8)*VISCOF
     TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
     FXX=DIJ*(TXXIJ*DLTYC_TMP-TXYIJ*DLTXC_TMP)
     FYY=DIJ*(TXYIJ*DLTYC_TMP-TYYIJ*DLTXC_TMP)

     UIJ1_TMP = UIJ1
     VIJ1_TMP = VIJ1
     UIJ2_TMP = UIJ2
     VIJ2_TMP = VIJ2

     UIJ_TMP=0.5_SP*(UIJ1+UIJ2)
     VIJ_TMP=0.5_SP*(VIJ1+VIJ2)
     UN_TMP=-UIJ_TMP*DLTYC_TMP + VIJ_TMP*DLTXC_TMP

!    ADD CONVECTIVE AND VISCOUS FLUXES

!    ACCUMULATE FLUX
     XADV_TMP=DIJ*UN_TMP*&
              ((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2_TMP                    &
              +(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1_TMP)*0.5_SP  
     YADV_TMP=DIJ*UN_TMP* &
              ((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2_TMP                    &
	      +(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1_TMP)*0.5_SP  

     IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
       XFLUX(IA)=XFLUX(IA)+(XADV_TMP+FXX*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
       YFLUX(IA)=YFLUX(IA)+(YADV_TMP+FYY*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
       XFLUX(IB)=XFLUX(IB)-(XADV_TMP+FXX*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
       YFLUX(IB)=YFLUX(IB)-(YADV_TMP+FYY*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
     ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
       XFLUX(IA)=XFLUX(IA)+(XADV_TMP+FXX*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
       YFLUX(IA)=YFLUX(IA)+(YADV_TMP+FYY*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
     ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN
       XFLUX(IB)=XFLUX(IB)-(XADV_TMP+FXX*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
       YFLUX(IB)=YFLUX(IB)-(YADV_TMP+FYY*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
     END IF 


!    ACCUMULATE BAROTROPIC FLUX

!ggao-------
!     VX1_TMP = REARTH * COS(VY(IENODE(I,1))*PI/180.0_SP) * COS(VX(IENODE(I,1))*PI/180.0_SP) &
!               * 2._SP /(1._SP+sin(VY(IENODE(I,1))*PI/180.0_SP))
!     VY1_TMP = REARTH * COS(VY(IENODE(I,1))*PI/180.0_SP) * SIN(VX(IENODE(I,1))*PI/180.0_SP) &
!               * 2._SP /(1._SP+sin(VY(IENODE(I,1))*PI/180.0_SP))

!     VX2_TMP = REARTH * COS(VY(IENODE(I,2))*PI/180.0_SP) * COS(VX(IENODE(I,2))*PI/180.0_SP) &
!               * 2._SP /(1._SP+sin(VY(IENODE(I,2))*PI/180.0_SP))
!     VY2_TMP = REARTH * COS(VY(IENODE(I,2))*PI/180.0_SP) * SIN(VX(IENODE(I,2))*PI/180.0_SP) &
!               * 2._SP /(1._SP+sin(VY(IENODE(I,2))*PI/180.0_SP))

!     DLTXC_TMP = VX2_TMP-VX1_TMP
!     DLTYC_TMP = VY2_TMP-VY1_TMP

!   ggao

     IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN

!   H1---->D1   ggao 0309/2006

      PSTX(IA)=PSTX(IA)-GRAV*D1(IA)*ELIJ*DLTYC_TMP/(2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
      PSTY(IA)=PSTY(IA)+GRAV*D1(IA)*ELIJ*DLTXC_TMP/(2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
      PSTX(IB)=PSTX(IB)+GRAV*D1(IB)*ELIJ*DLTYC_TMP/(2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
      PSTY(IB)=PSTY(IB)-GRAV*D1(IB)*ELIJ*DLTXC_TMP/(2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
     ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
      PSTX(IA)=PSTX(IA)-GRAV*D1(IA)*ELIJ*DLTYC_TMP/(2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
      PSTY(IA)=PSTY(IA)+GRAV*D1(IA)*ELIJ*DLTXC_TMP/(2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
     ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
      PSTX(IB)=PSTX(IB)+GRAV*D1(IB)*ELIJ*DLTYC_TMP/(2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
      PSTY(IB)=PSTY(IB)-GRAV*D1(IB)*ELIJ*DLTXC_TMP/(2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
     END IF

   END DO

   RETURN
   END SUBROUTINE ADVAVE_EDGE_XY
!==============================================================================|

              	 

!==============================================================================|
   SUBROUTINE ADVECTION_EDGE_XY(XFLUX,YFLUX)

   IMPLICIT NONE
!   REAL(SP), INTENT(OUT), DIMENSION(0:NT,KB) :: XFLUX,YFLUX
   REAL(SP) :: XFLUX(0:NT,KB),YFLUX(0:NT,KB)
   REAL(SP) :: DIJ
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: FXX,FYY,XADV,YADV,TXXIJ,TYYIJ,TXYIJ,UN
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP,TPA,TPB
   REAL(SP) :: XIJA,YIJA,XIJB,YIJB,UIJ,VIJ,UIJ1,VIJ1,UIJ2,VIJ2
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,IA,IB,J1,J2,K1,K2,K3,K4,K5,K6,K,II,J,I1,I2

   REAL(SP) :: UIJ1_TMP,VIJ1_TMP,UIJ2_TMP,VIJ2_TMP,TXXIJ_TMP,TYYIJ_TMP
   REAL(SP) :: XADV_TMP,YADV_TMP
   REAL(SP) :: UIJ_TMP,VIJ_TMP,UN_TMP
   REAL(SP) :: DLTXC_TMP,DLTYC_TMP
   REAL(SP) :: VX1_TMP,VX2_TMP,VY1_TMP,VY2_TMP
   REAL(SP) :: UIA,VIA,UIB,VIB,UK1,VK1,UK2,VK2,UK3,VK3,UK4,VK4,UK5,VK5,UK6,VK6
   REAL(SP) :: XIJC_TMP,YIJC_TMP,XCIA_TMP,YCIA_TMP,XCIB_TMP,YCIB_TMP

!------------------------------------------------------------------------------|

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

!
!--Initialize Variables--------------------------------------------------------|
!
   DO K = 1,KBM1
     DO II=1,NPE
       I=NPEDGE_LST(II)
       IA=IEC(I,1)
       IB=IEC(I,2)
       IF(CELL_NORTHAREA(IA) == 1)THEN
         XFLUX(IA,K) = 0.0_SP
         YFLUX(IA,K) = 0.0_SP
       END IF
       IF(CELL_NORTHAREA(IB) == 1)THEN	 
         XFLUX(IB,K) = 0.0_SP
         YFLUX(IB,K) = 0.0_SP
       END IF	 
     END DO
   END DO    

!
!--Loop Over Edges and Accumulate Fluxes-For Each Element----------------------|
!

   DO II=1,NPE
     I=NPEDGE_LST(II)
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)
     DIJ= 0.5_SP*(DT(J1)+DT(J2))

     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
     K4=NBE(IB,1)
     K5=NBE(IB,2)
     K6=NBE(IB,3)

     XIJC_TMP = REARTH * COS(YIJC(I)*PI/180.0_SP) * COS(XIJC(I)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YIJC(I)*PI/180.0_SP))
     YIJC_TMP = REARTH * COS(YIJC(I)*PI/180.0_SP) * SIN(XIJC(I)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YIJC(I)*PI/180.0_SP))
     XCIA_TMP = REARTH * COS(YC(IA)*PI/180.0_SP) * COS(XC(IA)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IA)*PI/180.0_SP))
     YCIA_TMP = REARTH * COS(YC(IA)*PI/180.0_SP) * SIN(XC(IA)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IA)*PI/180.0_SP))
     XCIB_TMP = REARTH * COS(YC(IB)*PI/180.0_SP) * COS(XC(IB)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IB)*PI/180.0_SP))
     YCIB_TMP = REARTH * COS(YC(IB)*PI/180.0_SP) * SIN(XC(IB)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IB)*PI/180.0_SP))
		  
     XIJA = XIJC_TMP-XCIA_TMP
     YIJA = YIJC_TMP-YCIA_TMP
     XIJB = XIJC_TMP-XCIB_TMP
     YIJB = YIJC_TMP-YCIB_TMP

     DO K=1,KBM1
       UIA = -V(IA,K)*COS(XC(IA)*PI/180.)-U(IA,K)*SIN(XC(IA)*PI/180.)
       VIA = -V(IA,K)*SIN(XC(IA)*PI/180.)+U(IA,K)*COS(XC(IA)*PI/180.)
       UIB = -V(IB,K)*COS(XC(IB)*PI/180.)-U(IB,K)*SIN(XC(IB)*PI/180.)
       VIB = -V(IB,K)*SIN(XC(IB)*PI/180.)+U(IB,K)*COS(XC(IB)*PI/180.)
       UK1 = -V(K1,K)*COS(XC(K1)*PI/180.)-U(K1,K)*SIN(XC(K1)*PI/180.)
       VK1 = -V(K1,K)*SIN(XC(K1)*PI/180.)+U(K1,K)*COS(XC(K1)*PI/180.)
       UK2 = -V(K2,K)*COS(XC(K2)*PI/180.)-U(K2,K)*SIN(XC(K2)*PI/180.)
       VK2 = -V(K2,K)*SIN(XC(K2)*PI/180.)+U(K2,K)*COS(XC(K2)*PI/180.)
       UK3 = -V(K3,K)*COS(XC(K3)*PI/180.)-U(K3,K)*SIN(XC(K3)*PI/180.)
       VK3 = -V(K3,K)*SIN(XC(K3)*PI/180.)+U(K3,K)*COS(XC(K3)*PI/180.)
       UK4 = -V(K4,K)*COS(XC(K4)*PI/180.)-U(K4,K)*SIN(XC(K4)*PI/180.)
       VK4 = -V(K4,K)*SIN(XC(K4)*PI/180.)+U(K4,K)*COS(XC(K4)*PI/180.)
       UK5 = -V(K5,K)*COS(XC(K5)*PI/180.)-U(K5,K)*SIN(XC(K5)*PI/180.)
       VK5 = -V(K5,K)*SIN(XC(K5)*PI/180.)+U(K5,K)*COS(XC(K5)*PI/180.)
       UK6 = -V(K6,K)*COS(XC(K6)*PI/180.)-U(K6,K)*SIN(XC(K6)*PI/180.)
       VK6 = -V(K6,K)*SIN(XC(K6)*PI/180.)+U(K6,K)*COS(XC(K6)*PI/180.)
       !!FORM THE LEFT FLUX
       COFA1=A1U_XY(IA,1)*UIA+A1U_XY(IA,2)*UK1   &
            +A1U_XY(IA,3)*UK2+A1U_XY(IA,4)*UK3
       COFA2=A2U_XY(IA,1)*UIA+A2U_XY(IA,2)*UK1   &
            +A2U_XY(IA,3)*UK2+A2U_XY(IA,4)*UK3
       COFA5=A1U_XY(IA,1)*VIA+A1U_XY(IA,2)*VK1   &
            +A1U_XY(IA,3)*VK2+A1U_XY(IA,4)*VK3
       COFA6=A2U_XY(IA,1)*VIA+A2U_XY(IA,2)*VK1   &
            +A2U_XY(IA,3)*VK2+A2U_XY(IA,4)*VK3
       UIJ1=UIA+COFA1*XIJA+COFA2*YIJA
       VIJ1=VIA+COFA5*XIJA+COFA6*YIJA

       !!FORM THE RIGHT FLUX
       COFA3=A1U_XY(IB,1)*UIB+A1U_XY(IB,2)*UK4   &
            +A1U_XY(IB,3)*UK5+A1U_XY(IB,4)*UK6
       COFA4=A2U_XY(IB,1)*UIB+A2U_XY(IB,2)*UK4   &
            +A2U_XY(IB,3)*UK5+A2U_XY(IB,4)*UK6
       COFA7=A1U_XY(IB,1)*VIB+A1U_XY(IB,2)*VK4   &
            +A1U_XY(IB,3)*VK5+A1U_XY(IB,4)*VK6
       COFA8=A2U_XY(IB,1)*VIB+A2U_XY(IB,2)*VK4   &
            +A2U_XY(IB,3)*VK5+A2U_XY(IB,4)*VK6
       UIJ2=UIB+COFA3*XIJB+COFA4*YIJB
       VIJ2=VIB+COFA7*XIJB+COFA8*YIJB

       !!    VISCOSITY COEFFICIENT EDGE
       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)
        
!       VISCOF = HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
       VISCOF = HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2) + FM1)/HPRNU

       VX1_TMP = REARTH * COS(VY(IENODE(I,1))*PI/180.0_SP) * COS(VX(IENODE(I,1))*PI/180.0_SP)       &
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,1))*PI/180.0_SP))
       VY1_TMP = REARTH * COS(VY(IENODE(I,1))*PI/180.0_SP) * SIN(VX(IENODE(I,1))*PI/180.0_SP)&
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,1))*PI/180.0_SP))

       VX2_TMP = REARTH * COS(VY(IENODE(I,2))*PI/180.0_SP) * COS(VX(IENODE(I,2))*PI/180.0_SP)&
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,2))*PI/180.0_SP))
       VY2_TMP = REARTH * COS(VY(IENODE(I,2))*PI/180.0_SP) * SIN(VX(IENODE(I,2))*PI/180.0_SP)&
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,2))*PI/180.0_SP))

       DLTXC_TMP = VX2_TMP-VX1_TMP
       DLTYC_TMP = VY2_TMP-VY1_TMP
       
       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX=DIJ*(TXXIJ*DLTYC_TMP-TXYIJ*DLTXC_TMP)
       FYY=DIJ*(TXYIJ*DLTYC_TMP-TYYIJ*DLTXC_TMP)

       UIJ1_TMP = UIJ1
       VIJ1_TMP = VIJ1
       UIJ2_TMP = UIJ2
       VIJ2_TMP = VIJ2
     
       UIJ_TMP=0.5_SP*(UIJ1+UIJ2)
       VIJ_TMP=0.5_SP*(VIJ1+VIJ2)
       UN_TMP=-UIJ_TMP*DLTYC_TMP + VIJ_TMP*DLTXC_TMP

       !!COMPUTE BOUNDARY FLUX AUGMENTERS   
       TPA = FLOAT(1-ISBC(I))*EPOR(IA)
       TPB = FLOAT(1-ISBC(I))*EPOR(IB)


       !!ACCUMULATE THE FLUX
       XADV_TMP=DIJ*UN_TMP*&
                ((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2_TMP                 &
	        +(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1_TMP)*0.5_SP
       YADV_TMP=DIJ*UN_TMP* &
                ((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2_TMP                 &
	          +(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1_TMP)*0.5_SP 
       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         XFLUX(IA,K)=XFLUX(IA,K)+XADV_TMP*TPA+(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IA)
         YFLUX(IA,K)=YFLUX(IA,K)+YADV_TMP*TPA+(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IA)
         XFLUX(IB,K)=XFLUX(IB,K)-XADV_TMP*TPB+(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IB)
         YFLUX(IB,K)=YFLUX(IB,K)-YADV_TMP*TPB+(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IB)
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         XFLUX(IA,K)=XFLUX(IA,K)+XADV_TMP*TPA+(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IA)
         YFLUX(IA,K)=YFLUX(IA,K)+YADV_TMP*TPA+(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IA)
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN 
         XFLUX(IB,K)=XFLUX(IB,K)-XADV_TMP*TPB-(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IB)
         YFLUX(IB,K)=YFLUX(IB,K)-YADV_TMP*TPB-(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IB)
       END IF
     END DO
   END DO

   RETURN
   END SUBROUTINE ADVECTION_EDGE_XY
!==============================================================================|

!==============================================================================!

   SUBROUTINE ADV_UV_EDGE_XY(XFLUX,YFLUX)

   IMPLICIT NONE
   REAL(SP) :: XFLUX(0:NT,KB),YFLUX(0:NT,KB)
   REAL(SP) :: PSTX_TM(0:NT,KB),PSTY_TM(0:NT,KB)
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: FXX,FYY,XADV,YADV,TXXIJ,TYYIJ,TXYIJ,UN
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP,TPA,TPB
   REAL(SP) :: XIJA,YIJA,XIJB,YIJB,UIJ,VIJ,UIJ1,VIJ1,UIJ2,VIJ2
   REAL(SP) :: SITA,DIJ,ELIJ,TMPA,TMPB,TMP,XFLUXV,YFLUXV
   REAL(SP) :: FACT,FM1,EXFLUX,ISWETTMP
   INTEGER  :: I,IA,IB,J1,J2,K1,K2,K3,K4,K5,K6,K,II,J,I1,I2

!   REAL(SP) :: TY
   REAL(SP) :: UIJ1_TMP,VIJ1_TMP,UIJ2_TMP,VIJ2_TMP,UIJ3_TMP,VIJ3_TMP
   REAL(SP) :: U_TMP,V_TMP,UF_TMP,VF_TMP
   REAL(SP) :: TXXIJ_TMP,TYYIJ_TMP
   REAL(SP) :: XADV_TMP,YADV_TMP,PSTX_TMP,PSTY_TMP
   REAL(SP) :: UIJ_TMP,VIJ_TMP,EXFLUX_TMP
   REAL(SP) :: DLTXC_TMP,DLTYC_TMP
   REAL(SP) :: VX1_TMP,VX2_TMP,VY1_TMP,VY2_TMP
   REAL(SP) :: UIA,VIA,UIB,VIB,UK1,VK1,UK2,VK2,UK3,VK3,UK4,VK4,UK5,VK5,UK6,VK6 
   REAL(SP) :: XIJC_TMP,YIJC_TMP,XCIA_TMP,YCIA_TMP,XCIB_TMP,YCIB_TMP

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
!  ggao

    IF(NPE== 0)  RETURN
!   ggao

   DO K = 1,KBM1
     DO II=1,NPE
       I=NPEDGE_LST(II)
       IA=IEC(I,1)
       IB=IEC(I,2)
       IF(CELL_NORTHAREA(IA) == 1)THEN
         XFLUX(IA,K) = 0.0_SP
         YFLUX(IA,K) = 0.0_SP
         PSTX_TM(IA,K) = 0.0_SP
         PSTY_TM(IA,K) = 0.0_SP
       END IF
       IF(CELL_NORTHAREA(IB) == 1)THEN
         XFLUX(IB,K) = 0.0_SP
         YFLUX(IB,K) = 0.0_SP
         PSTX_TM(IB,K) = 0.0_SP
         PSTY_TM(IB,K) = 0.0_SP
       END IF	 
     END DO
   END DO         

!
!-----Loop Over Edges and Accumulate Flux--------------------------------------!
!
   DO II=1,NPE
     I=NPEDGE_LST(II)
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)
     DIJ=0.5_SP*(DT(J1)+DT(J2))
     ELIJ=0.5_SP*(EGF(J1)+EGF(J2))
!JQI<	 
     ELIJ=ELIJ-0.5_SP*(EGF_EQI(J1)+EGF_EQI(J2))

!  if (ATMO_TIDE==1)   ELIJ=ELIJ-0.5_SP*(EGF_ATMO(J1)+EGF_ATMO(J2))
       
!JQI>                     


     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
     K4=NBE(IB,1)
     K5=NBE(IB,2)
     K6=NBE(IB,3)

     XIJC_TMP = REARTH * COS(YIJC(I)*PI/180.0_SP) * COS(XIJC(I)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YIJC(I)*PI/180.0_SP))
     YIJC_TMP = REARTH * COS(YIJC(I)*PI/180.0_SP) * SIN(XIJC(I)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YIJC(I)*PI/180.0_SP))
		  
     XCIA_TMP = REARTH * COS(YC(IA)*PI/180.0_SP) * COS(XC(IA)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IA)*PI/180.0_SP))
     YCIA_TMP = REARTH * COS(YC(IA)*PI/180.0_SP) * SIN(XC(IA)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IA)*PI/180.0_SP))
     XCIB_TMP = REARTH * COS(YC(IB)*PI/180.0_SP) * COS(XC(IB)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IB)*PI/180.0_SP))
     YCIB_TMP = REARTH * COS(YC(IB)*PI/180.0_SP) * SIN(XC(IB)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IB)*PI/180.0_SP))
		  
     XIJA = XIJC_TMP-XCIA_TMP
     YIJA = YIJC_TMP-YCIA_TMP
     XIJB = XIJC_TMP-XCIB_TMP
     YIJB = YIJC_TMP-YCIB_TMP

     DO K=1,KBM1
       UIA = -V(IA,K)*COS(XC(IA)*PI/180.)-U(IA,K)*SIN(XC(IA)*PI/180.)
       VIA = -V(IA,K)*SIN(XC(IA)*PI/180.)+U(IA,K)*COS(XC(IA)*PI/180.)
       UIB = -V(IB,K)*COS(XC(IB)*PI/180.)-U(IB,K)*SIN(XC(IB)*PI/180.)
       VIB = -V(IB,K)*SIN(XC(IB)*PI/180.)+U(IB,K)*COS(XC(IB)*PI/180.)
       UK1 = -V(K1,K)*COS(XC(K1)*PI/180.)-U(K1,K)*SIN(XC(K1)*PI/180.)
       VK1 = -V(K1,K)*SIN(XC(K1)*PI/180.)+U(K1,K)*COS(XC(K1)*PI/180.)
       UK2 = -V(K2,K)*COS(XC(K2)*PI/180.)-U(K2,K)*SIN(XC(K2)*PI/180.)
       VK2 = -V(K2,K)*SIN(XC(K2)*PI/180.)+U(K2,K)*COS(XC(K2)*PI/180.)
       UK3 = -V(K3,K)*COS(XC(K3)*PI/180.)-U(K3,K)*SIN(XC(K3)*PI/180.)
       VK3 = -V(K3,K)*SIN(XC(K3)*PI/180.)+U(K3,K)*COS(XC(K3)*PI/180.)
       UK4 = -V(K4,K)*COS(XC(K4)*PI/180.)-U(K4,K)*SIN(XC(K4)*PI/180.)
       VK4 = -V(K4,K)*SIN(XC(K4)*PI/180.)+U(K4,K)*COS(XC(K4)*PI/180.)
       UK5 = -V(K5,K)*COS(XC(K5)*PI/180.)-U(K5,K)*SIN(XC(K5)*PI/180.)
       VK5 = -V(K5,K)*SIN(XC(K5)*PI/180.)+U(K5,K)*COS(XC(K5)*PI/180.)
       UK6 = -V(K6,K)*COS(XC(K6)*PI/180.)-U(K6,K)*SIN(XC(K6)*PI/180.)
       VK6 = -V(K6,K)*SIN(XC(K6)*PI/180.)+U(K6,K)*COS(XC(K6)*PI/180.)

       COFA1=A1U_XY(IA,1)*UIA+A1U_XY(IA,2)*UK1   &
            +A1U_XY(IA,3)*UK2+A1U_XY(IA,4)*UK3
       COFA2=A2U_XY(IA,1)*UIA+A2U_XY(IA,2)*UK1   &
            +A2U_XY(IA,3)*UK2+A2U_XY(IA,4)*UK3
       COFA5=A1U_XY(IA,1)*VIA+A1U_XY(IA,2)*VK1   &
            +A1U_XY(IA,3)*VK2+A1U_XY(IA,4)*VK3
       COFA6=A2U_XY(IA,1)*VIA+A2U_XY(IA,2)*VK1   &
            +A2U_XY(IA,3)*VK2+A2U_XY(IA,4)*VK3

       UIJ1=UIA+COFA1*XIJA+COFA2*YIJA
       VIJ1=VIA+COFA5*XIJA+COFA6*YIJA

       COFA3=A1U_XY(IB,1)*UIB+A1U_XY(IB,2)*UK4   &
            +A1U_XY(IB,3)*UK5+A1U_XY(IB,4)*UK6
       COFA4=A2U_XY(IB,1)*UIB+A2U_XY(IB,2)*UK4   &
            +A2U_XY(IB,3)*UK5+A2U_XY(IB,4)*UK6
       COFA7=A1U_XY(IB,1)*VIB+A1U_XY(IB,2)*VK4   &
            +A1U_XY(IB,3)*VK5+A1U_XY(IB,4)*VK6
       COFA8=A2U_XY(IB,1)*VIB+A2U_XY(IB,2)*VK4   &
            +A2U_XY(IB,3)*VK5+A2U_XY(IB,4)*VK6

       UIJ2=UIB+COFA3*XIJB+COFA4*YIJB
       VIJ2=VIB+COFA7*XIJB+COFA8*YIJB

!-------ADD THE VISCOUS TERM & ADVECTION TERM---------------------------------!
!
       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)

!       VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON
       VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON/HPRNU

       VX1_TMP = REARTH * COS(VY(IENODE(I,1))*PI/180.0_SP) * COS(VX(IENODE(I,1))*PI/180.0_SP)       &
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,1))*PI/180.0_SP))
       VY1_TMP = REARTH * COS(VY(IENODE(I,1))*PI/180.0_SP) * SIN(VX(IENODE(I,1))*PI/180.0_SP)&
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,1))*PI/180.0_SP))

       VX2_TMP = REARTH * COS(VY(IENODE(I,2))*PI/180.0_SP) * COS(VX(IENODE(I,2))*PI/180.0_SP)&
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,2))*PI/180.0_SP))
       VY2_TMP = REARTH * COS(VY(IENODE(I,2))*PI/180.0_SP) * SIN(VX(IENODE(I,2))*PI/180.0_SP)&
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,2))*PI/180.0_SP))

       DLTXC_TMP = VX2_TMP-VX1_TMP
       DLTYC_TMP = VY2_TMP-VY1_TMP
       
       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX=DIJ*(TXXIJ*DLTYC_TMP-TXYIJ*DLTXC_TMP)
       FYY=DIJ*(TXYIJ*DLTYC_TMP-TYYIJ*DLTXC_TMP)


       UIJ1_TMP = UIJ1
       VIJ1_TMP = VIJ1
       UIJ2_TMP = UIJ2
       VIJ2_TMP = VIJ2

       UIJ_TMP=0.5_SP*(UIJ1+UIJ2)
       VIJ_TMP=0.5_SP*(VIJ1+VIJ2)
       EXFLUX_TMP = DIJ*(-UIJ_TMP*DLTYC_TMP + VIJ_TMP*DLTXC_TMP)

       !!CALCULATE BOUNDARY FLUX AUGMENTERS
       TPA = FLOAT(1-ISBC(I))*EPOR(IA)
       TPB = FLOAT(1-ISBC(I))*EPOR(IB)

       !!ACCUMULATE ADVECTIVE + DIFFUSIVE + BAROTROPIC PRESSURE GRADIENT TERMS
       XADV_TMP=EXFLUX_TMP*&
                ((1.0_SP-SIGN(1.0_SP,EXFLUX_TMP))*UIJ2_TMP                 &
                +(1.0_SP+SIGN(1.0_SP,EXFLUX_TMP))*UIJ1_TMP)*0.5_SP
       YADV_TMP=EXFLUX_TMP* &
                ((1.0_SP-SIGN(1.0_SP,EXFLUX_TMP))*VIJ2_TMP                 &
 	        +(1.0_SP+SIGN(1.0_SP,EXFLUX_TMP))*VIJ1_TMP)*0.5_SP
       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         XFLUX(IA,K)=XFLUX(IA,K)+XADV_TMP*TPA+(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IA)
         YFLUX(IA,K)=YFLUX(IA,K)+YADV_TMP*TPA+(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IA)
         XFLUX(IB,K)=XFLUX(IB,K)-XADV_TMP*TPB+(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IB)
         YFLUX(IB,K)=YFLUX(IB,K)-YADV_TMP*TPB+(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IB)
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         XFLUX(IA,K)=XFLUX(IA,K)+XADV_TMP*TPA+(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IA)
         YFLUX(IA,K)=YFLUX(IA,K)+YADV_TMP*TPA+(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IA)
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN 
         XFLUX(IB,K)=XFLUX(IB,K)-XADV_TMP*TPB-(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IB)
         YFLUX(IB,K)=YFLUX(IB,K)-YADV_TMP*TPB-(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IB)
       END IF

       VX1_TMP = REARTH * COS(VY(IENODE(I,1))*PI/180.0_SP) * COS(VX(IENODE(I,1))*PI/180.0_SP) &
                 * 2._SP /(1._SP+sin(VY(IENODE(I,1))*PI/180.0_SP))
       VY1_TMP = REARTH * COS(VY(IENODE(I,1))*PI/180.0_SP) * SIN(VX(IENODE(I,1))*PI/180.0_SP) &
                 * 2._SP /(1._SP+sin(VY(IENODE(I,1))*PI/180.0_SP))

       VX2_TMP = REARTH * COS(VY(IENODE(I,2))*PI/180.0_SP) * COS(VX(IENODE(I,2))*PI/180.0_SP) &
                 * 2._SP /(1._SP+sin(VY(IENODE(I,2))*PI/180.0_SP))
       VY2_TMP = REARTH * COS(VY(IENODE(I,2))*PI/180.0_SP) * SIN(VX(IENODE(I,2))*PI/180.0_SP) &
                 * 2._SP /(1._SP+sin(VY(IENODE(I,2))*PI/180.0_SP))

       DLTXC_TMP = VX2_TMP-VX1_TMP
       DLTYC_TMP = VY2_TMP-VY1_TMP

       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV*DT1(IA)*ELIJ*DLTYC_TMP/(2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
         PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV*DT1(IA)*ELIJ*DLTXC_TMP/(2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
         PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV*DT1(IB)*ELIJ*DLTYC_TMP/(2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
         PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV*DT1(IB)*ELIJ*DLTXC_TMP/(2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV*D1(IA)*ELIJ*DLTYC_TMP/(2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
         PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV*D1(IA)*ELIJ*DLTXC_TMP/(2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN       
         PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV*DT1(IB)*ELIJ*DLTYC_TMP/(2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
         PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV*DT1(IB)*ELIJ*DLTXC_TMP/(2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
       END IF
     END DO
   END DO

   DO K = 1,KBM1
     DO II=1,NP
       I=NP_LST(II)
       XFLUX(I,K)=XFLUX(I,K)+PSTX_TM(I,K)
       YFLUX(I,K)=YFLUX(I,K)+PSTY_TM(I,K)
     END DO
   END DO

!
!-------ADD VERTICAL CONVECTIVE FLUX, CORIOLIS TERM AND BAROCLINIC PG TERM----!
!
   DO II=1,NP
     I=NP_LST(II)
     IF(CELL_NORTHAREA(I) == 1)THEN
       DO K=1,KBM1
         IF(K == 1) THEN
           UIJ1_TMP = -V(I,K)*COS(XC(I)*PI/180.)-U(I,K)*SIN(XC(I)*PI/180.)
           VIJ1_TMP = -V(I,K)*SIN(XC(I)*PI/180.)+U(I,K)*COS(XC(I)*PI/180.)
           UIJ2_TMP = -V(I,K+1)*COS(XC(I)*PI/180.)-U(I,K+1)*SIN(XC(I)*PI/180.)
           VIJ2_TMP = -V(I,K+1)*SIN(XC(I)*PI/180.)+U(I,K+1)*COS(XC(I)*PI/180.)
           XFLUXV=-W(I,K+1)*(UIJ1_TMP*DZ(K+1)+UIJ2_TMP*DZ(K))/(DZ(K)+DZ(K+1))
           YFLUXV=-W(I,K+1)*(VIJ1_TMP*DZ(K+1)+VIJ2_TMP*DZ(K))/(DZ(K)+DZ(K+1))
         ELSE IF(K == KBM1) THEN
           UIJ1_TMP = -V(I,K)*COS(XC(I)*PI/180.)-U(I,K)*SIN(XC(I)*PI/180.)
           VIJ1_TMP = -V(I,K)*SIN(XC(I)*PI/180.)+U(I,K)*COS(XC(I)*PI/180.)
           UIJ2_TMP = -V(I,K-1)*COS(XC(I)*PI/180.)-U(I,K-1)*SIN(XC(I)*PI/180.)
           VIJ2_TMP = -V(I,K-1)*SIN(XC(I)*PI/180.)+U(I,K-1)*COS(XC(I)*PI/180.)
           XFLUXV= W(I,K)*(UIJ1_TMP*DZ(K-1)+UIJ2_TMP*DZ(K))/(DZ(K)+DZ(K-1))
           YFLUXV= W(I,K)*(VIJ1_TMP*DZ(K-1)+VIJ2_TMP*DZ(K))/(DZ(K)+DZ(K-1))
         ELSE
           UIJ1_TMP = -V(I,K)*COS(XC(I)*PI/180.)-U(I,K)*SIN(XC(I)*PI/180.)
           VIJ1_TMP = -V(I,K)*SIN(XC(I)*PI/180.)+U(I,K)*COS(XC(I)*PI/180.)
           UIJ2_TMP = -V(I,K-1)*COS(XC(I)*PI/180.)-U(I,K-1)*SIN(XC(I)*PI/180.)
           VIJ2_TMP = -V(I,K-1)*SIN(XC(I)*PI/180.)+U(I,K-1)*COS(XC(I)*PI/180.)
           UIJ3_TMP = -V(I,K+1)*COS(XC(I)*PI/180.)-U(I,K+1)*SIN(XC(I)*PI/180.)
           VIJ3_TMP = -V(I,K+1)*SIN(XC(I)*PI/180.)+U(I,K+1)*COS(XC(I)*PI/180.)
           XFLUXV= W(I,K)*(UIJ1_TMP*DZ(K-1)+UIJ2_TMP*DZ(K))/(DZ(K)+DZ(K-1))-  &
                   W(I,K+1)*(UIJ1_TMP*DZ(K+1)+UIJ3_TMP*DZ(K))/(DZ(K)+DZ(K+1))
           YFLUXV= W(I,K)*(VIJ1_TMP*DZ(K-1)+VIJ2_TMP*DZ(K))/(DZ(K)+DZ(K-1))-  &
                   W(I,K+1)*(VIJ1_TMP*DZ(K+1)+VIJ3_TMP*DZ(K))/(DZ(K)+DZ(K+1))
         END IF
         U_TMP = -V(I,K)*COS(XC(I)*PI/180.)-U(I,K)*SIN(XC(I)*PI/180.)
         V_TMP = -V(I,K)*SIN(XC(I)*PI/180.)+U(I,K)*COS(XC(I)*PI/180.)
         XFLUX(I,K)=XFLUX(I,K)+XFLUXV/DZ(K)*ART(I)&
                    +DRHOX(I,K)-COR(I)*V_TMP*DT1(I)*ART(I)
         YFLUX(I,K)=YFLUX(I,K)+YFLUXV/DZ(K)*ART(I)&
                    +DRHOY(I,K)+COR(I)*U_TMP*DT1(I)*ART(I)

       END DO
     END IF
   END DO

   RETURN
   END SUBROUTINE ADV_UV_EDGE_XY
!==============================================================================!

!==============================================================================|
   SUBROUTINE EXTEL_EDGE_XY(K,XFLUX)       

   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP) :: XFLUX(0:MT)
   REAL(SP) :: DIJ,UIJ,VIJ
   INTEGER  :: I,J,K,I1,IA,IB,JJ,J1,J2,II

   REAL(SP) :: UIJ_TMP,VIJ_TMP,EXFLUX_TMP
   REAL(SP) :: DLTXE_TMP,DLTYE_TMP
   REAL(SP) :: VX1_TMP,VX2_TMP,VY1_TMP,VY2_TMP

!----------INITIALIZE FLUX ARRAY ----------------------------------------------!

   DO II=1,NPCV
     I = NCEDGE_LST(II)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)
     IF(IA == NODE_NORTHPOLE)THEN
       XFLUX(IA) = 0.0_SP
     END IF
     IF(IB == NODE_NORTHPOLE)THEN  
       XFLUX(IB) = 0.0_SP
     END IF  
   END DO  

!---------ACCUMULATE FLUX BY LOOPING OVER CONTROL VOLUME HALF EDGES------------!

   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1  = NTRG(I)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)
     DIJ = D1(I1)

     UIJ = UA(I1)
     VIJ = VA(I1)

     IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
       UIJ_TMP = -VIJ*COS(XC(I1)*PI/180.)-UIJ*SIN(XC(I1)*PI/180.)
       VIJ_TMP = -VIJ*SIN(XC(I1)*PI/180.)+UIJ*COS(XC(I1)*PI/180.)
       
       VX1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * COS(XIJE(I,1)*PI/180.0_SP)&
                 * 2._SP /(1._SP+sin(YIJE(I,1)*PI/180.0_SP))
       VY1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * SIN(XIJE(I,1)*PI/180.0_SP)&
                 * 2._SP /(1._SP+sin(YIJE(I,1)*PI/180.0_SP))

       VX2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * COS(XIJE(I,2)*PI/180.0_SP)&
                 * 2._SP /(1._SP+sin(YIJE(I,2)*PI/180.0_SP))
       VY2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * SIN(XIJE(I,2)*PI/180.0_SP)&
                 * 2._SP /(1._SP+sin(YIJE(I,2)*PI/180.0_SP))

       DLTXE_TMP = VX2_TMP-VX1_TMP
       DLTYE_TMP = VY2_TMP-VY1_TMP
       
       EXFLUX_TMP = DIJ*(-UIJ_TMP*DLTYE_TMP+VIJ_TMP*DLTXE_TMP)
     END IF  
     
     IF(IA == NODE_NORTHPOLE) THEN    
       XFLUX(IA) = XFLUX(IA)-EXFLUX_TMP
     ELSE IF(IB == NODE_NORTHPOLE)THEN
       XFLUX(IB) = XFLUX(IB)+EXFLUX_TMP
     END IF    

   END DO

   RETURN
   END SUBROUTINE EXTEL_EDGE_XY
!==============================================================================|



!==============================================================================|

   SUBROUTINE EXTUV_EDGE_XY(K)       

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: K
   REAL(SP), DIMENSION(0:NT) :: RESX,RESY
   INTEGER  :: I,II

   REAL(SP) :: UARK_TMP,VARK_TMP,UAF_TMP,VAF_TMP,UA_TMP,VA_TMP

!ggao 0310/2006
   REAL(SP) :: WUSURF2_TMP,WVSURF2_TMP,WUBOT_TMP,WVBOT_TMP

!
!--ACCUMULATE RESIDUALS FOR EXTERNAL MODE EQUATIONS----------------------------|
!
   DO II=1,NP
     I=NP_LST(II)
     UA_TMP = -VA(I)*COS(XC(I)*PI/180.)-UA(I)*SIN(XC(I)*PI/180.)
     VA_TMP = -VA(I)*SIN(XC(I)*PI/180.)+UA(I)*COS(XC(I)*PI/180.)

!ggao
     WUSURF2_TMP = -WVSURF2(I)*COS(XC(I)*PI/180.)-WUSURF2(I)*SIN(XC(I)*PI/180.)
     WVSURF2_TMP = -WVSURF2(I)*SIN(XC(I)*PI/180.)+WUSURF2(I)*COS(XC(I)*PI/180.)

     WUBOT_TMP = -WVBOT(I)*COS(XC(I)*PI/180.)-WUBOT(I)*SIN(XC(I)*PI/180.)
     WVBOT_TMP = -WVBOT(I)*SIN(XC(I)*PI/180.)+WUBOT(I)*COS(XC(I)*PI/180.)
!ggao  0310/2006
     RESX(I) = ADX2D(I)+ADVUA(I)+DRX2D(I)+PSTX(I)                &
               -COR(I)*VA_TMP*D1(I)*ART(I)                       &
!               -(WUSURF2(I)+WUBOT(I))*ART(I)
               -(WUSURF2_TMP+WUBOT_TMP)*ART(I)
     RESY(I) = ADY2D(I)+ADVVA(I)+DRY2D(I)+PSTY(I)                &
               +COR(I)*UA_TMP*D1(I)*ART(I)                       &
!               -(WVSURF2(I)+WVBOT(I))*ART(I)
               -(WVSURF2_TMP+WVBOT_TMP)*ART(I)

!
!--UPDATE----------------------------------------------------------------------|
!
     UARK_TMP = -VARK(I)*COS(XC(I)*PI/180.)-UARK(I)*SIN(XC(I)*PI/180.)
     VARK_TMP = -VARK(I)*SIN(XC(I)*PI/180.)+UARK(I)*COS(XC(I)*PI/180.)

     UAF_TMP = (UARK_TMP*(H1(I)+ELRK1(I))              &
               -ALPHA_RK(K)*DTE*RESX(I)/ART(I))/(H1(I)+ELF1(I))
     VAF_TMP = (VARK_TMP*(H1(I)+ELRK1(I))              &
               -ALPHA_RK(K)*DTE*RESY(I)/ART(I))/(H1(I)+ELF1(I))
		       
     UAF(I)  = VAF_TMP*COS(XC(I)*PI/180.)-UAF_TMP*SIN(XC(I)*PI/180.)
     VAF(I)  = UAF_TMP*COS(XC(I)*PI/180.)+VAF_TMP*SIN(XC(I)*PI/180.)
     VAF(I)  = -VAF(I)			    

   END DO

   RETURN
   END SUBROUTINE EXTUV_EDGE_XY
!==============================================================================|

!==============================================================================|

   SUBROUTINE VERTVL_EDGE_XY(XFLUX)         

!------------------------------------------------------------------------------|
   IMPLICIT NONE 
   REAL(SP) :: XFLUX(MT,KBM1)
   REAL(SP) :: DIJ,UIJ,VIJ,UN,EXFLUX,TMP1
   INTEGER  :: I,K,IA,IB,I1 ,J,JJ,J1,J2,II

   REAL(SP) :: UIJ_TMP,VIJ_TMP,VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP
   REAL(SP) :: DLTXE_TMP,DLTYE_TMP,EXFLUX_TMP
!------------------------------------------------------------------------------|

!----------------------INITIALIZE FLUX-----------------------------------------!

   DO K=1,KBM1
     DO II=1,NPCV
       I = NCEDGE_LST(II)
       IA  = NIEC(I,1)
       IB  = NIEC(I,2)
       IF(IA == NODE_NORTHPOLE)THEN
         XFLUX(IA,K) = 0.0_SP
       END IF
       IF(IB == NODE_NORTHPOLE)THEN  
         XFLUX(IB,K) = 0.0_SP
       END IF  
     END DO  
   END DO
!----------------------ACCUMULATE FLUX-----------------------------------------!

   DO II=1,NPCV
     I=NCEDGE_LST(II)
     I1=NTRG(I)
     IA=NIEC(I,1)
     IB=NIEC(I,2)
     DIJ=DT1(I1)
     DO K=1,KBM1
       UIJ=U(I1,K)
       VIJ=V(I1,K)

       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         UIJ_TMP = -VIJ*COS(XC(I1)*PI/180.)-UIJ*SIN(XC(I1)*PI/180.)
         VIJ_TMP = -VIJ*SIN(XC(I1)*PI/180.)+UIJ*COS(XC(I1)*PI/180.)



       VX1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * COS(XIJE(I,1)*PI/180.0_SP)&
                 * 2._SP /(1._SP+sin(YIJE(I,1)*PI/180.0_SP))
       VY1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * SIN(XIJE(I,1)*PI/180.0_SP)&
                 * 2._SP /(1._SP+sin(YIJE(I,1)*PI/180.0_SP))

       VX2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * COS(XIJE(I,2)*PI/180.0_SP)&
                 * 2._SP /(1._SP+sin(YIJE(I,2)*PI/180.0_SP))
       VY2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * SIN(XIJE(I,2)*PI/180.0_SP)&
                 * 2._SP /(1._SP+sin(YIJE(I,2)*PI/180.0_SP))

         DLTXE_TMP = VX2_TMP-VX1_TMP
         DLTYE_TMP = VY2_TMP-VY1_TMP
       
         EXFLUX_TMP = DIJ*(-UIJ_TMP*DLTYE_TMP+VIJ_TMP*DLTXE_TMP)
       END IF  
     
       IF(IA == NODE_NORTHPOLE)THEN
         XFLUX(IA,K) = XFLUX(IA,K)-EXFLUX_TMP
       ELSE IF(IB == NODE_NORTHPOLE)THEN
         XFLUX(IB,K) = XFLUX(IB,K)+EXFLUX_TMP
       END IF
     END DO
   END DO

   RETURN
   END SUBROUTINE VERTVL_EDGE_XY
!==============================================================================|

!==============================================================================|
   SUBROUTINE ADV_S_XY(XFLUX,XFLUX_ADV,PSPX,PSPY,PSPXD,PSPYD,VISCOFF,K)               

!------------------------------------------------------------------------------|

   IMPLICIT NONE
   REAL(SP), DIMENSION(0:MT,KB)     :: XFLUX,XFLUX_ADV
   REAL(SP), DIMENSION(M)           :: PSPX,PSPY,PSPXD,PSPYD,VISCOFF
   REAL(SP), DIMENSION(3*(NT))      :: DTIJ 
   REAL(SP) :: XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2 
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF   
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
   REAL(SP) :: TXPI,TYPI

   REAL(SP) :: VX_TMP,VY_TMP,VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP,VX3_TMP,VY3_TMP
   REAL(SP) :: XI_TMP,YI_TMP,VXA_TMP,VYA_TMP,VXB_TMP,VYB_TMP
   REAL(SP) :: UIJ_TMP,VIJ_TMP,DLTXE_TMP,DLTYE_TMP,UVN_TMP,EXFLUX_TMP
   REAL(SP) :: PUPX_TMP,PUPY_TMP,PVPX_TMP,PVPY_TMP
   REAL(SP) :: PSPX_TMP,PSPY_TMP,PSPXD_TMP,PSPYD_TMP
   REAL(SP) :: U_TMP,V_TMP
   REAL(SP) :: X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2
!------------------------------------------------------------------------------!

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

!
!--Initialize Fluxes-----------------------------------------------------------!
!
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     IA = NIEC(I,1)
     IB = NIEC(I,2)
     IF(IA == NODE_NORTHPOLE)THEN
       XFLUX(IA,K) = 0.0_SP
       XFLUX_ADV(IA,K) = 0.0_SP
     ELSE IF(IB == NODE_NORTHPOLE)THEN  
       XFLUX(IB,K) = 0.0_SP
       XFLUX_ADV(IB,K) = 0.0_SP
     END IF  
   END DO  
     
!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1=NTRG(I)
     DTIJ(I)=DT1(I1)
   END DO

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!


   I = NODE_NORTHPOLE

!   ggao 
   IF(I==0)  RETURN
!  ggao

   PUPX_TMP=0.0_SP
   PUPY_TMP=0.0_SP
   PVPX_TMP=0.0_SP
   PVPY_TMP=0.0_SP

   DO J=1,NTVE(I)
     I1=NBVE(I,J)
     JTMP=NBVT(I,J)
     J1=JTMP+1-(JTMP+1)/4*3
     J2=JTMP+2-(JTMP+2)/4*3
       
     VX_TMP = REARTH * COS(VY(I)*PI/180.0_SP) * COS(VX(I)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(I)*PI/180.0_SP))
     VY_TMP = REARTH * COS(VY(I)*PI/180.0_SP) * SIN(VX(I)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(I)*PI/180.0_SP))
		     
     VX1_TMP= REARTH * COS(VY(NV(I1,J1))*PI/180.0_SP) * COS(VX(NV(I1,J1))*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J1))*PI/180.0_SP))
     VY1_TMP= REARTH * COS(VY(NV(I1,J1))*PI/180.0_SP) * SIN(VX(NV(I1,J1))*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J1))*PI/180.0_SP))
		     
     VX2_TMP= REARTH * COS(YC(I1)*PI/180.0_SP) * COS(XC(I1)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(YC(I1)*PI/180.0_SP))
     VY2_TMP= REARTH * COS(YC(I1)*PI/180.0_SP) * SIN(XC(I1)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(YC(I1)*PI/180.0_SP))
		     
     VX3_TMP= REARTH * COS(VY(NV(I1,J2))*PI/180.0_SP) * COS(VX(NV(I1,J2))*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J2))*PI/180.0_SP))
     VY3_TMP= REARTH * COS(VY(NV(I1,J2))*PI/180.0_SP) * SIN(VX(NV(I1,J2))*PI/180.0_SP) &
                    * 2._SP /(1._SP+SIN(VY(NV(I1,J2))*PI/180.0_SP))
		     
     X11=0.5_SP*(VX_TMP+VX1_TMP)
     Y11=0.5_SP*(VY_TMP+VY1_TMP)
     X22=VX2_TMP
     Y22=VX2_TMP
     X33=0.5_SP*(VX_TMP+VX3_TMP)
     Y33=0.5_SP*(VY_TMP+VY3_TMP)
     
     U_TMP = -V(I1,K)*COS(XC(I1)*PI/180.)-U(I1,K)*SIN(XC(I1)*PI/180.)
     V_TMP = -V(I1,K)*SIN(XC(I1)*PI/180.)+U(I1,K)*COS(XC(I1)*PI/180.)

     PUPX_TMP=PUPX_TMP+U_TMP*(Y11-Y33)
     PUPY_TMP=PUPY_TMP+U_TMP*(X33-X11)
     PVPX_TMP=PVPX_TMP+V_TMP*(Y11-Y33)
     PVPY_TMP=PVPY_TMP+V_TMP*(X33-X11)
   END DO

   PUPX_TMP=PUPX_TMP/ART1(I)
   PUPY_TMP=PUPY_TMP/ART1(I)
   PVPX_TMP=PVPX_TMP/ART1(I)
   PVPY_TMP=PVPY_TMP/ART1(I)
   TMP1=PUPX_TMP**2+PVPY_TMP**2
   TMP2=0.5_SP*(PUPY_TMP+PVPX_TMP)**2
   VISCOFF(I)=SQRT(TMP1+TMP2)*ART1(I)

   IF(K == KBM1) THEN
     AH_BOTTOM(I) = HORCON*(FACT*VISCOFF(I) + FM1)
   END IF

   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1=NTRG(I)
     IA=NIEC(I,1)
     IB=NIEC(I,2)
     
     IF((IA <= M .AND. IB <= M) .AND. I1 <= N)THEN
       XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
       YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))

       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         XI_TMP = REARTH * COS(YI*PI/180.0_SP) * COS(XI*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(YI*PI/180.0_SP))
         YI_TMP = REARTH * COS(YI*PI/180.0_SP) * SIN(XI*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(YI*PI/180.0_SP))

         VXA_TMP = REARTH * COS(VY(IA)*PI/180.0_SP) * COS(VX(IA)*PI/180.0_SP) &
                   * 2._SP /(1._SP+sin(VY(IA)*PI/180.0_SP))
         VYA_TMP = REARTH * COS(VY(IA)*PI/180.0_SP) * SIN(VX(IA)*PI/180.0_SP) &
                   * 2._SP /(1._SP+sin(VY(IA)*PI/180.0_SP))

         VXB_TMP = REARTH * COS(VY(IB)*PI/180.0_SP) * COS(VX(IB)*PI/180.0_SP) &
                   * 2._SP /(1._SP+sin(VY(IB)*PI/180.0_SP))
         VYB_TMP = REARTH * COS(VY(IB)*PI/180.0_SP) * SIN(VX(IB)*PI/180.0_SP) &
                   * 2._SP /(1._SP+sin(VY(IB)*PI/180.0_SP))

!         IF(IA == NODE_NORTHPOLE)THEN
         DXA=XI_TMP-VXA_TMP
         DYA=YI_TMP-VYA_TMP
!         ELSE IF(IB == NODE_NORTHPOLE)THEN
         DXB=XI_TMP-VXB_TMP
         DYB=YI_TMP-VYB_TMP
!	 END IF
!       END IF

        IF(IA == NODE_NORTHPOLE)THEN
	  PSPX_TMP=-PSPY(IB)*COS(VX(IB)*PI/180.)-PSPX(IB)*SIN(VX(IB)*PI/180.)
          PSPY_TMP=-PSPY(IB)*SIN(VX(IB)*PI/180.)+PSPX(IB)*COS(VX(IB)*PI/180.)
   
	  PSPXD_TMP=-PSPYD(IB)*COS(VX(IB)*PI/180.)-PSPXD(IB)*SIN(VX(IB)*PI/180.)
          PSPYD_TMP=-PSPYD(IB)*SIN(VX(IB)*PI/180.)+PSPXD(IB)*COS(VX(IB)*PI/180.)
   
          FIJ1=S1(IA,K)+DXA*PSPX(IA)+DYA*PSPY(IA)
          FIJ2=S1(IB,K)+DXB*PSPX_TMP+DYB*PSPY_TMP

          VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)

          TXX=0.5_SP*(PSPXD(IA)+PSPXD_TMP)*VISCOF
          TYY=0.5_SP*(PSPYD(IA)+PSPYD_TMP)*VISCOF
        ELSE IF(IB == NODE_NORTHPOLE)THEN
	  PSPX_TMP=-PSPY(IA)*COS(VX(IA)*PI/180.)-PSPX(IA)*SIN(VX(IA)*PI/180.)
          PSPY_TMP=-PSPY(IA)*SIN(VX(IA)*PI/180.)+PSPX(IA)*COS(VX(IA)*PI/180.)
   
	  PSPXD_TMP=-PSPYD(IA)*COS(VX(IA)*PI/180.)-PSPXD(IA)*SIN(VX(IA)*PI/180.)
          PSPYD_TMP=-PSPYD(IA)*SIN(VX(IA)*PI/180.)+PSPXD(IA)*COS(VX(IA)*PI/180.)
   
          FIJ1=S1(IA,K)+DXA*PSPX_TMP+DYA*PSPY_TMP
          FIJ2=S1(IB,K)+DXB*PSPX(IB)+DYB*PSPY(IB)

          VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)

          TXX=0.5_SP*(PSPXD_TMP+PSPXD(IB))*VISCOF
          TYY=0.5_SP*(PSPYD_TMP+PSPYD(IB))*VISCOF
        END IF

          FXX=-DTIJ(I)*TXX*DLTYE(I)
       FYY= DTIJ(I)*TYY*DLTXE(I)

!       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         UIJ_TMP = -V(I1,K)*COS(XC(I1)*PI/180.)-U(I1,K)*SIN(XC(I1)*PI/180.)
         VIJ_TMP = -V(I1,K)*SIN(XC(I1)*PI/180.)+U(I1,K)*COS(XC(I1)*PI/180.)
       
         VX1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * COS(XIJE(I,1)*PI/180.0_SP)
         VY1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * SIN(XIJE(I,1)*PI/180.0_SP)

         VX2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * COS(XIJE(I,2)*PI/180.0_SP)
         VY2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * SIN(XIJE(I,2)*PI/180.0_SP)

         DLTXE_TMP = VX2_TMP-VX1_TMP
         DLTYE_TMP = VY2_TMP-VY1_TMP
       
         FXX=-DTIJ(I)*TXX*DLTYE_TMP
         FYY= DTIJ(I)*TYY*DLTXE_TMP

         UVN_TMP = VIJ_TMP*DLTXE_TMP - UIJ_TMP*DLTYE_TMP
         EXFLUX_TMP = -UVN_TMP*DTIJ(I)*((1.0_SP+SIGN(1.0_SP,UVN_TMP))*FIJ2+   &
                      (1.0_SP-SIGN(1.0_SP,UVN_TMP))*FIJ1)*0.5_SP
       
         IF(IA == NODE_NORTHPOLE)THEN
           XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX_TMP+FXX+FYY
           XFLUX_ADV(IA,K)=XFLUX_ADV(IA,K)+EXFLUX_TMP
         ELSE IF(IB == NODE_NORTHPOLE)THEN
           XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX_TMP-FXX-FYY
           XFLUX_ADV(IB,K)=XFLUX_ADV(IB,K)-EXFLUX_TMP
         END IF
       END IF
     END IF  
   END DO

   RETURN
   END SUBROUTINE ADV_S_XY
!==============================================================================|

!==============================================================================|
   SUBROUTINE ADV_T_XY(XFLUX,XFLUX_ADV,PTPX,PTPY,PTPXD,PTPYD,VISCOFF,K)               

!------------------------------------------------------------------------------|

   IMPLICIT NONE
   REAL(SP), DIMENSION(0:MT,KB)     :: XFLUX,XFLUX_ADV
   REAL(SP), DIMENSION(M)           :: PTPX,PTPY,PTPXD,PTPYD,VISCOFF
   REAL(SP), DIMENSION(3*(NT))      :: DTIJ 
   REAL(SP) :: XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
   REAL(SP) :: TXPI,TYPI

   REAL(SP) :: VX_TMP,VY_TMP,VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP,VX3_TMP,VY3_TMP
   REAL(SP) :: XI_TMP,YI_TMP,VXA_TMP,VYA_TMP,VXB_TMP,VYB_TMP
   REAL(SP) :: UIJ_TMP,VIJ_TMP,DLTXE_TMP,DLTYE_TMP,UVN_TMP,EXFLUX_TMP
   REAL(SP) :: PUPX_TMP,PUPY_TMP,PVPX_TMP,PVPY_TMP
   REAL(SP) :: PTPX_TMP,PTPY_TMP,PTPXD_TMP,PTPYD_TMP
   REAL(SP) :: U_TMP,V_TMP
   REAL(SP) :: X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2
!------------------------------------------------------------------------------!

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

!
!--Initialize Fluxes-----------------------------------------------------------!
!
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     IA = NIEC(I,1)
     IB = NIEC(I,2)
     IF(IA == NODE_NORTHPOLE)THEN
       XFLUX(IA,K) = 0.0_SP
       XFLUX_ADV(IA,K) = 0.0_SP
     ELSE IF(IB == NODE_NORTHPOLE)THEN  
       XFLUX(IB,K) = 0.0_SP
       XFLUX_ADV(IB,K) = 0.0_SP
     END IF  
   END DO  
     
!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1=NTRG(I)
     DTIJ(I)=DT1(I1)
   END DO

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!

   I = NODE_NORTHPOLE


!   ggao
   IF(I==0)  RETURN
!  ggao

   PUPX_TMP=0.0_SP
   PUPY_TMP=0.0_SP
   PVPX_TMP=0.0_SP
   PVPY_TMP=0.0_SP

   DO J=1,NTVE(I)
     I1=NBVE(I,J)
     JTMP=NBVT(I,J)
     J1=JTMP+1-(JTMP+1)/4*3
     J2=JTMP+2-(JTMP+2)/4*3
       
     VX_TMP = REARTH * COS(VY(I)*PI/180.0_SP) * COS(VX(I)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(I)*PI/180.0_SP))
     VY_TMP = REARTH * COS(VY(I)*PI/180.0_SP) * SIN(VX(I)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(I)*PI/180.0_SP))
		     
     VX1_TMP= REARTH * COS(VY(NV(I1,J1))*PI/180.0_SP) * COS(VX(NV(I1,J1))*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J1))*PI/180.0_SP))
     VY1_TMP= REARTH * COS(VY(NV(I1,J1))*PI/180.0_SP) * SIN(VX(NV(I1,J1))*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J1))*PI/180.0_SP))
		     
     VX2_TMP= REARTH * COS(YC(I1)*PI/180.0_SP) * COS(XC(I1)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(YC(I1)*PI/180.0_SP))
     VY2_TMP= REARTH * COS(YC(I1)*PI/180.0_SP) * SIN(XC(I1)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(YC(I1)*PI/180.0_SP))
		     
     VX3_TMP= REARTH * COS(VY(NV(I1,J2))*PI/180.0_SP) * COS(VX(NV(I1,J2))*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J2))*PI/180.0_SP))
     VY3_TMP= REARTH * COS(VY(NV(I1,J2))*PI/180.0_SP) * SIN(VX(NV(I1,J2))*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J2))*PI/180.0_SP))
		     
     X11=0.5_SP*(VX_TMP+VX1_TMP)
     Y11=0.5_SP*(VY_TMP+VY1_TMP)
     X22=VX2_TMP
     Y22=VX2_TMP
     X33=0.5_SP*(VX_TMP+VX3_TMP)
     Y33=0.5_SP*(VY_TMP+VY3_TMP)
     
     U_TMP = -V(I1,K)*COS(XC(I1)*PI/180.)-U(I1,K)*SIN(XC(I1)*PI/180.)
     V_TMP = -V(I1,K)*SIN(XC(I1)*PI/180.)+U(I1,K)*COS(XC(I1)*PI/180.)

     PUPX_TMP=PUPX_TMP+U_TMP*(Y11-Y33)
     PUPY_TMP=PUPY_TMP+U_TMP*(X33-X11)
     PVPX_TMP=PVPX_TMP+V_TMP*(Y11-Y33)
     PVPY_TMP=PVPY_TMP+V_TMP*(X33-X11)
   END DO

   PUPX_TMP=PUPX_TMP/ART1(I)
   PUPY_TMP=PUPY_TMP/ART1(I)
   PVPX_TMP=PVPX_TMP/ART1(I)
   PVPY_TMP=PVPY_TMP/ART1(I)
   TMP1=PUPX_TMP**2+PVPY_TMP**2
   TMP2=0.5_SP*(PUPY_TMP+PVPX_TMP)**2
   VISCOFF(I)=SQRT(TMP1+TMP2)*ART1(I)

   IF(K == KBM1) THEN
     AH_BOTTOM(I) = HORCON*(FACT*VISCOFF(I) + FM1)
   END IF

   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1=NTRG(I)
     IA=NIEC(I,1)
     IB=NIEC(I,2)
     IF(IA <= M .AND. IB <= M .AND. I1 <= N)THEN
       XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
       YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))

       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         XI_TMP = REARTH * COS(YI*PI/180.0_SP) * COS(XI*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(YI*PI/180.0_SP))
         YI_TMP = REARTH * COS(YI*PI/180.0_SP) * SIN(XI*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(YI*PI/180.0_SP))

         VXA_TMP = REARTH * COS(VY(IA)*PI/180.0_SP) * COS(VX(IA)*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(VY(IA)*PI/180.0_SP))
         VYA_TMP = REARTH * COS(VY(IA)*PI/180.0_SP) * SIN(VX(IA)*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(VY(IA)*PI/180.0_SP))

         VXB_TMP = REARTH * COS(VY(IB)*PI/180.0_SP) * COS(VX(IB)*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(VY(IB)*PI/180.0_SP))
         VYB_TMP = REARTH * COS(VY(IB)*PI/180.0_SP) * SIN(VX(IB)*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(VY(IB)*PI/180.0_SP))

!        IF(IA == NODE_NORTHPOLE)THEN
           DXA=XI_TMP-VXA_TMP
           DYA=YI_TMP-VYA_TMP
!	 ELSE IF(IB == NODE_NORTHPOLE)THEN
           DXB=XI_TMP-VXB_TMP
           DYB=YI_TMP-VYB_TMP
!	 END IF
!       END IF

        IF(IA == NODE_NORTHPOLE)THEN
	  PTPX_TMP=-PTPY(IB)*COS(VX(IB)*PI/180.)-PTPX(IB)*SIN(VX(IB)*PI/180.)
          PTPY_TMP=-PTPY(IB)*SIN(VX(IB)*PI/180.)+PTPX(IB)*COS(VX(IB)*PI/180.)
   
	  PTPXD_TMP=-PTPYD(IB)*COS(VX(IB)*PI/180.)-PTPXD(IB)*SIN(VX(IB)*PI/180.)
          PTPYD_TMP=-PTPYD(IB)*SIN(VX(IB)*PI/180.)+PTPXD(IB)*COS(VX(IB)*PI/180.)
   
          FIJ1=T1(IA,K)+DXA*PTPX(IA)+DYA*PTPY(IA)
          FIJ2=T1(IB,K)+DXB*PTPX_TMP+DYB*PTPY_TMP

          VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)

          TXX=0.5_SP*(PTPXD(IA)+PTPXD_TMP)*VISCOF
          TYY=0.5_SP*(PTPYD(IA)+PTPYD_TMP)*VISCOF
        ELSE IF(IB == NODE_NORTHPOLE)THEN
	  PTPX_TMP=-PTPY(IA)*COS(VX(IA)*PI/180.)-PTPX(IA)*SIN(VX(IA)*PI/180.)
          PTPY_TMP=-PTPY(IA)*SIN(VX(IA)*PI/180.)+PTPX(IA)*COS(VX(IA)*PI/180.)
   
	  PTPXD_TMP=-PTPYD(IA)*COS(VX(IA)*PI/180.)-PTPXD(IA)*SIN(VX(IA)*PI/180.)
          PTPYD_TMP=-PTPYD(IA)*SIN(VX(IA)*PI/180.)+PTPXD(IA)*COS(VX(IA)*PI/180.)
   
          FIJ1=T1(IA,K)+DXA*PTPX_TMP+DYA*PTPY_TMP
          FIJ2=T1(IB,K)+DXB*PTPX(IB)+DYB*PTPY(IB)

          VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)

          TXX=0.5_SP*(PTPXD_TMP+PTPXD(IB))*VISCOF
          TYY=0.5_SP*(PTPYD_TMP+PTPYD(IB))*VISCOF
        END IF

!     IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         UIJ_TMP = -V(I1,K)*COS(XC(I1)*PI/180.)-U(I1,K)*SIN(XC(I1)*PI/180.)
         VIJ_TMP = -V(I1,K)*SIN(XC(I1)*PI/180.)+U(I1,K)*COS(XC(I1)*PI/180.)
       
         VX1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * COS(XIJE(I,1)*PI/180.0_SP)
         VY1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * SIN(XIJE(I,1)*PI/180.0_SP)

         VX2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * COS(XIJE(I,2)*PI/180.0_SP)
         VY2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * SIN(XIJE(I,2)*PI/180.0_SP)

         DLTXE_TMP = VX2_TMP-VX1_TMP
         DLTYE_TMP = VY2_TMP-VY1_TMP
       
         FXX=-DTIJ(I)*TXX*DLTYE_TMP
         FYY= DTIJ(I)*TYY*DLTXE_TMP

         UVN_TMP = VIJ_TMP*DLTXE_TMP - UIJ_TMP*DLTYE_TMP
         EXFLUX_TMP = -UVN_TMP*DTIJ(I)*((1.0_SP+SIGN(1.0_SP,UVN_TMP))*FIJ2+   &
                      (1.0_SP-SIGN(1.0_SP,UVN_TMP))*FIJ1)*0.5_SP
       
         IF(IA == NODE_NORTHPOLE)THEN
           XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX_TMP+FXX+FYY
           XFLUX_ADV(IA,K)=XFLUX_ADV(IA,K)+EXFLUX_TMP
         ELSE IF(IB == NODE_NORTHPOLE)THEN
           XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX_TMP-FXX-FYY
           XFLUX_ADV(IB,K)=XFLUX_ADV(IB,K)-EXFLUX_TMP
         END IF
       END IF
     END IF  
   END DO 

   RETURN
   END SUBROUTINE ADV_T_XY
!==============================================================================|

!==============================================================================|
!     CALCULATE THE BAROCLINIC PRESSURE GRADIENT IN SIGMA COORDINATES          |
!==============================================================================|

   SUBROUTINE BAROPG_XY(DRIJK1,DRIJK2) 

!==============================================================================|
   IMPLICIT NONE
   REAL(SP) :: DRIJK1(0:N,3,KBM1), DRIJK2(0:N,KBM1)
   REAL(SP) :: DIJ,DRHO1,DRHO2
   INTEGER  :: I,II,K,J,J1,J2,IJK
   REAL(SP) :: VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP
!==============================================================================|

   DO II = 1, NP
     I=NP_LST(II)
     DO K=1,KBM1
       DRHOX(I,K)=0.0_SP
       DRHOY(I,K)=0.0_SP
     END DO
   END DO

   DO II = 1, NP
     I=NP_LST(II)
     DO K=1,KBM1
        DO J = 1, 3
          J1=J+1-INT((J+1)/4)*3
          J2=J+2-INT((J+2)/4)*3
          IJK=NBE(I,J)
          DIJ=0.5_SP*(DT(NV(I,J1))+DT(NV(I,J2)))

          VY1_TMP=REARTH*COS(VY(NV(I,J1))*PI/180.0_SP)*SIN(VX(NV(I,J1))*PI/180.0_SP)
          VY2_TMP=REARTH*COS(VY(NV(I,J2))*PI/180.0_SP)*SIN(VX(NV(I,J2))*PI/180.0_SP)

          DRHO1=(VY1_TMP-VY2_TMP)*DRIJK1(I,J,K)*DT1(I)
          DRHO2=(VY1_TMP-VY2_TMP)*DIJ*DRIJK2(I,K)

          DRHOX(I,K)=DRHOX(I,K)+DRHO1+DRHO2

          VX1_TMP=REARTH*COS(VY(NV(I,J1))*PI/180.0_SP)*COS(VX(NV(I,J1))*PI/180.0_SP)
          VX2_TMP=REARTH*COS(VY(NV(I,J2))*PI/180.0_SP)*COS(VX(NV(I,J2))*PI/180.0_SP)

	  DRHO1=(VX2_TMP-VX1_TMP)*DRIJK1(I,J,K)*DT1(I)
          DRHO2=(VX2_TMP-VX1_TMP)*DIJ*DRIJK2(I,K)

          DRHOY(I,K)=DRHOY(I,K)+DRHO1+DRHO2

       END DO
     END DO
   END DO

   RETURN
   END SUBROUTINE BAROPG_XY
!==============================================================================|

!==============================================================================|
   SUBROUTINE SHAPE_COEF_XY

!----------------------------------------------------------------------!
!  This subrountine is used to calculate the coefficient for a linear  !
!  function on the x-y plane, i.e.:                                    !
!                     r(x,y;phai)=phai_c+cofa1*x+cofa2*y               !
!     innc(i)=0    cells on the boundary                               !
!     innc(i)=1    cells in the interior                               !
!----------------------------------------------------------------------!
     
   USE ALL_VARS
   IMPLICIT NONE
   REAL(DP) X1,X2,X3,Y1,Y2,Y3,DELT,AI1,AI2,AI3,BI1,BI2,BI3,CI1,CI2,CI3
   REAL(DP) DELTX,DELTY,TEMP1,ANG1,ANG2,B1,B2,ANGLE
   REAL(DP), ALLOCATABLE :: XC_TMP(:),YC_TMP(:),VX_TMP(:),VY_TMP(:)
   INTEGER  I,II,J,JJ,J1,J2
!
!---------------interior cells-----------------------------------------!
!
   ALLOCATE(A1U_XY(N,4)); A1U_XY = 0.0_SP
   ALLOCATE(A2U_XY(N,4)); A2U_XY = 0.0_SP
   ALLOCATE(AW0_XY(N,3)); AW0_XY = 0.0_SP
   ALLOCATE(AWX_XY(N,3)); AWX_XY = 0.0_SP
   ALLOCATE(AWY_XY(N,3)); AWY_XY = 0.0_SP
   
   ALLOCATE(XC_TMP(0:NT)); XC_TMP = 0.0_SP
   ALLOCATE(YC_TMP(0:NT)); YC_TMP = 0.0_SP
   ALLOCATE(VX_TMP(0:MT)); VX_TMP = 0.0_SP
   ALLOCATE(VY_TMP(0:MT)); VY_TMP = 0.0_SP
   
   DO I=1,NT
     XC_TMP(I) = REARTH * COS(YC(I)*PI/180.0_SP) * COS(XC(I)*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(YC(I)*PI/180.0_SP))
     YC_TMP(I) = REARTH * COS(YC(I)*PI/180.0_SP) * SIN(XC(I)*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(YC(I)*PI/180.0_SP))
   END DO		  

   DO I=1,MT
     VX_TMP(I) = REARTH * COS(VY(I)*PI/180.0_SP) * COS(VX(I)*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(VY(I)*PI/180.0_SP))
     VY_TMP(I) = REARTH * COS(VY(I)*PI/180.0_SP) * SIN(VX(I)*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(VY(I)*PI/180.0_SP))
   END DO		  

   DO I=1,N
     IF(ISBCE(I) == 0)THEN
       Y1 = YC_TMP(NBE(I,1))-YC_TMP(I)
       Y2 = YC_TMP(NBE(I,2))-YC_TMP(I)
       Y3 = YC_TMP(NBE(I,3))-YC_TMP(I)
       X1=XC_TMP(NBE(I,1))-XC_TMP(I)
       X2=XC_TMP(NBE(I,2))-XC_TMP(I)
       X3=XC_TMP(NBE(I,3))-XC_TMP(I)

       X1=X1/1000.0_SP
       X2=X2/1000.0_SP
       X3=X3/1000.0_SP
       Y1=Y1/1000.0_SP
       Y2=Y2/1000.0_SP
       Y3=Y3/1000.0_SP

       delt=(x1*y2-x2*y1)**2+(x1*y3-x3*y1)**2+(x2*y3-x3*y2)**2
       delt=delt*1000.

       a1u_XY(i,1)=(y1+y2+y3)*(x1*y1+x2*y2+x3*y3)- &
                (x1+x2+x3)*(y1**2+y2**2+y3**2)
       a1u_XY(i,1)=a1u_XY(i,1)/delt
       a1u_XY(i,2)=(y1**2+y2**2+y3**2)*x1-(x1*y1+x2*y2+x3*y3)*y1
       a1u_XY(i,2)=a1u_XY(i,2)/delt
       a1u_XY(i,3)=(y1**2+y2**2+y3**2)*x2-(x1*y1+x2*y2+x3*y3)*y2
       a1u_XY(i,3)=a1u_XY(i,3)/delt
       a1u_XY(i,4)=(y1**2+y2**2+y3**2)*x3-(x1*y1+x2*y2+x3*y3)*y3
       a1u_XY(i,4)=a1u_XY(i,4)/delt

       a2u_XY(i,1)=(x1+x2+x3)*(x1*y1+x2*y2+x3*y3)- &
                (y1+y2+y3)*(x1**2+x2**2+x3**2)
       a2u_XY(i,1)=a2u_XY(i,1)/delt
       a2u_XY(i,2)=(x1**2+x2**2+x3**2)*y1-(x1*y1+x2*y2+x3*y3)*x1
       a2u_XY(i,2)=a2u_XY(i,2)/delt
       a2u_XY(i,3)=(x1**2+x2**2+x3**2)*y2-(x1*y1+x2*y2+x3*y3)*x2
       a2u_XY(i,3)=a2u_XY(i,3)/delt
       a2u_XY(i,4)=(x1**2+x2**2+x3**2)*y3-(x1*y1+x2*y2+x3*y3)*x3
       a2u_XY(i,4)=a2u_XY(i,4)/delt
     end if

     x1=vx_TMP(nv(i,1))-xc_TMP(i)
     x2=vx_TMP(nv(i,2))-xc_TMP(i)
     x3=vx_TMP(nv(i,3))-xc_TMP(i)
     y1=vy_TMP(nv(i,1))-yc_TMP(i)
     y2=vy_TMP(nv(i,2))-yc_TMP(i)
     y3=vy_TMP(nv(i,3))-yc_TMP(i)


     ai1=y2-y3
     ai2=y3-y1
     ai3=y1-y2
     bi1=x3-x2
     bi2=x1-x3
     bi3=x2-x1
     ci1=x2*y3-x3*y2
     ci2=x3*y1-x1*y3
     ci3=x1*y2-x2*y1

     aw0_XY(i,1)=-ci1/2./art(i)
     aw0_XY(i,2)=-ci2/2./art(i)
     aw0_XY(i,3)=-ci3/2./art(i)
     awx_XY(i,1)=-ai1/2./art(i)
     awx_XY(i,2)=-ai2/2./art(i)
     awx_XY(i,3)=-ai3/2./art(i)
     awy_XY(i,1)=-bi1/2./art(i)
     awy_XY(i,2)=-bi2/2./art(i)
     awy_XY(i,3)=-bi3/2./art(i)
   end do

   DEALLOCATE(XC_TMP,YC_TMP,VX_TMP,VY_TMP)
   
   return
   end subroutine shape_coef_xy

!==============================================================================|

!==============================================================================|
!==============================================================================|
   SUBROUTINE ADV_Q_XY(XFLUX,PQPX,PQPY,PQPXD,PQPYD,VISCOFF,Q,UQ,VQ,K)               

!==============================================================================|
!   Calculate the Turbulent Kinetic Energy and Mixing Length Based on          |
!   The Mellor-Yamada Level 2.5 Turbulent Closure Model                        |
!==============================================================================|


!------------------------------------------------------------------------------|

   IMPLICIT NONE
   REAL(SP), DIMENSION(0:MT,KB)     :: XFLUX,Q
   REAL(SP), DIMENSION(M)           :: PQPX,PQPY,PQPXD,PQPYD,VISCOFF
   REAL(SP), DIMENSION(3*(NT))      :: DTIJ 
   REAL(SP) :: XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2 
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF   
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
   REAL(SP) :: TXPI,TYPI

   REAL(SP) :: VX_TMP,VY_TMP,VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP,VX3_TMP,VY3_TMP
   REAL(SP) :: XI_TMP,YI_TMP,VXA_TMP,VYA_TMP,VXB_TMP,VYB_TMP
   REAL(SP) :: UIJ_TMP,VIJ_TMP,DLTXE_TMP,DLTYE_TMP,UVN_TMP,EXFLUX_TMP
   REAL(SP) :: PUPX_TMP,PUPY_TMP,PVPX_TMP,PVPY_TMP
   REAL(SP) :: PQPX_TMP,PQPY_TMP,PQPXD_TMP,PQPYD_TMP
   REAL(SP) :: U_TMP,V_TMP
   REAL(SP) :: X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2

   REAL(SP), DIMENSION(0:NT,KB)    :: UQ,VQ
!------------------------------------------------------------------------------!

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

!
!--Initialize Fluxes-----------------------------------------------------------!
!
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     IA = NIEC(I,1)
     IB = NIEC(I,2)
     IF(IA == NODE_NORTHPOLE)THEN
       XFLUX(IA,K) = 0.0_SP
     ELSE IF(IB == NODE_NORTHPOLE)THEN  
       XFLUX(IB,K) = 0.0_SP
     END IF  
   END DO  
     
!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1=NTRG(I)
     DTIJ(I)=DT1(I1)
   END DO

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!


   I = NODE_NORTHPOLE

!   ggao
   IF(I==0)  RETURN
!  ggao

   if(i /= 0)then
   PUPX_TMP=0.0_SP
   PUPY_TMP=0.0_SP
   PVPX_TMP=0.0_SP
   PVPY_TMP=0.0_SP

   DO J=1,NTVE(I)
     I1=NBVE(I,J)
     JTMP=NBVT(I,J)
     J1=JTMP+1-(JTMP+1)/4*3
     J2=JTMP+2-(JTMP+2)/4*3
      
     VX_TMP = REARTH * COS(VY(I)*PI/180.0_SP) * COS(VX(I)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(I)*PI/180.0_SP))
     VY_TMP = REARTH * COS(VY(I)*PI/180.0_SP) * SIN(VX(I)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(I)*PI/180.0_SP))
		     
     VX1_TMP= REARTH * COS(VY(NV(I1,J1))*PI/180.0_SP) * COS(VX(NV(I1,J1))*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J1))*PI/180.0_SP))
     VY1_TMP= REARTH * COS(VY(NV(I1,J1))*PI/180.0_SP) * SIN(VX(NV(I1,J1))*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J1))*PI/180.0_SP))
		     
     VX2_TMP= REARTH * COS(YC(I1)*PI/180.0_SP) * COS(XC(I1)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(YC(I1)*PI/180.0_SP))
     VY2_TMP= REARTH * COS(YC(I1)*PI/180.0_SP) * SIN(XC(I1)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(YC(I1)*PI/180.0_SP))
		     
     VX3_TMP= REARTH * COS(VY(NV(I1,J2))*PI/180.0_SP) * COS(VX(NV(I1,J2))*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J2))*PI/180.0_SP))
     VY3_TMP= REARTH * COS(VY(NV(I1,J2))*PI/180.0_SP) * SIN(VX(NV(I1,J2))*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J2))*PI/180.0_SP))
		     
     X11=0.5_SP*(VX_TMP+VX1_TMP)
     Y11=0.5_SP*(VY_TMP+VY1_TMP)
     X22=VX2_TMP
     Y22=VX2_TMP
     X33=0.5_SP*(VX_TMP+VX3_TMP)
     Y33=0.5_SP*(VY_TMP+VY3_TMP)
     
     U_TMP = -VQ(I1,K)*COS(XC(I1)*PI/180.)-UQ(I1,K)*SIN(XC(I1)*PI/180.)
     V_TMP = -VQ(I1,K)*SIN(XC(I1)*PI/180.)+UQ(I1,K)*COS(XC(I1)*PI/180.)

     PUPX_TMP=PUPX_TMP+U_TMP*(Y11-Y33)
     PUPY_TMP=PUPY_TMP+U_TMP*(X33-X11)
     PVPX_TMP=PVPX_TMP+V_TMP*(Y11-Y33)
     PVPY_TMP=PVPY_TMP+V_TMP*(X33-X11)
   END DO

   PUPX_TMP=PUPX_TMP/ART1(I)
   PUPY_TMP=PUPY_TMP/ART1(I)
   PVPX_TMP=PVPX_TMP/ART1(I)
   PVPY_TMP=PVPY_TMP/ART1(I)
   TMP1=PUPX_TMP**2+PVPY_TMP**2
   TMP2=0.5_SP*(PUPY_TMP+PVPX_TMP)**2
   VISCOFF(I)=SQRT(TMP1+TMP2)*ART1(I)

   IF(K == KBM1) THEN
     AH_BOTTOM(I) = HORCON*(FACT*VISCOFF(I) + FM1)/HPRNU
   END IF
   endif
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1=NTRG(I)
     IA=NIEC(I,1)
     IB=NIEC(I,2)
     
     IF((IA <= M .AND. IB <= M) .AND. I1 <= N)THEN
       XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
       YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))

       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         XI_TMP = REARTH * COS(YI*PI/180.0_SP) * COS(XI*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(YI*PI/180.0_SP))
         YI_TMP = REARTH * COS(YI*PI/180.0_SP) * SIN(XI*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(YI*PI/180.0_SP))

         VXA_TMP = REARTH * COS(VY(IA)*PI/180.0_SP) * COS(VX(IA)*PI/180.0_SP) &
                   * 2._SP /(1._SP+sin(VY(IA)*PI/180.0_SP))
         VYA_TMP = REARTH * COS(VY(IA)*PI/180.0_SP) * SIN(VX(IA)*PI/180.0_SP) &
                   * 2._SP /(1._SP+sin(VY(IA)*PI/180.0_SP))

         VXB_TMP = REARTH * COS(VY(IB)*PI/180.0_SP) * COS(VX(IB)*PI/180.0_SP) &
                   * 2._SP /(1._SP+sin(VY(IB)*PI/180.0_SP))
         VYB_TMP = REARTH * COS(VY(IB)*PI/180.0_SP) * SIN(VX(IB)*PI/180.0_SP) &
                   * 2._SP /(1._SP+sin(VY(IB)*PI/180.0_SP))

!         IF(IA == NODE_NORTHPOLE)THEN
         DXA=XI_TMP-VXA_TMP
         DYA=YI_TMP-VYA_TMP
!         ELSE IF(IB == NODE_NORTHPOLE)THEN
         DXB=XI_TMP-VXB_TMP
         DYB=YI_TMP-VYB_TMP
!	 END IF
!       END IF

        IF(IA == NODE_NORTHPOLE)THEN
	  PQPX_TMP=-PQPY(IB)*COS(VX(IB)*PI/180.)-PQPX(IB)*SIN(VX(IB)*PI/180.)
          PQPY_TMP=-PQPY(IB)*SIN(VX(IB)*PI/180.)+PQPX(IB)*COS(VX(IB)*PI/180.)
   
	  PQPXD_TMP=-PQPYD(IB)*COS(VX(IB)*PI/180.)-PQPXD(IB)*SIN(VX(IB)*PI/180.)
          PQPYD_TMP=-PQPYD(IB)*SIN(VX(IB)*PI/180.)+PQPXD(IB)*COS(VX(IB)*PI/180.)
   
          FIJ1=Q(IA,K)+DXA*PQPX(IA)+DYA*PQPY(IA)
          FIJ2=Q(IB,K)+DXB*PQPX_TMP+DYB*PQPY_TMP

          VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)

          TXX=0.5_SP*(PQPXD(IA)+PQPXD_TMP)*VISCOF
          TYY=0.5_SP*(PQPYD(IA)+PQPYD_TMP)*VISCOF
        ELSE IF(IB == NODE_NORTHPOLE)THEN
	  PQPX_TMP=-PQPY(IA)*COS(VX(IA)*PI/180.)-PQPX(IA)*SIN(VX(IA)*PI/180.)
          PQPY_TMP=-PQPY(IA)*SIN(VX(IA)*PI/180.)+PQPX(IA)*COS(VX(IA)*PI/180.)
   
	  PQPXD_TMP=-PQPYD(IA)*COS(VX(IA)*PI/180.)-PQPXD(IA)*SIN(VX(IA)*PI/180.)
          PQPYD_TMP=-PQPYD(IA)*SIN(VX(IA)*PI/180.)+PQPXD(IA)*COS(VX(IA)*PI/180.)
   
          FIJ1=Q(IA,K)+DXA*PQPX_TMP+DYA*PQPY_TMP
          FIJ2=Q(IB,K)+DXB*PQPX(IB)+DYB*PQPY(IB)

          VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)/HPRNU

          TXX=0.5_SP*(PQPXD_TMP+PQPXD(IB))*VISCOF
          TYY=0.5_SP*(PQPYD_TMP+PQPYD(IB))*VISCOF
        END IF

          FXX=-DTIJ(I)*TXX*DLTYE(I)
       FYY= DTIJ(I)*TYY*DLTXE(I)

!       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         UIJ_TMP = -VQ(I1,K)*COS(XC(I1)*PI/180.)-UQ(I1,K)*SIN(XC(I1)*PI/180.)
         VIJ_TMP = -VQ(I1,K)*SIN(XC(I1)*PI/180.)+UQ(I1,K)*COS(XC(I1)*PI/180.)
       
         VX1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * COS(XIJE(I,1)*PI/180.0_SP)
         VY1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * SIN(XIJE(I,1)*PI/180.0_SP)

         VX2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * COS(XIJE(I,2)*PI/180.0_SP)
         VY2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * SIN(XIJE(I,2)*PI/180.0_SP)

         DLTXE_TMP = VX2_TMP-VX1_TMP
         DLTYE_TMP = VY2_TMP-VY1_TMP
       
         FXX=-DTIJ(I)*TXX*DLTYE_TMP
         FYY= DTIJ(I)*TYY*DLTXE_TMP

         UVN_TMP = VIJ_TMP*DLTXE_TMP - UIJ_TMP*DLTYE_TMP
         EXFLUX_TMP = -UVN_TMP*DTIJ(I)*((1.0_SP+SIGN(1.0_SP,UVN_TMP))*FIJ2+   &
                      (1.0_SP-SIGN(1.0_SP,UVN_TMP))*FIJ1)*0.5_SP
       
         IF(IA == NODE_NORTHPOLE)THEN
           XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX_TMP+FXX+FYY
         ELSE IF(IB == NODE_NORTHPOLE)THEN
           XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX_TMP-FXX-FYY
         END IF
       END IF
     END IF  
   END DO

   RETURN
   END SUBROUTINE ADV_Q_XY
!==============================================================================|

END MODULE MOD_NORTHPOLE
