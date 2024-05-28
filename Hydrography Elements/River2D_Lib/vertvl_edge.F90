!==============================================================================|
!   CALCULATE THE SIGMA COORDINATE VERTICAL VELOCITY FOR THE 3D MODE (omega)   |
!							                       |
!   DETERMINED FROM EQUATION:						       |
!   									       !
!   d/dt(D) + d/dx(uD) + d/dy(uD) = d/sigma(omega)                             !
!==============================================================================|

   SUBROUTINE VERTVL_EDGE         

!------------------------------------------------------------------------------|
   USE ALL_VARS
   USE BCS
   USE MOD_WD
   USE MOD_NORTHPOLE
 
   USE MOD_MEANFLOW

   IMPLICIT NONE 
   REAL(SP) :: XFLUX(MT,KBM1)
   REAL(SP) :: DIJ,UIJ,VIJ,UN,EXFLUX,TMP1
   INTEGER  :: I,K,IA,IB,I1 ,J,JJ,J1,J2
!------------------------------------------------------------------------------|

!----------------------INITIALIZE FLUX-----------------------------------------!

   XFLUX = 0.0_SP

!----------------------ACCUMULATE FLUX-----------------------------------------!

  DO I=1,NCV
     I1=NTRG(I)
     IA=NIEC(I,1)
     IB=NIEC(I,2)
     DIJ=DT1(I1)
     DO K=1,KBM1
       UIJ=US(I1,K)
       VIJ=VS(I1,K)
       EXFLUX=DIJ*(-UIJ*DLTYE(I)+VIJ*DLTXE(I))
       XFLUX(IA,K)=XFLUX(IA,K)-EXFLUX
       XFLUX(IB,K)=XFLUX(IB,K)+EXFLUX
     END DO
   END DO
  

  if  (SPHERICAL==1) then
  if  (NORTHPOLE==1) then 
   CALL VERTVL_EDGE_XY(XFLUX)
  endif
  endif
   
!-----------------------NULLIFY BOUNDARY FLUX----------------------------------!
! For "tide + meanflow"/"meanflow only" case, this part should be commented out;
! For "tide only" case, this part may be kept.
! However, the effect of this term is small from my experience.
if  (MEAN_FLOW==0)then
      DO I=1,M
        DO K=1,KBM1
          IF(ISONB(I) == 2) XFLUX(I,K)=0.0_SP  
        ENDDO
      ENDDO
! can be changed to (no IF statements)
!     DO I=1,IOBCN
!        DO K=1,KBM1
!           XFLUX(I_OBC_N(I),K)=0.0_SP
!        ENDDO
!     ENDDO
endif
!-----------------------FRESH WATER INFLOW-------------------------------------!

   IF(NUMQBC >= 1) THEN
     IF(INFLOW_TYPE == 'node') THEN
       DO J=1,NUMQBC
         JJ=INODEQ(J)
         DO K=1,KBM1
           XFLUX(JJ,K)=XFLUX(JJ,K)-QDIS(J)*VQDIST(J,K)/DZ(K)
         END DO
       END DO
     ELSE IF(INFLOW_TYPE == 'edge') THEN
       DO J=1,NUMQBC
         J1=N_ICELLQ(J,1)
         J2=N_ICELLQ(J,2)
         DO K=1,KBM1
           XFLUX(J1,K)=XFLUX(J1,K)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)/DZ(K)
           XFLUX(J2,K)=XFLUX(J2,K)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)/DZ(K)
         END DO
       END DO
     END IF
   END IF

  if  (MEAN_FLOW==1)then
   IF (nmfcell > 0) THEN
     DO I = 1, nmfcell
        J1= NODE_MFCELL(I,1)
        J2= NODE_MFCELL(I,2)
        DO K=1,KBM1
           XFLUX(J1,K) = XFLUX(J1,K) - MFQDIS(I)*RDISMF(I,1)*MFDIST(I,K)/DZ(K)
           XFLUX(J2,K) = XFLUX(J2,K) - MFQDIS(I)*RDISMF(I,2)*MFDIST(I,K)/DZ(K)
        END DO
     END DO
   END IF
  endif

!---IF NO FRESH WATER INFLOW, OMEGA IS ZERO AT FREE SURFACE AND BOTTOM---------!

   DO I=1,M
     WTS(I,1)=0.0_SP
     WTS(I,KB)=0.0_SP
   ENDDO


!--------------------------CALCULATE OMEGA-------------------------------------!

   DO I=1,M
    IF(ISWETNT(I)*ISWETN(I) == 1)THEN
     DO K=1,KBM1
       WTS(I,K+1)=WTS(I,K)+DZ(K)*(XFLUX(I,K)/ART1(I)+(EL(I)-ET(I))/DTI)
     END DO
    ELSE
     DO K=1,KBM1
       WTS(I,K+1)=0.0_SP
     END DO
    END IF
   END DO

!--------------------------ADJUST OMEGA----------------------------------------!
! IMPROVES MASS CONSERVATION

   DO I=1,M
     IF(ABS(WTS(I,KB)) > 1.0E-8_SP)THEN
!       IF(ISBCE(I) /= 2)THEN
       IF(mean_flow==0.and.ISONB(I) /= 2)THEN
! The effect of this comment is also small, at least from the experiment I did
         TMP1=ELF(I)*FLOAT(KBM1)-WTS(I,KB)*DTI/DZ(1)
         TMP1=TMP1/FLOAT(KBM1)
         DTFA(I)=TMP1+H(I)
         DO K=2,KB
           WTS(I,K)=WTS(I,K)-FLOAT(K-1)/FLOAT(KBM1)*WTS(I,KB)
         END DO
       END IF
     END IF
   END DO
!
!----TRANSFER OMEGA TO FACE CENTER---------------------------------------------!
!
   DO I=1,N
     DO K=1,KB
       W(I,K) = ONE_THIRD*(WTS(NV(I,1),K)+WTS(NV(I,2),K)+WTS(NV(I,3),K))
     END DO
   END DO

   DO I=1,N
     DO K=1,KB
       W(I,K) = FLOAT(ISWETC(I))*W(I,K)
     END DO
   END DO

   RETURN
   END SUBROUTINE VERTVL_EDGE
!==============================================================================|
