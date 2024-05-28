!==============================================================================!

   SUBROUTINE ADJUST2D3D(ADJUST_TYPE)
!==============================================================================|
!    ADJUST 3D VELOCITY USING DEFECT BETWEEN UPDATED AND CURRENT VERTICALLY    !
!    AVERAGED VELOCITIES						       !
! 									       !
!    FORMULA IS:							       !
!									       !
!      U_adjusted = U_orig + eps*(U_avg_new - U_avg_current)		       !
!      eps = 0 : no adjustment						       !
!      eps = 1 : full adjustment					       !
!==============================================================================|
   USE MOD_WD
   USE ALL_VARS
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: ADJUST_TYPE
   INTEGER :: I,K
   REAL(SP), PARAMETER :: EPS = 1.0_SP
   REAL(SP) :: UAC,VAC,UTMP,VTMP
!==============================================================================!


   SELECT CASE(ADJUST_TYPE)

   CASE(1)
   DO I=1,NT
     UAC    = SUM(U(I,1:KBM1)*DZ(1:KBM1))
     VAC    = SUM(V(I,1:KBM1)*DZ(1:KBM1))
     U(I,1:KBM1) = U(I,1:KBM1) + EPS*(UA(I) - UAC) 
     V(I,1:KBM1) = V(I,1:KBM1) + EPS*(VA(I) - VAC) 
   END DO

   CASE(2)
   UARD = UARD/FLOAT(ISPLIT)
   VARD = VARD/FLOAT(ISPLIT)
   UARDS = UARDS/FLOAT(ISPLIT)
   VARDS = VARDS/FLOAT(ISPLIT)

                                                                                                                         
   DO I=1,NT
     IF(ISWETCT(I)*ISWETC(I) == 1)THEN
       UTMP = 0.0_SP ; VTMP = 0.0_SP
       DO K=1,KBM1
         UTMP = UTMP + U(I,K)*DZ(K)
         VTMP = VTMP + V(I,K)*DZ(K)
       END DO
       UTMP = UTMP*DT1(I)
       VTMP = VTMP*DT1(I)
       DO K=1,KBM1
         U(I,K) = U(I,K) - (UTMP-UARD(I))/DT1(I)
         V(I,K) = V(I,K) - (VTMP-VARD(I))/DT1(I)
       END DO
     END IF
   END DO

   DO I=1,NT
     UTMP = 0.0_SP ; VTMP = 0.0_SP
     DO K=1,KBM1
       UTMP = UTMP + U(I,K)*DZ(K)
       VTMP = VTMP + V(I,K)*DZ(K)
     END DO
     UTMP = UTMP*DT1(I)
     VTMP = VTMP*DT1(I)
     DO K=1,KBM1
       US(I,K) = U(I,K) - (UTMP-UARDS(I))/DT1(I)
       VS(I,K) = V(I,K) - (VTMP-VARDS(I))/DT1(I)
     END DO
   END DO

   END SELECT

   RETURN
   END SUBROUTINE ADJUST2D3D
!==============================================================================|
