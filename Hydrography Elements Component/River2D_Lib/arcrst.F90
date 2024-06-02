!==============================================================================|
!   DUMP DATA FILE FOR RESTART                                                 |
!==============================================================================|

   SUBROUTINE ARCRST            

!------------------------------------------------------------------------------|

   USE ALL_VARS

   USE MOD_WQM

   USE MOD_DYE

   USE MOD_EQUITIDE
   IMPLICIT NONE
   REAL(SP), ALLOCATABLE :: RTP(:)
   INTEGER I,K,ME,NPC,N1
!==============================================================================|
   
   ME = MYID ; NPC = NPROCS 
   ALLOCATE(RTP(N)) ; RTP = 0.0_SP

   IF(MSR)THEN
     OPEN(1,FILE='restart',FORM='UNFORMATTED',STATUS='UNKNOWN')
     REWIND(1)
     WRITE(1) IINT
   END IF

   IF(SERIAL)THEN
     WRITE(1) ((U(I,K),    K=1,KB),I=0,N)
     WRITE(1) ((V(I,K),    K=1,KB),I=0,N)
     WRITE(1) ((W(I,K),    K=1,KB),I=0,N)
     WRITE(1) ((Q2(I,K),   K=1,KB),I=0,M)
     WRITE(1) ((Q2L(I,K),  K=1,KB),I=0,M)
     WRITE(1) ((L(I,K)  ,  K=1,KB),I=0,M)
     WRITE(1) ((S(I,K),    K=1,KB),I=0,N)
     WRITE(1) ((T(I,K),    K=1,KB),I=0,N)
     WRITE(1) ((RHO(I,K),  K=1,KB),I=0,N)
     WRITE(1) ((TMEAN(I,K),K=1,KB),I=0,N)
     WRITE(1) ((SMEAN(I,K),K=1,KB),I=0,N)
     WRITE(1) ((RMEAN(I,K),K=1,KB),I=0,N)

     WRITE(1) ((S1(I,K),    K=1,KB),I=1,M)
     WRITE(1) ((T1(I,K),    K=1,KB),I=1,M)
     WRITE(1) ((RHO1(I,K),  K=1,KB),I=1,M)
     WRITE(1) ((TMEAN1(I,K),K=1,KB),I=1,M)
     WRITE(1) ((SMEAN1(I,K),K=1,KB),I=1,M)
     WRITE(1) ((RMEAN1(I,K),K=1,KB),I=1,M)

     WRITE(1) ((KM(I,K),K=1,KB),I=1,M)
     WRITE(1) ((KH(I,K),K=1,KB),I=1,M)
     WRITE(1) ((KQ(I,K),K=1,KB),I=1,M)

     WRITE(1) (UA(I), I=0,N)
     WRITE(1) (VA(I), I=0,N)

     WRITE(1) (EL1(I), I=1,N)
     WRITE(1) (ET1(I), I=1,N)
     WRITE(1) (H1(I),  I=1,N)
     WRITE(1) (D1(I),  I=1,N)
     WRITE(1) (DT1(I), I=1,N)
     WRITE(1) (RTP(I), I=1,N)

     WRITE(1) (EL(I), I=1,M)
     WRITE(1) (ET(I), I=1,M)
     WRITE(1) (H(I),  I=1,M)
     WRITE(1) (D(I),  I=1,M)
     WRITE(1) (DT(I), I=1,M)

    if  (EQUI_TIDE==1) WRITE(1) (EL_EQI(I), I=1,M)



    if  (WATER_QUALITY==1)then
     DO N1 = 1, NB
       WRITE(1) ((WQM(I,K,N1),K=1,KB),I=1,M)
     END DO
    endif
!QXU
    if (DYE_RELEASE==1)then
     IF(IINT.GT.IINT_SPE_DYE_B) THEN
     WRITE(1) ((DYE(I,K),K=1,KB),I=1,M)
     WRITE(1) ((DYEMEAN(I,K),K=1,KB),I=1,M)
     ENDIF
    endif
!QXU


   END IF
   IF(MSR) CLOSE(1)
   DEALLOCATE(RTP)

   RETURN
   END SUBROUTINE ARCRST
!==============================================================================|



