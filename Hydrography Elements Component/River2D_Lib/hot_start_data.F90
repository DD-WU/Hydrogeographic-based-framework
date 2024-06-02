
!==============================================================================|
!   READ IN RESTART DATA FILE AND RESTART                                      |
!==============================================================================|

   SUBROUTINE HOT_START_DATA      

!------------------------------------------------------------------------------|

   USE ALL_VARS

   USE MOD_WQM
   USE MOD_DYE


   USE MOD_EQUITIDE
   IMPLICIT NONE
   INTEGER :: I,K,N1
   REAL(SP), DIMENSION(N) ::RTP 
!==============================================================================|
!  NOTE: TO MAINTAIN COMPATIBILITY WITH PREVIOUS FVCOM, ARRAY DTF1 (NO LONGER  !
!  USED) IS READ IN AS RTP                                                     !
!==============================================================================|



   IF(SERIAL)THEN
     REWIND(INRES)
     READ(INRES) IINT
     READ(INRES) ((U(I,K),K=1,KB),I=0,N)
     READ(INRES) ((V(I,K),K=1,KB),I=0,N)
     READ(INRES) ((W(I,K),K=1,KB),I=0,N)
   READ(INRES) ((Q2(I,K),K=1,KB),I=0,M)
     READ(INRES) ((Q2L(I,K),K=1,KB),I=0,M)
     READ(INRES) ((L(I,K),K=1,KB),I=0,M)
     READ(INRES) ((S(I,K),K=1,KB),I=0,N)
     READ(INRES) ((T(I,K),K=1,KB),I=0,N)
     READ(INRES) ((RHO(I,K),K=1,KB),I=0,N)
     READ(INRES) ((TMEAN(I,K),K=1,KB),I=0,N)
     READ(INRES) ((SMEAN(I,K),K=1,KB),I=0,N)
     READ(INRES) ((RMEAN(I,K),K=1,KB),I=0,N)

     READ(INRES) ((S1(I,K),K=1,KB),I=1,M)
     READ(INRES) ((T1(I,K),K=1,KB),I=1,M)
     READ(INRES) ((RHO1(I,K),K=1,KB),I=1,M)
     READ(INRES) ((TMEAN1(I,K),K=1,KB),I=1,M)
     READ(INRES) ((SMEAN1(I,K),K=1,KB),I=1,M)
     READ(INRES) ((RMEAN1(I,K),K=1,KB),I=1,M)
     READ(INRES) ((KM(I,K),K=1,KB),I=1,M)
     READ(INRES) ((KH(I,K),K=1,KB),I=1,M)
     READ(INRES) ((KQ(I,K),K=1,KB),I=1,M)

     READ(INRES) (UA(I), I=0,N)
     READ(INRES) (VA(I), I=0,N)
     READ(INRES) (EL1(I), I=1,N)
     READ(INRES) (ET1(I), I=1,N)
     READ(INRES) (H1(I), I=1,N)
     READ(INRES) (D1(I), I=1,N)
     READ(INRES) (DT1(I), I=1,N)
!    READ(INRES) (DTF1(I), I=1,N)
     READ(INRES) (RTP(I), I=1,N)

     READ(INRES) (EL(I), I=1,M)
     READ(INRES) (ET(I), I=1,M)
     READ(INRES) (H(I), I=1,M)
     READ(INRES) (D(I), I=1,M)
     READ(INRES) (DT(I), I=1,M)
      if(EQUI_TIDE==1)then
     READ(INRES) (EL_EQI(I), I=1,M)
     endif
    if (WATER_QUALITY==1)then
     DO N1=1,NB
       READ(INRES) ((WQM(I,K,N1),K=1,KB),I=1,M)
     END DO
    endif
!QXU{
    if (DYE_RELEASE==1)then
       IF(IINT.GT.IINT_SPE_DYE_B) THEN
         READ(INRES) ((DYE(I,K),K=1,KB),I=1,M)
         READ(INRES) ((DYEMEAN(I,K),K=1,KB),I=1,M)
       ENDIF
    endif
!QXU}     

     CLOSE(INRES)

   END IF

!--Set Turbulent Macro-Scale


   CALL N2E3D(KM,KM1)

   RETURN
   END SUBROUTINE HOT_START_DATA
!==============================================================================|
