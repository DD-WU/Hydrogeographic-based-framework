!==============================================================================|

   SUBROUTINE RHO_MEAN
!==============================================================================|
   USE ALL_VARS
   IMPLICIT NONE

   INTEGER, PARAMETER :: KBB=601
   INTEGER, PARAMETER :: KBBM1=KBB-1

   REAL(SP) DELTZ
   REAL(SP), DIMENSION(KBBM1)   :: PHY_Z   !Depth(m) in every standary Z levels 
   REAL(SP), DIMENSION(KBBM1)   :: RHOZTMP !density in standary Z levels 
   REAL(SP), DIMENSION(M,KBBM1) :: RHOZ    !density in standary Z levals,in Nodes 
   REAL(SP), DIMENSION(KBBM1)   :: RHOA    !density mean in standary Z levals
   
   REAL(SP), DIMENSION(KBM1)    :: ZM      !Depth (m) in every sigma levels for giving node
   REAL(SP), DIMENSION(KBM1)    :: RHOS    !density for giving node
   INTEGER ICOUNT,I,K

!--CALCULATE Z-LEVELS TO MAX DEPTH---------------------------------------------|
   
   DELTZ=HMAX/FLOAT(KBBM1)

   DO K=1,KBBM1
     PHY_Z(K)=(0.5_SP-FLOAT(K))*DELTZ
   END DO

!--LINEARLY INTERPOLATE TO OBTAIN DENSITY VALUES AT Z LEVELS-------------------|
   
   DO I=1,M
     DO K=1,KBM1
       ZM(K)=ZZ(K)*DT(I)+EL(I)
       RHOS(K)= RHO1(I,K)                   !=DBLE(RHO1(I,K))
     END DO
     
     CALL SINTER_P(ZM,RHOS,PHY_Z,RHOZTMP,KBM1,KBBM1)

     DO K=1,KBBM1
       RHOZ(I,K)=RHOZTMP(K)
     END DO
   END DO

!--DO THE AVERAGE OVER Z_levels 

   DO K=1,KBBM1
      RHOA(K)=0.0
      ICOUNT=0
      DO I=1,M
         IF(-H(I).LE.PHY_Z(K)) THEN
           ICOUNT = ICOUNT + 1
           RHOA(K)=RHOA(K)+RHOZ(I,K)
         END IF
      ENDDO

      IF(ICOUNT.GE.1) then
         RHOA(K) = RHOA(K)/float(ICOUNT)
      ELSE 
         RHOA(K) = RHOA(K-1) 
      END IF
!        write(*,'("K,RHOA= ",f10.6)') K,RHOA(K)
   END DO

!--LINEARLY INTERPOLATE TO OBTAIN DENSITY VALUES AT SIGMA LEVELS-------------------|
        

   DO I=1,M
      DO K=1,KBM1
         ZM(K)=ZZ(K)*DT(I)+EL(I)
      END DO

      CALL SINTER(PHY_Z,RHOA,ZM,RHOS,KBBM1,KBM1)
      DO K=1,KBM1
         RMEAN1(I,K)=RHOS(K)                  !=REAL(RHOS(K))
      END DO
   END DO
   DO I=1,N
      DO K=1,KBM1
         RMEAN(I,K) = (RMEAN1(NV(I,1),K) + RMEAN1(NV(I,2),K) + RMEAN1(NV(I,3),K))/3.0_SP
      END DO
   END DO
   
   RETURN
   END
   
