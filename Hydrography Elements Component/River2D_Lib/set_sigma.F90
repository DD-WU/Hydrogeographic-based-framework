!==============================================================================|
! This program is used to set up the sigma coordinate in the vertical.         !
!								               !
! sigma levels are determined by a formula of                                  !
!                      sigma(k)=-[(k-1)/(kb-1)]^k11                            !
!    p_sigma=1: uniform sigma layers                                           !
!    p_sigma=2: layers satisfying a parabolic function with high               !
!               vertical resolution near the surface and bottom.               !
!    p_sigma can be used any real number                                       !
!									       !
!  calculates: z(kb) sigma levels					       !
!  calculates: dz(kb-1) delta between sigma levels		               !
!  calculates: zz(kb-1) intra-sigma levels				       !
!  calculates: dzz(kb-2) delta between intra-sigma levels		       !
!==============================================================================|

   SUBROUTINE SET_SIGMA           

!==============================================================================|
   USE ALL_VARS
   IMPLICIT NONE
   INTEGER :: K
!==============================================================================|
!pengfei 2006/2/1 
   IF(P_SIGMA > 1 .AND. MOD(KB,2) == 0)THEN
     IF(MSR) PRINT*, 'kb shoude be an odd number,stop ....'
     CALL PSTOP
   END IF
!pengfei 2006/2/1
!--------  SET SIGMA LEVELS  --------------------------------------------------!  
!orginal formula to set sigma   
   IF(P_SIGMA == 1)THEN
     DO K=1,KB
       Z(K) = -((K-1)/FLOAT(KB-1))**P_SIGMA 
     END DO
   ELSE
!pengfei 2006/2/1
     DO K=1,(KB+1)/2
       Z(K) = -((K-1)/FLOAT((KB+1)/2-1))**P_SIGMA/2 
     END DO
     DO K=(KB+1)/2+1,KB
       Z(K) = ((KB-K)/FLOAT((KB+1)/2-1))**P_SIGMA/2-1.0
     END DO
!pengfei 2006/2/1
   END IF
!---------COMPUTE SIGMA DERIVATIVES AND INTRA SIGMA LEVELS---------------------!
   
   DO K=1,KB-1
     DZ(K)  = Z(K)-Z(K+1)
     ZZ(K)  = .5_SP*(Z(K)+Z(K+1))
   END DO
   ZZ(KB) = 2.0_SP*ZZ(KB-1)-ZZ(KB-2)

   DO K=1,KBM2
     DZZ(K) = ZZ(K)-ZZ(K+1)
   END DO
   DZZ(KBM1) = 0.0_SP
   DZ(KB)    = 0.0_SP

!----------OUTPUT VALUES-TO INFOFILE-------------------------------------------!

   IF(MSR)THEN
     WRITE(IPT,*  )'!'
     WRITE(IPT,*  )'!'
     WRITE(IPT,*)'!                SIGMA LAYER INFO     '
     WRITE(IPT,70)
     DO K=1,KB
       WRITE(IPT,80) K,Z(K),ZZ(K),DZ(K),DZZ(K)
     END DO
     WRITE(IPT,*  )'!'
   END IF

!----------FORMAT STATEMENTS---------------------------------------------------!

70 FORMAT(2x,'k',13x,'z',11x,'zz',11x,'dz',11x,'dzz')
80 FORMAT(' ',i5,4f13.8)

   RETURN
   END SUBROUTINE SET_SIGMA
!==============================================================================|
