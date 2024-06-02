!==============================================================================|
!     this subroutine is used to calculate the true temperature                !
!     and salinity by including vertical diffusion implicitly.                 !
!==============================================================================|

   SUBROUTINE VDIF_TS_GOM(NCON1,F)                

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE BCS
   IMPLICIT NONE
   INTEGER :: I,K,NCON1,J,KI
   REAL(DP) :: TMP,TMP1,TMP2,TMP3,QTMP,GW,ZDEP,FKH,UMOLPR
   REAL(SP), DIMENSION(0:MT,KB)  :: F
   REAL(DP), DIMENSION(M,KB)     :: FF,AF,CF,VHF,VHPF,RAD
   REAL(DP), DIMENSION(M)        :: KHBOTTOM,WFSURF,SWRADF

   UMOLPR = UMOL*1.E0_DP

!------------------------------------------------------------------------------!
!                                                                              !
!        the following section solves the equation                             !
!         dti*(kh*f')'-f=-fb                                                   !
!                                                                              !
!------------------------------------------------------------------------------!


   DO K = 2, KBM1
     DO I = 1, M
!       FKH=0.0_SP
!       DO J=1,NTVE(I)
!         FKH=FKH+KH(NBVE(I,J),K)
!       END DO
!       FKH=FKH/FLOAT(NTVE(I))       
       FKH = KH(I,K)
    
       IF(K == KBM1) THEN
         KHBOTTOM(I)=FKH
       END IF

       AF(I,K-1)=-DTI*(FKH+UMOLPR)/(DZ(K-1)*DZZ(K-1)*D(I)*D(I))
       CF(I,K)=-DTI*(FKH+UMOLPR)/(DZ(K)*DZZ(K-1)*D(I)*D(I))
     END DO
   END DO


!------------------------------------------------------------------------------!
!     the net heat flux input.                                                 !
!     the method shown below can be used when we turn off the                  !
!     body force in subroutine advt. be sure this method could                 !
!     cause the surface overheated if the vertical resolution                  !
!     is not high enough.                                                      !
!------------------------------------------------------------------------------!

   IF(NCON1 == 1) THEN
       
   IF(H_TYPE == 'flux_h')THEN
     DO I=1,M
       WFSURF(I)=WTSURF(I)
       SWRADF(I)=SWRAD(I)
     END DO
     DO K = 1, KB
       DO I = 1, M
           ZDEP = Z(K)*D(I)

!LIU       SPECIFIED PARAMETER FOR GOM CASE (JUNE 18, 2003)
           IF(D(I) < 100.0_SP) THEN
             RHEAT=0.78_SP
             ZETA1=1.4_SP
             ZETA2=6.3_SP
            ELSE IF(D(I) > 150.0_SP) THEN
             RHEAT=0.58_SP
             ZETA1=0.35_SP
             ZETA2=23.0_SP
            ELSE
             TMP1=(D(I)-100.0_SP)/(150.0_SP-100.0_SP)
             RHEAT=0.78_SP+(0.58_SP-0.78_SP)*TMP1
             ZETA1=1.4_SP+ (0.35_SP-1.4_SP )*TMP1
             ZETA2=6.3_SP+ (23.0_SP-6.3_SP )*TMP1
           END IF
!LIU       END (JUNE 18, 2003)         

           RAD(I,K)=SWRADF(I)*(RHEAT*EXP(ZDEP/ZETA1)+(1-RHEAT)*EXP(ZDEP/ZETA2))
       END DO
     END DO

   ELSE !! H_TYPE = 'body_h'   
     RAD    = 0.0_SP
     WFSURF = 0.0_SP
     SWRADF = 0.0_SP
   END IF 
        
   ELSE IF(NCON1 == 2) THEN

     DO I = 1, M
       SWRADF(I)= 0.0_SP
       WFSURF(I)=0.0_SP
       DO K=1,KB
         RAD(I,K)=0.0_SP
       END DO
     END DO

   ELSE
     PRINT*,'NCON NOT CORRECT IN VDIF_TS',NCON1
     CALL PSTOP
   END IF

!------------------------------------------------------------------------------!
!   surface bcs; wfsurf                                                        !
!------------------------------------------------------------------------------!

   DO I = 1, M
     VHF(I,1) = AF(I,1) / (AF(I,1)-1.)
     VHPF(I,1) = -DTI *(WFSURF(I)-SWRADF(I) &
                 +RAD(I,1)-RAD(I,2)) / (-DZ(1)*D(I)) - F(I,1)
     VHPF(I,1) = VHPF(I,1) / (AF(I,1)-1.)
   END DO

   DO K = 2, KBM2
     DO I = 1, M
       VHPF(I,K)=1./ (AF(I,K)+CF(I,K)*(1.-VHF(I,K-1))-1.)
       VHF(I,K) = AF(I,K) * VHPF(I,K)
       VHPF(I,K) = (CF(I,K)*VHPF(I,K-1)-DBLE(F(I,K)) &
                   +DTI*(RAD(I,K)-RAD(I,K+1))/(D(I)*DZ(K)))*VHPF(I,K)
     END DO
   END DO


   DO  K = 1, KBM1
     DO  I = 1, M
       FF(I,K) = F(I,K)
     END DO
   END DO

   DO I = 1, M
     IF (ISONB(I) /= 2) THEN
       TMP1=PFPXB(I)*COS(SITA_GD(I))+PFPYB(I)*SIN(SITA_GD(I))
       TMP2=AH_BOTTOM(I)*PHPN(I)
       TMP3=KHBOTTOM(I)+UMOLPR+AH_BOTTOM(I)*PHPN(I)*PHPN(I)
       TMP=TMP1*TMP2/TMP3*(KHBOTTOM(I)+UMOLPR)
! --- Huang change
       IF (TMP1 > 0.0_SP) TMP=0.0_SP
!       !also try TMP = 0.0_SP
! change end

       GW=0.0_SP
       IF(IBFW > 0) THEN
         DO J=1,IBFW
           IF(I == NODE_BFW(J).AND.(NCON1 == 2)) THEN
             QTMP=-(F(I,KBM1)*D(I)*DZ(KBM1)*BFWDIS(J))/ &
                   (D(I)*DZ(KBM1)*ART1(I)+BFWDIS(J))
             GW=DTI/D(I)/DZ(KBM1)*QTMP
             TMP=0.0_SP
           END IF
         END DO
       END IF
       
       FF(I,KBM1) = ((CF(I,KBM1)*VHPF(I,KBM2)-FF(I,KBM1)-GW &
               +DTI*(RAD(I,KBM1)-RAD(I,KB)-TMP)/(D(I)*DZ(KBM1))) &
                /(CF(I,KBM1)*(1.-VHF(I,KBM2))-1.))
     END IF
   END DO

   DO  K = 2, KBM1
     KI = KB - K
     DO  I = 1, M
       IF (ISONB(I) /= 2) THEN
         FF(I,KI) = (VHF(I,KI)*FF(I,KI+1)+VHPF(I,KI))
       END IF
     END DO
   END DO

   DO  K = 1, KBM1
     DO I=1,M
       F(I,K) = FF(I,K)
     END DO
   END DO

   RETURN
   END SUBROUTINE VDIF_TS_GOM
!==============================================================================|
