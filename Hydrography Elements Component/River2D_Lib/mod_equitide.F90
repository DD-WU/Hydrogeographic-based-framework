!==============================================================================!
!
!==============================================================================!

MODULE MOD_EQUITIDE
  
   
 !  USE MOD_ROTATION   
   USE CONTROL
   IMPLICIT NONE
   SAVE


   REAL(SP), ALLOCATABLE :: ELF_EQI(:),ELRK_EQI(:),EL_EQI(:),EGF_EQI(:)
   REAL(SP), ALLOCATABLE :: PHI(:),LAMDA(:)
   REAL(SP), ALLOCATABLE :: BATA(:)
   
!--The amplitudes of equilibrium tides----------------------------------------------
   REAL(SP),PARAMETER, DIMENSION(6) :: &

!                          S2           M2           N2           K1          P1     
           APT_EQI = (/0.112841_SP, 0.242334_SP, 0.046398_SP, 0.141565_SP, 0.00_SP, &
!                          O1	   
	      	       0.100514_SP/)  
!--The frequncy of equilibrium tides------------------------------------------------
   REAL(SP),DIMENSION(6) :: FREQ_EQI
!   REAL(SP),PARAMETER, DIMENSION(6) :: &
!                           S2              M2              N2              K1            
!           FREQ_EQI = (/1.454441E-4_SP, 1.405189E-4_SP, 1.378791E-4_SP, 0.729221E-4_SP, &
	        
!                           P1              O1
!	                0.000000E-4_SP, 0.675981E-4_SP/)  
!--The Love number k-----------------------------------------------------------------
   REAL(SP),PARAMETER, DIMENSION(6) :: &
!                        S2        M2        N2        K1        P1        O1
           K_LOVE = (/0.302_SP, 0.302_SP, 0.302_SP, 0.256_SP, 0.000_SP, 0.298_SP/)   
!--The Love number h -----------------------------------------------------------------
   REAL(SP),PARAMETER, DIMENSION(6) :: &
!                        S2        M2        N2        K1        P1        O1
           H_LOVE = (/0.602_SP, 0.602_SP, 0.602_SP, 0.520_SP, 0.000_SP, 0.603_SP/)   	   

   INTEGER   :: APT_FACT_EQUI(6)
   
   CONTAINS


!==========================================================================|
!==========================================================================|
   SUBROUTINE ALLOCATE_EQUI
   USE ALL_VARS
   IMPLICIT NONE
   
   INTEGER :: I
   
   ALLOCATE(ELF_EQI(0:MT)); ELF_EQI = ZERO
   ALLOCATE(ELRK_EQI(0:MT)); ELRK_EQI = ZERO
   ALLOCATE(EL_EQI(0:MT)); EL_EQI = ZERO
   ALLOCATE(EGF_EQI(0:MT)); EGF_EQI = ZERO
   ALLOCATE(PHI(0:MT))   ; PHI     = ZERO
   ALLOCATE(LAMDA(0:MT))   ; LAMDA   = ZERO
   ALLOCATE(BATA(6))      ; BATA    = ZERO
   
  if  (SPHERICAL==1)THEN
   IF(SERIAL)THEN
     PHI   = YG
     LAMDA = XG
   END IF  

  else
   IF(MSR) PRINT*,"THE EQUILIBRIUM TIDE HAS NOT BEEN ADDED IN THE ",     &
                  "NON SPHERICAL COORDINATE"
   CALL PSTOP		    
  endif     

   DO I=1,6
     BATA(I) = 1.0_SP+K_LOVE(I)-H_LOVE(I)
   END DO  
          
   DO I=1,6
     FREQ_EQI(I) = PI2/PERIOD(I)
     print*,freq_eqi(i),period(i),'tide'
   END DO  

   RETURN
   END SUBROUTINE ALLOCATE_EQUI
!==========================================================================|


!==========================================================================|
   SUBROUTINE ELEVATION_EQUI

!--------------------------------------------------------------------------|
!  Surface Elevation of EQUILIBRIUM TIDE                                   |
!--------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_OBCS
   IMPLICIT NONE

   INTEGER :: I,J
   REAL(SP):: TIME1
   REAL(SP):: FORCE,PHAI_IJ

   TIME1 = TIME * 86400.0_SP
!
!-Julian: Set Elevation Based on Linear Interpolation Between Two Data Times-|
!
   IF(S_TYPE == 'julian')THEN
!  not finish yet
   END IF

!
!-Non-Julian: Set Elevation of Equilibrium Tide -----------------------------|
!

   IF(S_TYPE == 'non-julian')THEN
     DO I = 1, MT
       FORCE = 0.0_SP
       DO J = 1,6
         IF(J <= 3)THEN
           PHAI_IJ = LAMDA(I)*PI2/360.0_SP
           FORCE = BATA(J)*APT_EQI(J)*APT_FACT_EQUI(J)*COS(PHI(I)*PI2/360.0_SP)**2                  &
	           *COS(FREQ_EQI(J)*TIME1+2.0_SP*PHAI_IJ) + FORCE
	 ELSE
           PHAI_IJ = LAMDA(I)*PI2/360.0_SP
           FORCE = BATA(J)*APT_EQI(J)*APT_FACT_EQUI(J)*SIN(2.0_SP*PHI(I)*PI2/360.0_SP)              &
	           *COS(FREQ_EQI(J)*TIME1+PHAI_IJ) + FORCE
	 END IF  
       END DO
       FORCE = FORCE
       ELF_EQI(I) = FORCE * RAMP
     END DO

   END IF

   RETURN
   END SUBROUTINE ELEVATION_EQUI
!============================================================================|
!============================================================================|

END MODULE MOD_EQUITIDE
