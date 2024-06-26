!==============================================================================|
!  MODULE CONTAINING SUBROUTINES USED TO SET UP MOMENTUM BALANCE OUTPUT        |
!==============================================================================|

MODULE MOD_BALANCE_2D
   USE MOD_PREC
   USE MOD_INP
   USE CONTROL
   IMPLICIT NONE
   SAVE
   LOGICAL  :: OUT_BALANCE                 !!TRUE IF MOMENTUM BALANCE CHECHACTIVE
   INTEGER :: NUM_BALANCE,IOMOB
   
   INTEGER,  ALLOCATABLE :: NO_CELL(:)    !!CELL NO FOR OUTPUT MOMENTUM BALANCE
   REAL(SP), ALLOCATABLE :: ADFXA(:)      
   REAL(SP), ALLOCATABLE :: ADFYA(:)

   REAL(SP), ALLOCATABLE :: ADVUA2(:)   !!ADVECTION TERM   
   REAL(SP), ALLOCATABLE :: ADVVA2(:) 
   
   REAL(SP), ALLOCATABLE :: ADFX2(:)    
   REAL(SP), ALLOCATABLE :: ADFY2(:)
   
   REAL(SP), ALLOCATABLE :: DRX2D2(:)   !!BAROCLINIC PRESURE GRADENT FORCE
   REAL(SP), ALLOCATABLE :: DRY2D2(:)
   
   REAL(SP), ALLOCATABLE :: CORX2(:)    !!CORIOLIS FORCE TERM
   REAL(SP), ALLOCATABLE :: CORY2(:)
      
   REAL(SP), ALLOCATABLE :: PSTX2(:)    !!BAROTROPIC PRESURE GRSDENT FORCE
   REAL(SP), ALLOCATABLE :: PSTY2(:)
   
   REAL(SP), ALLOCATABLE :: ADX2D2(:)   !!DIFFUSION TERM (GX,GY)
   REAL(SP), ALLOCATABLE :: ADY2D2(:)
   
   
   REAL(SP), ALLOCATABLE :: WUSURBF2(:)  !!STRESS TERM
   REAL(SP), ALLOCATABLE :: WVSURBF2(:)
   
   REAL(SP), ALLOCATABLE :: DUDT2(:)
   REAL(SP), ALLOCATABLE :: DVDT2(:)
   
   REAL(SP), ALLOCATABLE :: DIVX2D2(:)
   REAL(SP), ALLOCATABLE :: DIVY2D2(:)
   REAL(SP), ALLOCATABLE :: DEDT2(:)
   
   
   REAL(SP), ALLOCATABLE :: ADVUA2_AVE(:)
   REAL(SP), ALLOCATABLE :: ADVVA2_AVE(:)
   
   REAL(SP), ALLOCATABLE :: ADFX2_AVE(:)
   REAL(SP), ALLOCATABLE :: ADFY2_AVE(:)
   
   REAL(SP), ALLOCATABLE :: DRX2D2_AVE(:)
   REAL(SP), ALLOCATABLE :: DRY2D2_AVE(:)
   
   REAL(SP), ALLOCATABLE :: CORX2_AVE(:)
   REAL(SP), ALLOCATABLE :: CORY2_AVE(:)
   
   REAL(SP), ALLOCATABLE :: PSTX2_AVE(:)
   REAL(SP), ALLOCATABLE :: PSTY2_AVE(:)
      
   REAL(SP), ALLOCATABLE :: ADX2D2_AVE(:)
   REAL(SP), ALLOCATABLE :: ADY2D2_AVE(:)
   
   REAL(SP), ALLOCATABLE :: WUSURBF2_AVE(:)
   REAL(SP), ALLOCATABLE :: WVSURBF2_AVE(:)
   
   REAL(SP), ALLOCATABLE :: DUDT2_AVE(:)
   REAL(SP), ALLOCATABLE :: DVDT2_AVE(:)
   
!===================================================================================|
   CONTAINS   !!INCLUDED SUBROUTINES FOLLOW
!===================================================================================|
    SUBROUTINE ALLOC_BALANCE_VARS
    USE LIMS

    ALLOCATE(ADFXA(0:NT))        ;ADFXA     = ZERO
    ALLOCATE(ADFYA(0:NT))        ;ADFYA     = ZERO
    
    ALLOCATE(ADVUA2(0:NT))       ;ADVUA2    = ZERO
    ALLOCATE(ADVVA2(0:NT))       ;ADVVA2    = ZERO
    ALLOCATE(ADFX2(0:NT))        ;ADFX2     = ZERO
    ALLOCATE(ADFY2(0:NT))        ;ADFY2     = ZERO
    ALLOCATE(DRX2D2(0:NT))       ;DRX2D2    = ZERO
    ALLOCATE(DRY2D2(0:NT))       ;DRY2D2    = ZERO
    ALLOCATE(CORX2(0:NT))        ;CORX2     = ZERO
    ALLOCATE(CORY2(0:NT))        ;CORY2     = ZERO
    ALLOCATE(PSTX2(0:NT))        ;PSTX2     = ZERO
    ALLOCATE(PSTY2(0:NT))        ;PSTY2     = ZERO
    ALLOCATE(ADX2D2(0:NT))       ;ADX2D2    = ZERO
    ALLOCATE(ADY2D2(0:NT))       ;ADY2D2    = ZERO
    ALLOCATE(WUSURBF2(0:NT))     ;WUSURBF2  = ZERO 
    ALLOCATE(WVSURBF2(0:NT))     ;WVSURBF2  = ZERO 
    ALLOCATE(DUDT2(0:NT))        ;DUDT2     = ZERO 
    ALLOCATE(DVDT2(0:NT))        ;DVDT2     = ZERO
     
    ALLOCATE(DIVX2D2(0:NT))      ;DIVX2D2   = ZERO 
    ALLOCATE(DIVY2D2(0:NT))      ;DIVY2D2   = ZERO 
    ALLOCATE(DEDT2(0:NT))        ;DEDT2     = ZERO  
   RETURN
   END SUBROUTINE ALLOC_BALANCE_VARS

   SUBROUTINE SET_BALANCE_PARAM
   USE MOD_PREC
   USE CONTROL
   IMPLICIT NONE
   INTEGER  INTVEC(150),ISCAN,KTEMP
   CHARACTER(LEN=120) :: FNAME
   FNAME = trim(casename1)
!------------------------------------------------------------------------------|
!     "OUT_BALANCE"   !! 
!------------------------------------------------------------------------------|     
   ISCAN = SCAN_FILE(TRIM(FNAME),"OUT_BALANCE",LVAL = OUT_BALANCE)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING OUT_BALANCE: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF
!------------------------------------------------------------------------------|
!     "NUM_BALANCE"   !! 
!------------------------------------------------------------------------------|  
   ISCAN = SCAN_FILE(FNAME,"NUM_BALANCE",ISCAL = NUM_BALANCE)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING NUM_BALANCE: ',ISCAN
     IF(ISCAN == -2)THEN
       WRITE(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',TRIM(FNAME)
     END IF
     CALL PSTOP 
   END IF
   
!------------------------------------------------------------------------------|
!     "NO_CELL"   !! 
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"NO_CELL",IVEC =INTVEC ,NSZE = KTEMP)
   IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING NO_CELL: ',ISCAN
     CALL PSTOP
   END IF
   IF(MSR)THEN
     IF(KTEMP /= NUM_BALANCE)THEN
       WRITE(*,*)'NUMBER OF SPECIFIED NO_CELL IS NOT EQUAL TO NUM_BALANCE' 
       WRITE(*,*)'NUM_BALANCE: ',NUM_BALANCE
       WRITE(*,*)'NO_CELL: ',INTVEC
     END IF
   END IF
  
   ALLOCATE(NO_CELL(NUM_BALANCE)) ; NO_CELL=0
   NO_CELL(1:NUM_BALANCE)= INTVEC(1:NUM_BALANCE)

!==============================================================================|
!            SCREEN REPORT OF SET MOMENTUM BALANCE OUT VARIABlES               !
!==============================================================================|
   IF(MSR) THEN  
     WRITE(IPT,*) '!                                                   !'     
     WRITE(IPT,*) '!------SPECIFY MOMENTUM BALANCE OUT VARIABlES-------!'     
     WRITE(IPT,*) '!                                                   !'     
     WRITE(IPT,*) '!  # OUT_BALANCE         :',OUT_BALANCE
     WRITE(IPT,*) '!  # NUM_BALANCE         :',NUM_BALANCE
     WRITE(IPT,*) '!  # NO_CELL             :',NO_CELL
   END IF
   
   
   
   RETURN
   END SUBROUTINE SET_BALANCE_PARAM

!
!  out time series of momentum balance terms
!   
   SUBROUTINE OUT_TIMESERIES_BALANCE
   USE MOD_PREC
   USE ALL_VARS

   IMPLICIT NONE
   INTEGER I
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ADVUA2TMP,  ADVVA2TMP,  ADFX2TMP, ADFY2TMP  
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: DRX2D2TMP,  DRY2D2TMP,  CORX2TMP, CORY2TMP  
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: PSTX2TMP,   PSTY2TMP,   ADX2D2TMP,ADY2D2TMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: WUSURBF2TMP,WVSURBF2TMP,DUDT2TMP, DVDT2TMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: DIVX2D2TMP,    DIVY2D2TMP,    DEDT2TMP
   IF(SERIAL)THEN
     WRITE(IOMOB,'(i6,150(19E13.5,2X))') IINT, &
       (ADVUA2(NO_CELL(I)),  ADVVA2(NO_CELL(I)),&
        ADFX2(NO_CELL(I)),   ADFY2(NO_CELL(I)),& 
        DRX2D2(NO_CELL(I)),  DRY2D2(NO_CELL(I)),&
	CORX2(NO_CELL(I)),   CORY2(NO_CELL(I)),&  
        PSTX2(NO_CELL(I)),   PSTY2(NO_CELL(I)),&
	ADX2D2(NO_CELL(I)),  ADY2D2(NO_CELL(I)),&
        WUSURBF2(NO_CELL(I)),WVSURBF2(NO_CELL(I)),&
	DUDT2(NO_CELL(I)),   DVDT2(NO_CELL(I)),&
	DIVX2D2(NO_CELL(I)),DIVY2D2(NO_CELL(I)),&
	DEDT2(NO_CELL(I)),I=1,NUM_BALANCE)
   ENDIF 

   RETURN
 
  END SUBROUTINE OUT_TIMESERIES_BALANCE
   
   
END MODULE  MOD_BALANCE_2D
