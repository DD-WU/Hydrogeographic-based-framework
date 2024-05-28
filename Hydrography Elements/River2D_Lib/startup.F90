!==============================================================================|
!   Begin Restart Run From Specified Time                                      |
!==============================================================================|

   SUBROUTINE STARTUP             

!------------------------------------------------------------------------------|

   USE ALL_VARS

! for water quality model  JQI
   USE MOD_WQM
! end for water quality model  JQI

   USE MOD_WD
!QXU{  
   USE MOD_DYE
!QXU}
   USE BCS
   IMPLICIT NONE

!==============================================================================|
   

!
!--Set Water Depth-Using Bathymetry and Free Surface Elevation-----------------!
!

   CALL WATER_DEPTH

     IF(WET_DRY_ON) CALL SET_WD_DATA
!
!--Set up Temperature, Salinity, and Turbulence Quantity Fields----------------!
! 
      
   IF((RESTART == 'cold_start').AND.(S_TYPE == 'non-julian'))THEN

     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    COLD_START'
     IF(MSR)WRITE(IPT,*)  '!  S_TYPE                :    NON-JULIAN'
     CALL INITIAL_TS
     IF(WET_DRY_ON) CALL SET_WD_DATA
     CALL INITIAL_QQL
    if  (WATER_QUALITY==1) CALL INITIAL_WQM

!QXU{
  if (DYE_RELEASE==1)   CALL INITIAL_DYE
  
!QXU}

   ELSE IF((RESTART=='cold_start').AND.(S_TYPE=='julian'))THEN

     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    COLD_START'
     IF(MSR)WRITE(IPT,*)  '!  S_TYPE                :    JULIAN'
     CALL INITIAL_TS
     CALL INITIAL_UVEL
     IF(WET_DRY_ON) CALL SET_WD_DATA
     CALL INITIAL_QQL
    if  (WATER_QUALITY==1)then
     CALL INITIAL_WQM
    endif
!QXU{
  if (DYE_RELEASE==1)then
     CALL INITIAL_DYE
  endif
!QXU}
   !CALL WATER_DEPTH !wanghaocheng: after judge dry and wet, el and d has changed, why not recaculate water depth?
   ELSE IF((RESTART=='hot_cold_s').AND.(S_TYPE=='julian'))THEN
          
     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    HOT_COLD_S'
     IF(MSR)WRITE(IPT,*)  '!  S_TYPE                :    JULIAN'
     CALL HOT_START_DATA
     CALL INITIAL_TS
    if (WATER_QUALITY==1)then
     CALL INITIAL_WQM
    endif
     IF(WET_DRY_ON) CALL WD_READ 
!QXU{
  if  (DYE_RELEASE==1)then
     CALL INITIAL_DYE
  endif
!QXU}

   ELSE IF(RESTART == 'hot_start') THEN
     
     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    HOT_START'
     CALL HOT_START_DATA
     IF(MSR)WRITE(IPT,*)  '!  RESTART DATA          :    READ     '
     IF(WET_DRY_ON) CALL WD_READ 
!QXU{
  if  (DYE_RELEASE==1)then
     CALL INITIAL_DYE
  endif
!QXU}

   ELSE
         
     PRINT*,'RESTAR AND S_TYPE DEFINITION NOT CORRECT'
     PRINT*,'RESTAR==',RESTART
     PRINT*,'S_TYPE==',S_TYPE
     CALL PSTOP
         
   END IF


!
!--Set Values in the Halos-----------------------------------------------------!
! 

   IF(SERIAL)RETURN

   RETURN
   END SUBROUTINE STARTUP
!==============================================================================|


!==============================================================================|
!   Exchange All Flow Variables                                                |
!==============================================================================|

   SUBROUTINE EXCHANGE_ALL 

!------------------------------------------------------------------------------|

   USE ALL_VARS


!QXU}
   IMPLICIT NONE

!==============================================================================|
   



   RETURN
   END SUBROUTINE EXCHANGE_ALL
!==============================================================================|
