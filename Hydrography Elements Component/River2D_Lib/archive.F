!==============================================================================|
!     ARCHIVE THE MODEL RESULTS                                                |
!==============================================================================|

   SUBROUTINE ARCHIVE

!==============================================================================|
#  if defined (NETCDF_IO)
   USE MOD_NCDIO
   USE MOD_NCDAVE
#  endif
#  if defined (WET_DRY)
   USE MOD_WD
#  endif   
!QXU{
#  if defined (BALANCE_2D)
   USE MOD_BALANCE_2D
#  endif
!qxu}
#  if defined (SEDIMENT)
   USE MOD_SED
#  endif
   USE CONTROL  
   IMPLICIT NONE
!==============================================================================|

!
!--DUMP MEDM FILE EVERY IRECORD  ITERATIONS (NO DUMP IF IRECORD  = 0)----------!
!
   IF(IRECORD /= 0)THEN
     IF(MOD(IINT,IRECORD) == 0)THEN
       CALL OUT_BINARY(IINT/IRECORD)
     END IF
   END IF
!
!----------------------DUMP SMS FILE-------------------------------------------!
!
   IF(IDMPSMS /= 0)THEN
     IF(MOD(IINT,IDMPSMS) == 0)THEN
       IF(MSR)WRITE(IPT,*)  '!  DUMPING               :    SMS FILE'    
       CALL OUT_SMS_ONE
     END IF
   END IF
!
!--DUMP RESTART FILE EVERY IRESTART ITERATIONS (NO DUMP IF IRESTART = 0)-------!
!
   IF(IRESTART /= 0)THEN
     IF(MOD(IINT,IRESTART) == 0)THEN
       IF(MSR)WRITE(IPT,*)  '!  DUMPING               :    RESTART FILE'
       CALL ARCRST

# if defined(WET_DRY)
       IF(WET_DRY_ON) CALL WD_DUMP
# endif       
# if defined(SEDIMENT)
       IF(SEDIMENT_ON) CALL ARCHIVE_SED
# endif       
     END IF
   END IF

!
!--FLOW FIELD AVERAGING AND OUTPUT OF FLOWFIELD AVERAGES-----------------------!
!
   IF(AVGE_ON) CALL OUT_AVGE
#  if defined (NETCDF_IO)
   IF(CDF_OUT_AVE) CALL AVGE_FIELDS(IINT) 
#  endif

!
!--NETCDF OUTPUT---------------------------------------------------------------!
!
#  if defined (NETCDF_IO)
   IF(CDF_INT /= 0 .AND. CDF_OUT)THEN
     IF(MOD(IINT,CDF_INT)==0) CALL OUT_NETCDF
   END IF
#  endif

!QXU{
#  if defined (BALANCE_2D)
   IF(OUT_BALANCE)THEN
     CALL OUT_TIMESERIES_BALANCE
   ENDIF
#  endif
!QXU}
   RETURN
   END SUBROUTINE ARCHIVE
!==============================================================================|
