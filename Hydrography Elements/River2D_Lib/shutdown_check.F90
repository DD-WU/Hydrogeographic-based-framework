
!==============================================================================|
!  CHECK DEPTH ARRAY FOR NAN.  SHUTDOWN IF FOUND                               |
!==============================================================================|

   SUBROUTINE SHUTDOWN_CHECK 

!==============================================================================|
   USE ALL_VARS
   IMPLICIT NONE
   REAL(DP) :: SBUF,RBUF  
   INTEGER  :: IERR
!==============================================================================|

   !Collect Depth Average to Master Processor
   SBUF = SUM(DBLE(D1(1:N)))
   RBUF = SBUF


   !Halt FVCOM if Depth Average = NaN          
   IF(RBUF /= RBUF)THEN
     IF(MSR)THEN
       WRITE(*,*)'NON FINITE DEPTH FOUND'
       WRITE(*,*)'FVCOM MODEL HAS BECOME UNSTABLE'
       WRITE(*,*)'HALTING'
     END IF
     CALL PSTOP
   END IF


   RETURN
   END SUBROUTINE SHUTDOWN_CHECK
