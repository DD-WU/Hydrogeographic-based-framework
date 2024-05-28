!==============================================================================|
!   Implement Fresh Water Boundary Condition on Bottom (Groundwater)           |
!==============================================================================|

   SUBROUTINE BCOND_BFW           

!------------------------------------------------------------------------------|
   USE ALL_VARS
   USE BCS
   IMPLICIT NONE
   REAL(SP) :: FACT,UFACT
   INTEGER  :: L1,L2,IERR 
!------------------------------------------------------------------------------|

   IF(IBFW <= 0)RETURN

   CALL BRACKET(BFW_TM,THOUR,L1,L2,FACT,UFACT,IERR)
   BFWDIS(:)=UFACT*BFWQDIS(:,L1)+FACT*BFWQDIS(:,L2)
   BFWDIS   =BFWDIS*RAMP

   RETURN
   END SUBROUTINE BCOND_BFW
!==============================================================================|
