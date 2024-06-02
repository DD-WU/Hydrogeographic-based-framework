!==============================================================================|
!   READ IN STATIC WATER DEPTH AND CALCULATE RELATED QUANTITIES                |
!                                                                              |
!   INPUTS: H(NNODE) BATHYMETRIC DEPTH AT NODES				       |
!   INITIALIZES: D(NNODE) DEPTH AT NODES				       |
!   INITIALIZES: DT(NNODE) ???					               |
!   INITIALIZES: H1(NNODE) BATHYMETRIC DEPTH AT ELEMENTS		       |
!   INITIALIZES: D1(NNODE) DEPTH AT NODES                		       | 
!   INITIALIZES: DT1(NNODE) ??                                   	       |
!==============================================================================|

   SUBROUTINE WATER_DEPTH         

!------------------------------------------------------------------------------|
   USE ALL_VARS
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP) :: TEMP
   INTEGER  :: I,K,J1,J2
!------------------------------------------------------------------------------|

!
!  ADJUST STATIC HEIGHT AND CALCULATE DYNAMIC DEPTHS (D) AND (DT)
!
   H  = H + DJUST
   D  = H + EL
   DT = H + ET

!
!  ADJUST HEIGHT ON OUTER BOUNDARY
!
   IF(IOBCN > 0) THEN
     DO I=1,IOBCN
       J1=I_OBC_N(I)
       J2=NEXT_OBC(I)
       H(J1)=H(J2)
       D(J1)=D(J2)
       DT(J1)=DT(J2)
     END DO
   END IF
          
! 
!  CALCULATE FACE-CENTERED VALUES OF BATHYMETRY AND DEPTH
!
   DO I=1,NT    
     H1(I)  = (H(NV(I,1))+H(NV(I,2))+H(NV(I,3)))/3.0_SP
     D1(I)  = H1(I)+EL1(I)
     DT1(I) = H1(I)+ET1(I)
   END DO


   RETURN
   END SUBROUTINE WATER_DEPTH
!==============================================================================|
