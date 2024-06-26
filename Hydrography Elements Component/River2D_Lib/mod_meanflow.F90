MODULE MOD_MEANFLOW
   USE ALL_VARS
   USE MOD_PREC
   USE MOD_TYPES
   IMPLICIT NONE
   SAVE
   INTEGER              :: INMF,INTCELL,INTNODE,INTELEL,INTUV
   INTEGER              :: nmfcell_GL, nmfcell
   INTEGER, ALLOCATABLE :: MF_GL2LOC(:)
   INTEGER, ALLOCATABLE :: I_MFCELL_GL(:),I_MFCELL_N(:)
   REAL(SP),ALLOCATABLE :: DMFQDIS(:,:),MFQDIS(:),MFDIST(:,:)
   REAL(SP),ALLOCATABLE :: ANGLEMF(:),MFAREA(:),VLCTYMF(:)
   TYPE(BC)             :: MF_TM           !!TIME MAP FOR MEAN FLOW DATA
   REAL(SP),ALLOCATABLE :: RDISMF(:,:)
   INTEGER ,ALLOCATABLE :: NODE_MFCELL(:,:)

   CONTAINS

! we still need to consider the case in which MEAN FLOW bring in/take out T & S
!==============================================================================|
!  READ IN MEAN FLOW OPEN BOUNDARY FLUX (m^3/s^1) TIME SERIES                  |
!==============================================================================|

   SUBROUTINE READ_MEANFLOW

!------------------------------------------------------------------------------!
     INTEGER              :: k,i,j,i1,i2,i3,ii,NCNT,itemp,IERR
     INTEGER, ALLOCATABLE :: temp1(:),temp2(:)
     REAL(SP),ALLOCATABLE :: RTEMP1(:,:),RTEMP2(:,:)
     REAL(SP)             :: ttemp

     REWIND(INMF)
     READ(INMF,*) nmfcell_GL

     nmfcell = 0
  IF (nmfcell_GL > 0) THEN

     ALLOCATE(I_MFCELL_GL(nmfcell_GL))
     DO I=1,nmfcell_GL
        READ(INMF,*)I_MFCELL_GL(I)
     ENDDO

!----Read in Mean Flow Flux Vertical Distribution---------------------
     ALLOCATE(RTEMP1(nmfcell_GL,KBM1))
     DO I = 1, nmfcell_GL
       READ(INMF,*) J,(RTEMP1(I,K),K = 1,KBM1)
     END DO

!----Read in Time Dependent DataSets ---------------------------------
       READ(INMF,*) itemp
       MF_TM%NTIMES = itemp
       MF_TM%LABEL  = "open boundary mean flow flux"
       ALLOCATE(MF_TM%TIMES(itemp))
       ALLOCATE(RTEMP2(nmfcell_GL,itemp))
       DO I = 1, itemp
         READ(INMF,*) ttemp
         MF_TM%TIMES(I) = ttemp
         READ(INMF,*) (RTEMP2(J,I),J = 1,nmfcell_GL)
         WRITE(IOPRT,*) ttemp
         WRITE(IOPRT,*) (RTEMP2(J,I),J = 1,nmfcell_GL)
       END DO
       CLOSE(INMF)

!
!---Map to Local Domain----------------------------------------

     IF(SERIAL)THEN
       nmfcell = nmfcell_GL
       ALLOCATE(I_MFCELL_N(nmfcell))
       I_MFCELL_N = I_MFCELL_GL
       ALLOCATE(MFDIST(nmfcell,kbm1))
       MFDIST = RTEMP1
       ALLOCATE(DMFQDIS(nmfcell,MF_TM%NTIMES))
       DMFQDIS = RTEMP2
     END IF



     DEALLOCATE(RTEMP1,RTEMP2)

  ELSE  ! if statement end for nmfcell_GL > 0
    close(INMF)
  END IF




   RETURN
   END SUBROUTINE READ_MEANFLOW
!==============================================================================|


!==============================================================================|
!  SET METRICS FOR MEAN FLOW BOUNDARY CONDITIONS       			       |
!==============================================================================|

   SUBROUTINE SET_BNDRY_MEANFLOW     

!------------------------------------------------------------------------------!

   USE BCS
   USE MOD_SPHERICAL

   USE MOD_OBCS

   IMPLICIT NONE
   REAL(DP)  DX12,DY12,ATMP1,HTMP
   INTEGER I,J,I1,I2,J1,J2,II,ITMP,JTMP
   REAL(DP) X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE

!------------------------------------------------------------------------------!

   IF(nmfcell > 0)THEN

     ALLOCATE(ANGLEMF(nmfcell),MFAREA(nmfcell),VLCTYMF(nmfcell),MFQDIS(nmfcell))
     ALLOCATE(NODE_MFCELL(nmfcell,2),RDISMF(nmfcell,2))

     DO I=1,nmfcell
       II=I_MFCELL_N(I)
       IF(ISBCE(II) /= 2) THEN
          PRINT*, 'NO.',I,'MEAN FLOW CELL'
          PRINT*, 'IS NOT A OPEN BOUNDARY ONE'
          CALL PSTOP
       END IF
       ITMP=0
       DO J=1,3
         IF(NBE(II,J) == 0) THEN
           JTMP=J
           ITMP=ITMP+1
         END IF
       END DO
       IF(ITMP /= 1) THEN
         PRINT*, 'NO OPEN BOUNDARY OR MORE THAN ONE OPEN BOUNDARY'
         PRINT*, 'IN NO.',I,'MEAN FLOW CELL'
         CALL PSTOP
       END IF
       J1=JTMP+1-INT((JTMP+1)/4)*3
       J2=JTMP+2-INT((JTMP+2)/4)*3
       I1=NV(II,J1)
       I2=NV(II,J2)
         
       NODE_MFCELL(I,1)=I1
       NODE_MFCELL(I,2)=I2

       HTMP=0.5_SP*(H(I1)+H(I2))     ! may be a problem here, should be replaced dy D
       DY12=VY(I1)-VY(I2)
     if (SPHERICAL==1)THEN
       X1_DP = VX(I2)
       Y1_DP = VY(I2)
       X2_DP = VX(I1)
       Y2_DP = VY(I1)
       CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
       DX12 = SIDE

       DY12 = TPI*DY12
      else
       DX12=VX(I1)-VX(I2)
      endif
       ATMP1=ATAN2(DY12,DX12)
       MFAREA(I)=SQRT(DX12**2+DY12**2)*HTMP    ! for spherical coordinates is Phthagolean Theorem still valid?
       ANGLEMF(I)=ATMP1+3.1415927/2.0
       RDISMF(I,1)=ART1(I1)/(ART1(I1)+ART1(I2))
       RDISMF(I,2)=ART1(I2)/(ART1(I1)+ART1(I2))
     END DO
   END IF

   RETURN
   END SUBROUTINE SET_BNDRY_MEANFLOW
!==============================================================================|

!==============================================================================|
!  INTERPOLATION MEAN FLOW OPEN BOUNDARY FLUX (m^3/s^1) TIME SERIES            |
!==============================================================================|

   SUBROUTINE BCOND_MEANFLOW
!
!------------------------------------------------------------------------------!

   USE BCS
   USE MOD_OBCS

   INTEGER  L1,L2,IERR
   REAL(SP) :: FACT,UFACT

   IF(nmfcell > 0)THEN
     CALL BRACKET(MF_TM,THOUR,L1,L2,FACT,UFACT,IERR)
     MFQDIS(:) = UFACT*DMFQDIS(:,L1) + FACT*DMFQDIS(:,L2)
     MFQDIS    = MFQDIS*RAMP
   END IF

   RETURN
   END SUBROUTINE BCOND_MEANFLOW

END MODULE MOD_MEANFLOW
