!==============================================================================|
!  this subroutine is used to specify the output files used for the            !
!  graphics developed by the Ocean Ecosystem Dynamics Laboratory at            !
!  SMAST/UMASSD.                                                               !
!==============================================================================|

   SUBROUTINE SECTINF             

!------------------------------------------------------------------------------|

   USE ALL_VARS

   IMPLICIT NONE
   INTEGER :: I,K,J,N1,N2
   CHARACTER(LEN=80)   :: TEMPDIR
   REAL(SP), ALLOCATABLE,DIMENSION(:,:) :: FTEMP1,FTEMP2,FTEMP3
   INTEGER , ALLOCATABLE,DIMENSION(:,:) :: NTEMP1,NTEMP2
   INTEGER , ALLOCATABLE,DIMENSION(:,:) :: NVTEMP,NBETEMP

!==============================================================================|
   
   TEMPDIR = TRIM(OUTDIR)//"/medm/"
!
!-----------------OUTPUT SURFACE NODE COORDINATES------------------------------|
!
   IF(MSR)OPEN(1,FILE=TRIM(TEMPDIR)//'xy_node.dat',STATUS='unknown') ; REWIND(1)

   IF(SERIAL)THEN
     DO I=1,M
       WRITE(1,*) VX(I)+VXMIN,VY(I)+VYMIN
     END DO
   END IF



   IF(MSR)CLOSE(1)

!
!------------------OUTPUT SURFACE ELEMENT COORDINATES--------------------------|
!

   IF(MSR)OPEN(1,FILE=TRIM(TEMPDIR)//'xy_cell.dat',STATUS='unknown') ; REWIND(1)

   IF(SERIAL)THEN
     DO I=1,N
       WRITE(1,*) XC(I)+VXMIN,YC(I)+VYMIN
     END DO
   END IF



   IF(MSR)CLOSE(1)

!
!------------------OUTPUT EDGES AND VERTICES FOR EACH ELEMENT------------------|
!
   IF(MSR)OPEN(1,FILE=TRIM(TEMPDIR)//'mesh.inf',STATUS='unknown')  ; REWIND(1)
   IF(MSR)WRITE(1,*) 'nbe(i,j),j=1,3; nv(i,j),j=1,3'

   IF(SERIAL)THEN
     DO I=1,N
       WRITE(1,100) (NBE(I,J),J=1,3),(NV(I,J),J=1,3)
     END DO
   END IF
  


   IF(MSR)CLOSE(1)

!
!------------------SHAPE FACTORS-----------------------------------------------|
!
   IF(MSR)OPEN(1,FILE=TRIM(TEMPDIR)//'shape.inf',STATUS='unknown') ; REWIND(1)
   IF(MSR)WRITE(1,*) 'au(i,j),j=1,4; av(i,j),j=1,4'

   IF(SERIAL)THEN
     DO I=1,N
       WRITE(1,200) (A1U(I,J),J=1,4),(A2U(I,J),J=1,4)
     END DO
   END IF



   IF(MSR)CLOSE(1)

!
!-------------------LINEAR FUNCTIONS-------------------------------------------|
!
   IF(MSR)OPEN(1,FILE=TRIM(TEMPDIR)//'awxcof.inf',STATUS='unknown') ; REWIND(1)
   IF(MSR)WRITE(1,*) 'aw0(i,j),awx(i,j),awy(i,j),j=1,3'

   IF(SERIAL)THEN
     DO I=1,N
       WRITE(1,300) (AW0(I,J),AWX(I,J),AWY(I,J),J=1,3)
     END DO
   END IF




   IF(MSR)CLOSE(1)

!
!------------------OUTPUT DEPTH AT NODE POINTS---------------------------------!
!
   IF(MSR)OPEN(1,FILE=TRIM(TEMPDIR)//'depth.xy',STATUS='unknown') ; REWIND(1)
   IF(MSR)WRITE(1,*) 'scat2d'
   IF(MSR)WRITE(1,500) MGL,1

   IF(SERIAL)THEN
     DO I=1,M
       WRITE(1,'(3E20.10)') VX(I)+VXMIN,VY(I)+VYMIN,H(I)
     END DO
   END IF



   IF(MSR)CLOSE(1)

!
!------------------OUTPUT SIGMA DISTRIBUTION-----------------------------------!
!
   IF(MSR)THEN
     OPEN(1,FILE=TRIM(TEMPDIR)//'sigma.dat',STATUS='unknown') ; REWIND(1)
     DO K=1,KB
       WRITE(1,'(I10,E20.10)') K,Z(K)
     END DO
     CLOSE(1)
   END IF

!
!-------------------------FORMATTING-------------------------------------------!
!

100 FORMAT(6I10)
200 FORMAT(8E18.8)
300 FORMAT(9E18.8)
500 FORMAT('xyd ',i10,' depth ',i3,' h')

   RETURN
   END SUBROUTINE SECTINF
!==============================================================================|
