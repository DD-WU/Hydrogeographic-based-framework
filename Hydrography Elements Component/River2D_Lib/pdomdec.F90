!==============================================================================|
!  SET UP LOCAL PHYSICAL DOMAIN (CONNECTIVITY/MESH)                            |
!==============================================================================|

   SUBROUTINE PDOMDEC

!==============================================================================!
   USE ALL_VARS
   USE MOD_SPHERICAL
   IMPLICIT NONE
   INTEGER I,EGL,J,IERR,I1,I2,N_SPONGE
   REAL(SP), ALLOCATABLE :: CORRG(:),CORR(:)
   REAL(SP), ALLOCATABLE :: R_SPG(:),C_SPG(:) 
   INTEGER, ALLOCATABLE  :: N_SPG(:)
   REAL(SP)  TEMP,DTMP,C_SPONGE
   INTEGER K,ITMP
   REAL(DP) VX1,VY1,VX2,VY2,VX3,VY3,EVX12,EVX13,EVX23,&
            EVY12,EVY13,EVY23,EVXY,VX12,VY12,VX23,VY23,VX31,VY31,&
	    X1_DP,Y1_DP,X2_DP,Y2_DP,DTMP_DP

!==============================================================================|
!  GENERATE LOCAL NODE CONNECTIVITY (NV) FROM GLOBAL NODE CONNECTIVITY (NVG)   |
!  USING LOCAL TO GLOBAL MAPPING FOR INTERIOR ELEMENTS (EGID)                  |
!  AND LOCAL TO GLOBAL MAPPING FOR HALO ELEMENTS (HE_LST)                      |
!==============================================================================|

   IF(SERIAL) NV = NVG



!==============================================================================|
!   SET UP LOCAL MESH (HORIZONTAL COORDINATES)                                 |
!==============================================================================|


!--------------READ IN X AND Y GLOBAL COORDINATES AT NODES---------------------!

   ALLOCATE(XG(0:MGL),YG(0:MGL)) ; XG = 0.0_SP ; YG = 0.0_SP
   DO I=1,MGL
     READ(INGRD,*)J,XG(I),YG(I)
  if (SPHERICAL==1)then
     IF(XG(I) < 0.0) XG(I) = XG(I) + 360.0
!     IF(YG(I) < 0.0) YG(I) = YG(I) + 360.0
  endif
   END DO
   CLOSE(INGRD)

!--------------CALCULATE GLOBAL MINIMUMS AND MAXIMUMS--------------------------!

  if  (SPHERICAL==1)then
   VXMIN = 0.0_SP ; VXMAX = MAXVAL(XG(1:MGL))
   VYMIN = 0.0_SP ; VYMAX = MAXVAL(YG(1:MGL))
  else
   VXMIN = MINVAL(XG(1:MGL)) ; VXMAX = MAXVAL(XG(1:MGL))
   VYMIN = MINVAL(YG(1:MGL)) ; VYMAX = MAXVAL(YG(1:MGL))
  endif

!--------------SHIFT GRID TO UPPER RIGHT CARTESIAN-----------------------------!

   XG = XG - VXMIN
   YG = YG - VYMIN
   XG(0) = 0.0_SP ; YG(0) = 0.0_SP

!--------------CALCULATE GLOBAL ELEMENT CENTER GRID COORDINATES----------------!

   ALLOCATE(XCG(0:NGL),YCG(0:NGL)) ; XCG = 0.0_SP ; YCG = 0.0_SP
   DO I=1,NGL   
     XCG(I)  = (XG(NVG(I,1)) + XG(NVG(I,2)) + XG(NVG(I,3)))/3.0_SP
     YCG(I)  = (YG(NVG(I,1)) + YG(NVG(I,2)) + YG(NVG(I,3)))/3.0_SP
   END DO
 if  (SPHERICAL==1)then
   DO I=1,NGL
     VX1=XG(NVG(I,1))
     VY1=YG(NVG(I,1))
     VX2=XG(NVG(I,2))
     VY2=YG(NVG(I,2))
     VX3=XG(NVG(I,3))
     VY3=YG(NVG(I,3))

     DO 56 K=1,1000000
!JQI< 
       EVX12=VX2-VX1
       EVX13=VX3-VX1
       EVX23=VX3-VX2

       IF(EVX12 >  180.0_SP)THEN
         EVX12 = -360.0_SP+EVX12
       ELSE IF(EVX12 < -180.0_SP)THEN
         EVX12 =  360.0_SP+EVX12
       END IF
       IF(EVX13 >  180.0_SP)THEN
	 EVX13 = -360.0_SP+EVX13
       ELSE IF(EVX13 < -180.0_SP)THEN
	 EVX13 =  360.0_SP+EVX13
       END IF
       IF(EVX23 >  180.0_SP)THEN
         EVX23 = -360.0_SP+EVX23
       ELSE IF(EVX23 < -180.0_SP)THEN
         EVX23 =  360.0_SP+EVX23
       END IF
!JQI>	    
       EVX12=ABS(EVX12)
       EVX13=ABS(EVX13)
       EVX23=ABS(EVX23)

       EVY12=ABS(VY2-VY1)
       EVY13=ABS(VY3-VY1)
       EVY23=ABS(VY3-VY2)

       EVXY=1.E-10_SP

       IF((EVX12 < EVXY) .AND.(EVX13 < EVXY) .AND. (EVX23 < EVXY) &
          .AND.(EVY12 < EVXY) .AND. (EVY13 < EVXY)                &
          .AND.(EVY23 < EVXY))THEN
         XCG(I)=VX1
         YCG(I)=VY1
         GOTO 57
       ELSE
         CALL ARCC(VX1,VY1,VX2,VY2,VX12,VY12)
         CALL ARCC(VX2,VY2,VX3,VY3,VX23,VY23)
         CALL ARCC(VX3,VY3,VX1,VY1,VX31,VY31)

         VX1=VX12
         VY1=VY12
         VX2=VX23
         VY2=VY23
         VX3=VX31
         VY3=VY31
       END IF
56     CONTINUE
57     CONTINUE
     END DO
 endif

     XCG(0) = 0.0_SP ; YCG(0) = 0.0_SP


!--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!

     IF(SERIAL)THEN
       VX = XG
       VY = YG
     END IF



!==============================================================================|
!   SET UP LOCAL MESH (BATHYMETRIC DEPTH)                                      |
!==============================================================================|

!--------------READ IN BATHYMETRY----------------------------------------------!

     ALLOCATE(HG(0:MGL))  ; HG = 0.0_SP
     DO I=1,MGL
       READ(INDEP,*) TEMP,TEMP,HG(I)
     END DO
     CLOSE(INDEP)


!--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!

     IF(SERIAL) H = HG+0.6                      !  85的转化到平海



!--------------CALCULATE EXTREMUMS---------------------------------------------!

     HMAX = MAXVAL(ABS(HG(1:MGL)))
     HMIN = MINVAL(HG(1:MGL))

!==============================================================================|
!   SET UP LOCAL CORIOLIS FORCE                                                |
!==============================================================================|

!--------------READ IN CORIOLIS PARAMETER--------------------------------------!

     ALLOCATE(CORRG(0:MGL))  ; CORRG = 0.0_SP
!  MHB:ZHAO  ADJUST FOR DIFFERENT CORIOLIS FILE FORMAT
  if  (SPHERICAL==1) then
     CORRG=YG
  else
     IF(CASENAME == "mhb")THEN
       DO I=1,MGL
         READ(INCOR,*) TEMP,CORRG(I)
       END DO
     ELSE
       DO I=1,MGL
         READ(INCOR,*) TEMP,TEMP,CORRG(I)
       END DO

     END IF
     CLOSE(INCOR)
 endif

!--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!
     ALLOCATE(CORR(0:MT)) ; CORR = 0.0_SP
     IF(SERIAL) CORR = CORRG



!==============================================================================|
!   COMPUTE FACE CENTER VALUES FOR GRID, DEPTH, AND CORIOLIS PARAMETER         |
!==============================================================================|

 if  (SPHERICAL==1)then
     IF(SERIAL) XC = XCG; YC = YCG


     COR = YC
     DO I=1,NT
       H1(I)  = SUM( H(NV(I,1:3)))/3.0_SP
       COR(I) = 2.*7.292e-5_SP*SIN(COR(I)*2.0_SP*3.14159_SP/360.0_SP)
     END DO

 else
     DO I=1,NT
!       XC(I)  = SUM(VX(NV(I,1:3)))/3.0
       XC(I)  = (VX(NV(I,1)) + VX(NV(I,2)) + VX(NV(I,3)))/3.0_SP
       YC(I)  = (VY(NV(I,1)) + VY(NV(I,2)) + VY(NV(I,3)))/3.0_SP
!       YC(I)  = SUM(VY(NV(I,1:3)))/3.0
       H1(I)  = SUM( H(NV(I,1:3)))/3.0_SP
       COR(I) = CORR(NV(I,1)) + CORR(NV(I,2)) + CORR(NV(I,3))
       COR(I) = COR(I)/3.0_SP
!       COR(I) = SUM(CORR(NV(I,1:3)))/3.0
       COR(I) = 2.*7.292e-5_SP*SIN(COR(I)*2.0_SP*3.14159_SP/360.0_SP)!自己加的？
     END DO
  endif

!==============================================================================|
!   COMPUTE SPONGE LAYER FOR OPEN BOUNDARY DAMPING                             |
!==============================================================================|

!--READ NUMBER OF SPONGE NODES AND ALLOCATE ARRAYS-----------------------------|

     READ(INSPO,*) N_SPONGE
     IF(N_SPONGE > 0 )THEN

     ALLOCATE( N_SPG(N_SPONGE) , R_SPG(N_SPONGE) , C_SPG(N_SPONGE) )

!--READ IN INDICES OF SPONGE NODES --------------------------------------------|

     DO I=1,N_SPONGE
       READ(INSPO,*) N_SPG(I),R_SPG(I),C_SPG(I)
     END DO
     CLOSE(INSPO)



!--SET SPONGE PARAMETERS-------------------------------------------------------|

     CC_SPONGE = 0.0_SP

     DO I=1,NT
       DO I1=1,N_SPONGE
         I2=N_SPG(I1)
 if (SPHERICAL==1)then
         X1_DP=XC(I)
         Y1_DP=YC(I)
         X2_DP=XG(I2)
         Y2_DP=YG(I2)
         CALL ARC(X1_DP,Y1_DP,X2_DP,Y2_DP,DTMP_DP)
         DTMP=DTMP_DP/R_SPG(I1)
 else
         DTMP=(XC(I)-XG(I2))**2+(YC(I)-YG(I2))**2
         DTMP=SQRT(DTMP)/R_SPG(I1)
 endif
         IF(DTMP <= 1.) THEN
           C_SPONGE=C_SPG(I1)*(1.-DTMP)
           CC_SPONGE(I)=MAX(C_SPONGE,CC_SPONGE(I))
         END IF
       END DO
     END DO

     DEALLOCATE(N_SPG,R_SPG,C_SPG)

   END IF !! N_SPONGE > 0

   IF(MSR)WRITE(IPT,*)'!  # SPONGE LAYER SET BY :',N_SPONGE

!==============================================================================|
!   WRITE TO SMS GRID FILE WHILE GLOBAL VALUES EXIST                           |
!==============================================================================|

   IF(MSR)THEN
     WRITE(IOSMSD,*)'scat2d'
     WRITE(IOSMSD,*)'xyd ',MGL,' dep ',1,' dep '
     DO I=1,MGL
       WRITE(IOSMSD,*) XG(I),YG(I),HG(I)
     END DO
     CLOSE(IOSMSD)
   END IF
   DEALLOCATE(CORR,CORRG)

   RETURN
   END SUBROUTINE PDOMDEC
!==============================================================================|
