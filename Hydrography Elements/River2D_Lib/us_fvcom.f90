!==============================================================================!
!  VERSION 2.6 
!==============================================================================!

    PROGRAM US_FVCOM

!==============================================================================!
!                                                                              !
!                            USG-FVCOM                                         !
!    The Unstructured Grid Finite Volume Coastal Ocean Model                   !
!                                                                              !
!    The USG-FVCOM (publically called FVCOM) was developed by Drs. Changsheng  !
!  Chen and Hedong Liu at the Marine Ecosystem Dynamics Modeling Laboratory    !
!  at the School of Marine Science and Technology (SMAST), University of       !
!  Massachusetts-Dartmouth (UMASSD) and Dr. Robert C. Beardsley at the         !
!  Department of Physical Oceanography, Woods Hole Oceanographic Institution   !
!  (WHOI). This code was rewritten in Fortran 90/2K, modularized, rendered     ! 
!  somewhat understandable, and parallelized by Geoff Cowles at SMAST/UMASSD.  !
!                                                                              !
!    The Development was initially supported by the Georgia Sea Grant College  !
!  Program for the study of the complex dynamic system in Georgia estuaries.   !
!  The code improvement has been supported by Dr. Chen's research grants       !
!  received from NSF and NOAA Coastal Ocean research grants and SMAST-NASA     !
!  fishery research grants. The University of Massachusetts-Dartmouth          !
!                                                                              !
!    FVCOM is a three dimensional,time dependent,primitive equations,          !
!  coastal ocean circulation model. The model computes the momentum,           !
!  continuity, temperature, salinity, and density equations and is closed      !
!  physically and mathematically using the Mellor and Yamada level-2.5         !
!  turbulent closure submodel. The irregular bottom slope is represented       !
!  using a sigma-coordinate transformation, and the horizontal grids           !
!  comprise unstructured triangular cells. The finite-volume method (FVM)      !
!  used in this model combines the advantages of the finite-element            !
!  method (FEM) for geometric flexibility and the finite-difference            !
!  method (FDM) for simple discrete computation. Current, temperature,         !
!  and salinity in the model are computed in the integral form of the          !
!  equations, which provides a better representation of the conservation       !
!  laws for mass, momentum, and heat in the coastal region with complex        !
!  geometry.                                                                   !
!                                                                              !
!    All users should read this agreement carefully.  A user, who receives any !  
!  version of the source code of FVCOM, must accept all the terms and          !
!  conditions of this agreement and also agree that this agreement is like any !
!  written negotiated agreement signed by you. You may be required to have     !
!  another written agreement directly with Dr. Changsheng Chen at SMAST/UMASS-D!
!  or Dr. Brian Rothschild, Director of the SMAST/UMASS-D that supplements     !
!  all or portions of this agreement. Dr. Changsheng Chen, leader of the       !
!  FVCOM development team, owns all intellectual property rights to the        !
!  software. The University of Massachusetts-Dartmouth and the Georgia Sea     !
!  Grant Program share the copyright of the software. All copyrights are       !
!  reserved. Unauthorized reproduction and re-distribution of this program     !
!  are expressly prohibited. This program is only permitted for use in         !
!  non-commercial academic research and education.  Commercial use must be     !
!  approved by Dr. Chen and is subject to a license fee. Registration is       !
!  required for all new users.  Users should realize that this model software  !
!  is a research product without any warranty. Users must take full            !
!  responsibility for any mistakes that might be caused by any inappropriate   !
!  modification of the source code.  Modification is not encouraged for users  !
!  who do not have a deep understanding of the finite-volume numerical methods !
!  used in FVCOM. Contributions made to correcting and modifying the programs  !
!  will be credited, but will not affect copyrights. No duplicate              !
!  configurations of FVCOM are allowed in the same geographical region,        !
!  particularly in the regions where FVCOM has been already been applied.      !
!  Users who want to use FVCOM in a region that the SMAST/UMASS Marine         !
!  Ecosystem Dynamics Modeling (MEDM) group (led by Dr. Chen) is working on    !
!  must request permission from Dr. Chen. No competition is allowed in the     !
!  same region using FVCOM, especially with Dr. Chen脙聲s group. FVCOM has been   !
!  validated for many standard model test cases.  Users are welcome to do any  !
!  further model validation experiments. These experiments shall be carried    !
!  out in collaboration with the SMAST/UMASSD model development team. To avoid !
!  or reduce deriving any incorrect conclusions due to an inappropriate use of !
!  FVCOM, users are required to contact the scientific leaders of the FVCOM    !
!  development team (Dr. Chen at SMAST/UMASS-D and Dr. Beardsley at WHOI)      !
!  before any formal publications are prepared for model validation.           !
!                                                                              !
!    For public use, all users should name this model as "FVCOM". In any       !
!  publications with the use of FVCOM, acknowledgement must be included. The   !
!  rationale behind this FVCOM distribution policy is straightforward.  New    !
!  researchers and educators who want to use FVCOM and agree to the above      !
!  requirements get free access to the latest version of FVCOM and the         !
!  collective wisdom and experience of the FVCOM development team and existing !
!  users. Problems arising in new FVCOM applications, both related to          !
!  conceptual as well as numerical and coding issues, can be shared with the   !
!  development team and other users who can work together on physics and code  !
!  improvements that over time will lead to a better FVCOM.                    !
!                                                                              !
!    FVCOM has been developed to date with state and federal funding with the  !
!  idea that FVCOM will become a 脙聮community脙聯 model as new users start to use   !
!  the model and its scientific usefulness and numerical accuracy and          !
!  efficiency continue to improve.  The FVCOM distribution policy is designed  !
!  to encourage this transition while maintaining a central core group         !
!  responsible for overall FVCOM development direction, implementing 脙聮official脙聯!
!  code improvements, and maintaining well tested and documented updated code  !
!  versions.                                                                   !       
!                                                                              !
!                                                                              !                                                                           !
!  External forces used to drive this model:                                   !
!                                                                              !
!  1) Tidal amplitudes and phases at open boundaries (initial designs          !
!         include 6 tidal consituents, more can be added as needed);           !
!  2) Wind Stress [3 ways: a) uniform wind speed and direction, b)             !
!         spatially distributed wind velocity field, and c)the MM5 model-out   !
!         wind fields]                                                         !
!  3) Surface heat flux [3 ways: a) uniform heat flux, b) spatially            !
!         distributed heat flux, and c) the MM5-output heat flux fields        !
!         All the surface net heat flux and short-wave radiation are needed    !
!         in the input file                                                    ! 
!  4) River discharges: specify the location and discharge volume,             !
!         temperature, and salinity                                            !                      
!  5) Groundwater input: currently diffused bottom flux only                   !
!                                                                              !
!  Initial conditions:                                                         !
!                                                                              !
!  The model can be prognostically run for both barotropic and baroclinic      !
!  cases.                                                                      !
!                                                                              !
!  Tidal forcing can be added into the system with zero velocity               !
!  field at initial or specified the 3-D tidal initial velocity field          !
!  using the model-predicted harmonic tidal currents.                          !
!                                                                              !
!  Initial fields of temperature and salinity needed to be specified           !
!  by using either climatological field, real-time observed field or           !
!  idealized functions. The model has included Gregorian time for the          !
!  time simulation for tidal currents.                                         !
!                                                                              !
!  For the purpose of interdisciplinary studies, biological, chemical, and     !
!  sediment suspension models are available for FVCOM.  These submodels are    !
!  directly driven by the FVCOM physical model. A description of these         !
!  submodels follows.                                                          !
!                                                                              !
!  Generalized biological modules-a software platform that allows users to     !                                                        !              
!             build his own biological model in FVCOM                          !
!                                                                              !
!  NPZ model--a 3 component nutrient-phytoplankton-zooplankton model           !
!                                                                              !
!  NPZD model--an 8 component nutrient-phytolankton-zooplankton-detritus       !
!              model;                                                          !
!                                                                              !
!  NPZDB-model-a 9 phosphorus-controlled component nutrient-                   !
!               phytoplankton-zooplankton-detritus-bacteria model;             !
!                                                                              !
!  Water quality model with inclusion of the benthic process                   !
!                                                                              !
!  Sediment model--A new module that was developed by Dr. Geoff Cowles         !
!                                                                              !
!  Lagrangian particle tracking:                                               !
!                                                                              !
!  A bilinear interpolation scheme is used to determine the particle           !
!  velocity for the Lagrangian particle tracking. A random walk process        !
!  also could be included with a specified function related to horizontal      !
!  and vertical diffusion coefficients                                         !
!                                                                              !
!  Key reference:                                                              !
!                                                                              !
!   Chen, C., H. Liu, and R. C. Beardsley, 2003. An unstructured grid,         !
!       finite-volume, three-dimensional, primitive equations ocean            !
!       model: application to coastal ocean and estuaries, Journal             !
!       of Atmospheric and Oceanic Technology,  20, 159-186.                   !
!                                                                              !
!                                                                              !
!                                                                              !
!  Please direct criticisms and suggestions to                                 !
!                                                                              !
!               Changsheng Chen                                                !
!               School for Marine Science and Technology                       !
!               University of Massachusetts-Dartmouth                          !
!               New Bedford, MA 02742                                          !
!               Phone: 508-910-6388, Fax: 508-910-6371                         !
!               E-mail: c1chen@umassd.edu                                      !
!               Web: http://fvcom.smast.umassd.edu                             !
!                                                                              !
! What are new for version 2.5?                                                !
!    1) A new spherical coordinate version is added with an accurate treatment !
!       of the north pole (Arctic Ocedan) by Dr. Chen et al. at UMASSD         !
!    2) Spherical and x-y coordinate versions was merged into a single code    !
!       with users' choice for either coordiante for their application         !
!    3) Multiple choices of open radiation boundary conditions are added       !
!    4) General turbulence modules were linked by Dr. Cowles                   !
!    5) The selection for a 2-D barotrophic application is added               !
!    6) bugs in paralleziation and wet/dry point treatments are corrected      !
! For more detailed information, please read the upgrade user manaul           !
!                                                                              !
! What will be included in version 2.4 (under validation tests)                !
!    1) Generalized Biological Modules (developed by Dr. Tian and Chen)        !
!    2) 3-D sediment model(developed by Dr. Cowles)                            !
!    3) Reduced Kalman Filter, Ensemble Kalman Filter, and Ensemble Tranistion !
!       Kalman Filter (implented into FVCOM with UMASSD and MIT groups led     !
!       by Chen and Roziili                                                    !
!    4) Full nonlinear ice models (developed orignally to FVCOM by Dr.Dupont   !
!       and modified to be implemented to a parallelized new spherical         !
!       coordinate version by the MEDM group at UMASSD                         !
!                                                                              !   
! Enjoy!                                                                       !
!==============================================================================!

!==============================================================================!
!  INCLUDE MODULES                                                             !
!==============================================================================!
    USE ALL_VARS
    USE MOD_CLOCK
    !  if (DATA_ASSIM==1)USE MOD_ASSIM
    USE MOD_SPHERICAL
    USE MOD_WD
!   if (NETCDF_IO==1)then
!   USE MOD_NCDIO
!   USE MOD_NCDAVE
!   endif
    USE BCS
    USE PROBES  
    USE MOD_OBCS   
    !            USE MOD_LAG           NETCDF_IO
    !  if(GOTM==1)
    !  USE MOD_GOTM
    !  endif
!  if (MEAN_FLOW==1)then
    USE MOD_MEANFLOW
    USE MOD_OBCS2
    USE MOD_OBCS3
!  endif
    USE MOD_WQM
    USE MOD_NORTHPOLE
    USE MOD_EQUITIDE
!  SEDIMENT BEGIN
    USE MOD_SED
    USE MOD_DYE
    USE MOD_TSOBC
    USE MOD_BALANCE_2D
    use wave
!QXU}
!------------------------------------------------------------------------------|
    IMPLICIT NONE
    LOGICAL  :: FEXIST
    REAL(sp) :: TMP1,TMP,UTMP,VTMP,TTIME,  uvmax
    INTEGER  :: I,K,J,IERR,N1,J1,J2,JN
    REAL(sp), ALLOCATABLE :: FTEMP(:),FTEMP2(:)
    CHARACTER(LEN=13) :: TSTRING
!------------------------------------------------------------------------------|
! New Open Boundary Condition ----1
!  SEDIMENT END
!JQI<added on 07/14/04
  
!>JQI 07/14/04
!QXU{
!==============================================================================!
!  FVCOM VERSION                                                               !
!==============================================================================!

    FVCOM_VERSION     = 'FVCOM_2.6'
    FVCOM_WEBSITE     = 'http://fvcom.smast.umassd.edu, http://codfish.smast.umassd.edu'
!#  if defined (NETCDF_IO)
!   INSTITUTION       = 'School for Marine Science and Technology'
!   NETCDF_TIMESTRING = 'seconds after 00:00:00'
! #  endif
!==============================================================================!
!   SETUP PARALLEL ENVIRONMENT                                                 !
!==============================================================================!
    SERIAL = .TRUE. 
    PAR    = .FALSE. 
    MSR    = .TRUE.
    MYID   = 1
    NPROCS = 1
!    CALL GET_CASENAME(CASENAME)
!-----------------------------------------------------
!  giving then casename
!-----------------------------------------------------------
    INPDIR="..\run\"
    OUTDIR="..\output\"
    casename='zjk'
    casename1= trim(INPDIR)//trim(casename)//"_run.dat" 
   
    INQUIRE(FILE=casename1,EXIST=FEXIST)
        IF(.NOT.FEXIST)THEN
        WRITE(IPT,*)'FILE ',TRIM(casename)//"_run.dat",' DOES NOT EXIST'
        WRITE(IPT,*)'STOPPING...'
        CALL PSTOP
        END IF
   
!==============================================================================!
!   SETUP MODEL RUN                                                            !
!==============================================================================!
!--------------------------------------------------------------------------
!        SPHERICAL          SELECT SPHERICAL COORDINATES FOR INTEGRATION
!                           DEFAULT: CARTESIAN
!                           UNCOMMENT TO SELECT SPHERICAL COORDINATES
!--------------------------------------------------------------------------
    spherical=0
!--------------------------------------------------------------------------
!        FLOODYING/DRYING   INCLUDE WET/DRY TREATMENT OF DOMAIN
!                           CAN BE ACTIVATED/DEACTIVATED AT RUN TIME WITH
!                           INPUT FILE CONTROL.  (SEE exa_run.dat) FILE
!                           DEFAULT: NO FLOODYING/DRYING INCLUDED
!                           UNCOMMENT TO INCLUDE FLOODYING/DRYING
!--------------------------------------------------------------------------
	    WET_DRY=1
!--------------------------------------------------------------------------
!        WATER_QUALITY      INCLUDE EPA WATER QUALITY MOD
!                           CAN BE ACTIVATED/DEACTIVATED AT RUN TIME WITH
!                           VARIABLE WQM_ON IN INPUT FILE
!                           DEFAULT: NO WATER QUALITY MODEL
!                           UNCOMMENT TO INCLUDE WATER QUALITY MODEL
!--------------------------------------------------------------------------
	    WATER_QUALITY=0
!--------------------------------------------------------------------------
!        DATA_ASSIMILATION  INCLUDE NUDGING BASED DATA ASSIMILATION FOR
!                           CURRENT/TEMP/SALINITY/SST
!                           CAN BE ACTIVATED/DEACTIVATED AT RUN TIME WITH 
!                           INPUT FILE CONTROL.  (SEE exa_run.dat) FILE 
!                           DEFAULT: NO DATA ASSIMILATION INCLUDED 
!                           UNCOMMENT TO INCLUDE DATA ASSIMILATION 
!--------------------------------------------------------------------------
	DATA_ASSIM=0
!------------------------------------------------------------------------------
!       SOLID BOUNDARY     IF GCN, NO GHOST CELL
!                          IF GCY1, GHOST CELL IS SYMMETRIC RELATIVE TO BOUNDARY
!                         CELL EDGE
!                          IF GCY2, GHOST CELL IS SYMMETRIC RELATIVE TO MIDDLE
!                          POINT OF THE BOUNDARY CELL EDGE
!   !!!!!! ONLY ONE OF THE FLAGS BELOW CAN BE AND MUST BE CHOSEN
!---------------------------------------------------------------------------
	GCN=1
    GCY1=0
    GCY2=0
!---------------------------------------------------------------------
    EQUI_TIDE=0
    ATMO_TIDE=0
!--------------------------------------------------------------------------
!       ARCTIC OCEAN INCLUDED (If you chose this flag, FLAG_2 should be 
!                              selected)
!--------------------------------------------------------------------------
    NORTHPOLE=0
!--------------------------------------------------------------------------
!        Using A fully multidimensional positive definite advection
!        transport algorithm with small implicit diffusion. 
!        Based on Smolarkiewicz, P. K; Journal of Computational
!        Physics, 54, 325-362, 1984
!--------------------------------------------------------------------------	
	MPDATA=0
!--------------------------------------------------------------------------
!         Run Two-D Barotropic Mode Only
!--------------------------------------------------------------------------
    TWO_D_MODEL=1
!--------------------------------------------------------------------------
!         Output 2-D Momentum Balance Checking
!--------------------------------------------------------------------------  
	BALANCE_2D=0
!--------------------------------------------------------------------------
!           open boundary T/S time series nudging
!---------------------------------------------------------------------------
	TS_OBC=0
!--------------------------------------------------------------------------
!           OPEN BOUNDARY FORCING TYPE 
!           DEFAULT: OPEN BOUNDARY NODE WATER ELEVATION FORCING
!           UNCOMMENT TO SELECT BOTH OPEN BOUNDARY NODE WATER ELEVATION
!           FORCING AND OPEN BOUNDARY VOLUME TRANSPORT FORCING 
!---------------------------------------------------------------------------
	MEAN_FLOW=0
!--------------------------------------------------------------------------
!           dye release                                                                                                                            
!---------------------------------------------------------------------------  
    DYE_RELEASE=0
	sediment=0
    iwave=0
    icyclone=0
!
!  READ PARAMETERS CONTROLLING MODEL RUN
!
    CALL DATA_RUN
!
!  READ PARAMETERS CONTROLLING OBC FOR TEMPERATURE AND SALINITY
!
    CALL TSOBC_TYPE
!
!  READ PARAMETERS CONTROLLING DATA ASSIMILATION
!
!  if (DATA_ASSIM==1)   CALL SET_ASSIM_PARAM
!
!  READ PARAMETERS CONTROLLING WET_DRY TREATMENT
!
    CALL SET_WD_PARAM
!QXU{
    if (BALANCE_2D==1)then
    CALL SET_BALANCE_PARAM
    endif
!
!  READ PARAMETERS CONTROLLING DYE RELEASE
!
    if (DYE_RELEASE==1)then
    CALL SET_DYE_PARAM
    endif   
!
!  READ PARAMETERS CONTROLLING T/S OBC Series NUDGING
!
    if  (TS_OBC==1)then
    CALL SET_TSOBC_PARAM
    endif   
!QXU}
!  OPEN INPUT/OUTPUT FILES
    CALL IOFILES
!
!  DETERMINE NUMBER OF ELEMENTS AND NODES IN THE MODEL
!
    CALL GETDIM
!
!  READ PARAMETERS CONTROLLING WATER QUALITY MODEL 
!
    if (WATER_QUALITY==1)then
    CALL GET_WQMPAR
    endif
!
!  READ PARAMETERS CONTROLLING NETCDF OUTPUT
!
!
!  DECOMPOSE DOMAIN BY ELEMENTS USING METIS
!
!
!  MAP OPEN BOUNDARY CONDITION NODES TO LOCAL DOMAIN
!
    CALL BCMAP           !  开边界位置，类型
!
!  INPUT AND SETUP BOUNDARY FORCING (HEAT/RIVERS/WIND/etc)
!
    CALL BCS_FORCE           ! 边界条件数据
!
!  INPUT WATER QUALITY MODEL VARIABLES
!
    if (WATER_QUALITY==1)then
    IF(WQM_ON)CALL BCS_FORCE_WQM
    endif
!QXU{
!
!  INPUT T/S OBC Series NUDGING VARIABLES
!
    IF(TSOBC_ON) CALL READ_TSOBC
!QXU}
!
!  ALLOCATE FLOWFIELD VARIABLES
!
    CALL ALLOC_VARS
!
!  ALLOCATE ELM1 AND ELM2 FOR ORLANSKI RADIATION OPEN BOUNDARY CONDITION
!
    CALL ALLOC_OBC_DATA

!
!  ALLOCATE AND WET/DRY CONTROL ARRAYS
!

    CALL ALLOC_WD_DATA

!
!  ALLOCATE WATER QUALITY MODEL VARIABLES
!
    if (WATER_QUALITY==1)then
    IF(WQM_ON)CALL ALLOC_WQM_VARS
    endif

!
!  ALLOCATE SPHERICAL COORDINATE SYSTEM VARS
!
    if (SPHERICAL==1)then
    CALL ALLOC_SPHERE_VARS
    endif

!QXU{
!
!  ALLOCATE MOMENTUM BALANCE CHEKING VARS
!
    if  (BALANCE_2D==1)then
    CALL ALLOC_BALANCE_VARS
    endif
!
!  ALLOCATE DYE VARS
!
    if (DYE_RELEASE==1)then
    CALL ALLOC_VARS_DYE
    endif
!QXU}

!
!  SHIFT GRID/CORIOLIS/BATHYMETRY TO LOCAL DOMAIN
!
    CALL PDOMDEC      ! 读网格

!
!  ALLOCATE EQUILIBRIUM AND ATMOSPHERIC TIDE VARS
!
!JQI<
    if (EQUI_TIDE==1)then
    CALL ALLOCATE_EQUI
    endif  


!
!  SET UP GRID METRICS (FLUX EDGES/CONTROL VOLUMES/ETC)
!
    CALL TRIANGLE_GRID_EDGE      !Set up fluxes and control Volumes
    CALL SET_SIGMA               !Build Sigma Coordinate
    CALL CELL_AREA               !Calculate Element and Control Volume Areas
    if (GCN==1)then
    CALL SHAPE_COEF_GCN          !Calc Shape Coefficients for Flux Construction
    else
    CALL SHAPE_COEF_GCY          !Calc Shape Coefficients for Flux Construction
    endif

!JQI<for north pole
    if (SPHERICAL==1)then
    if (NORTHPOLE==1)then
    CALL SHAPE_COEF_XY
!  FIND THE NODE NUMBER OF NORTH POLE IF EXIST
    CALL FIND_NORTHPOLE
    CALL FIND_CELLSIDE
    endif   
    endif
!JQI>


!      ---------- by  zhangzhuo************************************
	    if (icyclone==1)then
        call  cyclone 
        end if 

        if (iwave==1) then
        call ALLOC_VARS_WAVE

		call readwave(1)

!		do i=1, 110
!        call readwave(2)
!        enddo
		call wavechazhi
        call wavepump
	
        endif
!
!  SETUP OPEN BOUNDARY METRICS                                   
!
    CALL SETUP_OBC 

    CALL SET_BNDRY               !Boundary Condition Metrics  inflow condtion
!
!  INITIALIZE FLOWFIELD  --> [T,S,U,V,EL,D,Q2,Q2L]
!
    CALL STARTUP
!
!  INITIALIZE GOTM
!

!
!  INITIALIZE SEDIMENT MODEL
!
    if(sediment==1) then
    IF(SEDIMENT_ON)CALL SETUP_SED(RESTART_SED)
    endif   
    !   call tecplot
    !  stop
!
!  CALCULATE DEPTH HORIZONTAL DERIVATIVES
!
    CALL DEPTH_GRADIENT
! 
!  GENERATE SECTION INFORMATION FILES
!
!   CALL SECTINF
!

    !  call tecplot

!

!
!  INITIALIZE VISUALIZATION SERVER
!

!  INITIALIZE DATA ASSIMILATION VARIABLES
!
!  if  (DATA_ASSIM==1)then
!   IF(CURRENT_ASSIM) CALL SET_CUR_ASSIM_DATA
!   IF(SST_ASSIM)     CALL SET_SST_ASSIM_DATA
!   IF(TS_ASSIM)      CALL SET_TS_ASSIM_DATA
!  endif
!
!  INITIALIZE TIME SERIES OBJECTS 
!
    !  CALL SET_PROBES      

!
!  REPORT STATISTICS ON INITIAL VALUES
!
    CALL REPORT('INITIAL VALUE INFORMATION')

!
!  CALCULATE INTERNAL TIME STEP AND SET INTEGRATION LIMITS
!
    ISTART=IINT+1

    CALL START_CLOCK

!
!  REPORT INTEGRATION INITIAL TIMES 
!
    IF(MSR)CALL REPORT_SIMTIME

!
!  CALCULATE INTERVALS FOR SST DATA ASSIMILATION 
!
!  if (DATA_ASSIM==1)then
!   CALL SET_ASSIM_INTERVALS
! endif

 
       
! New Open Boundary Condition ----2
    if  (MEAN_FLOW==1)then
        CALL FIND_OBSIDE
        CALL ALLOC_OBC3_DATA
        CALL SETUP_OBC3
        CALL ALLOC_OBC2_DATA
    endif


    open(112, file=TRIM(OUTDIR)//"/"//trim(casename)//'_topo.dat')
    !  BEGIN MAIN LOOP OVER PHYSICAL TIME STEPS                                    |
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DO IINT=ISTART,IEND
            TIME =DTI*FLOAT(IINT)/86400.0_SP
            THOUR=DTI*FLOAT(IINT)/3600.0_SP
        !---------------------------renew wave parameter------------------------------------------------------
            if(iwave==1.and.thour>=27.and.thour<72.and.mod(int(thour*3600),900)==0) then
                call readwave(2)                       ! 这块不具有热启动功能   thour>=27
	        if(mod(iint,360)==0)then
	            call wavechazhi
        !	 call wavestress   
	            call wavepump
	            endif
    
            endif

        !----SET RAMP FACTOR TO EASE SPINUP--------------------------------------------!
            RAMP=1.0_SP
            IF(IRAMP /= 0) RAMP=TANH(FLOAT(IINT)/FLOAT(IRAMP))
        !    IF(IRAMP /= 0) RAMP=TANH(FLOAT(IINT-ISTART+1)/FLOAT(IRAMP))

        !----SET UP WATER QUALITY MODEL COEFFICIENTS-----------------------------------!
            if  (WATER_QUALITY==1)then
                IF(WQM_ON)THEN
                TIME_R=MOD(IINT*DTI/3600.0_SP-14.743_SP, 24.0_SP)+6.0_SP
                CALL WQMCONST
                END IF
            endif



            if  (TWO_D_MODEL==0)then
                CALL ADJUST2D3D(1)

        !----SPECIFY THE SOLID BOUNDARY CONDITIONS OF U&V INTERNAL MODES---------------!
            if  (GCN==1)then
                CALL BCOND_GCN(5)
            else
                CALL BCOND_GCY(5)
            endif


            !----SPECIFY THE SURFACE FORCING OF INTERNAL MODES-----------------------------!
            if  (GCN==1)then
                CALL BCOND_GCN(8)
            else
                CALL BCOND_GCY(8)
            endif

        ! New Open Boundary Condition ----3
            if (MEAN_FLOW==1)then
                CALL BCOND_MEANFLOW
                CALL BCOND_TIDE_3D
        !     CALL BCOND_NG_3D       !change BKI
        !     CALL BCOND_NG_2D       !change BKI
                CALL BCOND_BKI_3D(1)  
        !     CALL BCOND_BKI_2D(4)

                CALL FLUX_OBC3D_2
                CALL FLUX_OBC2D
            endif
            endif    !end defined (TWO_D_MODEL)
        !----SPECIFY BOTTOM FRESH WATER INPUT BOUNDARY CONDITION-----------------------!
                CALL BCOND_BFW
        !----SPECIFY THE BOTTOM ROUGHNESS AND CALCULATE THE BOTTOM STRESSES------------!
                CALL BOTTOM_ROUGHNESS
        !
        !==============================================================================!
        !  CALCULATE DISPERSION (GX/GY) AND BAROCLINIC PRESSURE GRADIENT TERMS         !
        !==============================================================================!
            if  (TWO_D_MODEL==0)then

            if  (GCN==1)then
                CALL ADVECTION_EDGE_GCN(ADVX,ADVY)          !Calculate 3-D Adv/Diff       !
            else
                CALL ADVECTION_EDGE_GCY(ADVX,ADVY)          !Calculate 3-D Adv/Diff       !
            endif
                IF(.NOT. BAROTROPIC)THEN                    !Barotropic Flow ?            !
        !QXU{
        !Calculate the rho mean
                IF(IRHO_MEAN > 0) THEN
                IF(MOD(IINT,IRHO_MEAN) == 0) CALL RHO_MEAN    
                END IF
        !QXU}     
                IF(C_BAROPG == 'sigma')     CALL BAROPG      !Sigma Level Pressure Gradient!
                IF(C_BAROPG == 's_levels') CALL PHY_BAROPG  !Z Level Pressure Gradient    !
                END IF                                      !                             !
                                                            !                             !
                ADX2D = 0.0_SP ; ADY2D = 0.0_SP             !Initialize GX/GY Terms       !
                DRX2D = 0.0_SP ; DRY2D = 0.0_SP             !Initialize BCPG for Ext Mode !

                DO K=1,KBM1
                DO I=1, N
                    ADX2D(I)=ADX2D(I)+ADVX(I,K)*DZ(K)
                    ADY2D(I)=ADY2D(I)+ADVY(I,K)*DZ(K)
                    DRX2D(I)=DRX2D(I)+DRHOX(I,K)*DZ(K)
                    DRY2D(I)=DRY2D(I)+DRHOY(I,K)*DZ(K)
                END DO
                END DO
                if  (GCN==1)then
                CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
            else
                CALL ADVAVE_EDGE_GCY(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
            endif
                ADX2D = ADX2D - ADVUA                       !Subtract to Form GX
                ADY2D = ADY2D - ADVVA                       !Subtract to Form GY
        !----INITIALIZE ARRAYS USED TO CALCULATE AVERAGE UA/E  OVER EXTERNAL STEPS-----!
                UARD = 0.0_SP
                VARD = 0.0_SP
                EGF  = 0.0_SP
            if (EQUI_TIDE==1) then
                EGF_EQI = 0.0_SP
            endif
                UARDS = 0.0_SP
                VARDS = 0.0_SP
                IF(IOBCN > 0) THEN
                    UARD_OBCN(1:IOBCN)=0.0_SP
                END IF
            endif
        !end defined (TWO_D_MODEL)
            if  (BALANCE_2D==1)then
            ADVUA2    = 0.0_SP
            ADVVA2    = 0.0_SP
            ADFX2     = 0.0_SP
            ADFY2     = 0.0_SP
            DRX2D2    = 0.0_SP
            DRY2D2    = 0.0_SP
            CORX2     = 0.0_SP
            CORY2     = 0.0_SP
            PSTX2     = 0.0_SP
            PSTY2     = 0.0_SP
            ADX2D2    = 0.0_SP
            ADY2D2    = 0.0_SP
            WUSURBF2  = 0.0_SP
            WVSURBF2  = 0.0_SP 
            DUDT2     = 0.0_SP 
            DVDT2     = 0.0_SP
            DIVX2D2   = ZERO 
            DIVY2D2   = ZERO 
            DEDT2     = ZERO  
            endif

        ! New Open Boundary Condition ----4
            if (MEAN_FLOW==1)then
                CALL ZERO_OBC3
            endif   
        !==============================================================================!
        !  LOOP OVER EXTERNAL TIME STEPS                                               !
        !==============================================================================!
            DO IEXT=1,ISPLIT
                TIME  =(DTI*FLOAT(IINT-1)+DTE*FLOAT(IEXT))/86400.0_SP
                THOUR1=(DTI*FLOAT(IINT-1)+DTE*FLOAT(IEXT))/3600.0_SP
        !
        !----- USE RAMP VARIABLE TO EASE MODEL STARTUP---------------------------------!
        !
        !       TMP1 = FLOAT(IINT-ISTART)+FLOAT(IEXT)/FLOAT(ISPLIT)
                TMP1 = FLOAT(IINT-1)+FLOAT(IEXT)/FLOAT(ISPLIT)
                RAMP = 1.0_SP
                IF(IRAMP /= 0) RAMP = TANH(TMP1/FLOAT(IRAMP))
        !
        !------SURFACE BOUNDARY CONDITIONS FOR EXTERNAL MODEL--------------------------!
        !
                if  (GCN==1)then
                CALL BCOND_GCN(9)
                else
                CALL BCOND_GCY(9)
                endif
        !
        !------SAVE VALUES FROM CURRENT TIME STEP--------------------------------------!
        !
                ELRK1 = EL1
                ELRK  = EL
                UARK  = UA
                VARK  = VA

                if  (EQUI_TIDE==1)then
                ELRK_EQI = EL_EQI
                endif

        ! New Open Boundary Condition ----5
            if (MEAN_FLOW==1)then
                ELRKT  = ELT
                ELRKP  = ELP
                UARKNT = UANT   !change BKI
                VARKNT = VANT   !change BKI
                UARKN  = UAN    !change BKI
                VARKN  = VAN    !change BKI
            endif 
!
!------BEGIN MAIN LOOP OVER EXTERNAL MODEL 4 STAGE RUNGE-KUTTA INTEGRATION-----!
!
            DO K=1,4
                TIMERK = TIME + (ALPHA_RK(K)-1.)*DTE/86400.0_SP
! New Open Boundary Condition ----6
                if  (MEAN_FLOW==1)then
                    CALL BCOND_TIDE_2D
            !       CALL BCOND_NG_2D        ! change BKI
                    CALL FLUX_OBN2D(K)
                    CALL FLUX_OBC2D
                endif 
!FREE SURFACE AMPLITUDE UPDATE  --> ELF
                CALL EXTEL_EDGE(K)
        !JQI<
                if  (EQUI_TIDE==1)then
                    CALL ELEVATION_EQUI
	            ELF_EQI = ELRK_EQI +ALPHA_RK(K)*(ELF_EQI-ELRK_EQI) 
                endif
        !JQI>
! New Open Boundary Condition ----7
                if  (MEAN_FLOW==1)then
                    IF (ntidenode > 0) THEN
                        DO I = 1, ntidenode     ! need to calculate NEXT_OBC column of ELPF
                            J = I_TIDENODE_N(I)
                            ELPF(I) = ELF(J) - ELTF(I)
                        END DO
                    END IF

                    CALL EXTELPF_EDGE(K)

                    IF (IOBCN > 0) THEN
                        DO I = 1, IOBCN
                            J = I_OBC_N(I)
                            J1= I_OBC_NODE(J)
                            ELF(J) = ELTF(J1) + ELPF(J1)
                        END DO
                    END IF
                else                                             !  加潮位边界
                    if  (GCN==1)then
                        CALL BCOND_GCN(1)
                    else
                        CALL BCOND_GCY(1)
                    endif
                    DO I=1,IBCN(1)
                    JN = OBC_LST(1,I)
                    J=I_OBC_N(JN)
                    ELF(J)=ELRK(J)+ALPHA_RK(K)*(ELF(J)-ELRK(J))
                    END DO
                endif
            CALL N2E2D(ELF,ELF1)
            IF(WET_DRY_ON)CALL WET_JUDGE
        !JQI<added on 07/14/04
            CALL FLUX_OBN(K)
        !>JQI
        !CALCULATE ADVECTIVE, DIFFUSIVE, AND BAROCLINIC MODES --> UAF ,VAF
            if  (GCN==1)then
                CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
            else
                CALL ADVAVE_EDGE_GCY(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
            endif
            CALL EXTUV_EDGE(K)
            if  (GCN==1)then
                CALL BCOND_GCN(2)
            else
                CALL BCOND_GCY(2)
            endif
        !UPDATE WATER SURFACE ELEVATION
        !JQI<added on 07/14/04
            CALL ASSIGN_ELM1_TO_ELM2
        !>JQI 07/14/04
            EL  = ELF
            EL1 = ELF1
        !JQI<	 
            if(EQUI_TIDE==1)then
                EL_EQI = ELF_EQI
            endif
        !JQI>                     
        !!INTERPOLATE DEPTH FROM NODE-BASED TO ELEMENT-BASED VALUES
            CALL N2E2D(EL,EL1)
    !UPDATE DEPTH AND VERTICALLY AVERAGED VELOCITY FIELD
            D   = H + EL
            D1  = H1 + EL1
            UA  = UAF
            VA  = VAF
            DTFA = D
    ! New Open Boundary Condition ----8
            if  (MEAN_FLOW==1) then
                ELT = ELTF
                UAT = UATF
                VAT = VATF
    ! The part below is equivalent to do both NODE MATCH and EXCHANGE for ELPF

            ELP = ELPF
                IF(K == 4 .and. nmfcell > 0)THEN
                DO I = 1, nmfcell
                    J = I_MFCELL_N(I)
                    OBC2D_X_TIDE(I) = OBC2D_X_TIDE(I) + UANT(I) * D1(J)
                    OBC2D_Y_TIDE(I) = OBC2D_Y_TIDE(I) + VANT(I) * D1(J)
                ENDDO
                END IF
                CALL BCOND_BKI_2D(K)            ! change BKI
            endif


    !!ENSURE ALL CELLS ARE WET IN NO FLOOD/DRY CASE  
    !CALL DEPTH_CHECK
    !EXCHANGE ELEMENT-BASED VALUES ACROSS THE INTERFACE
            if (TWO_D_MODEL==0)then
    !SAVE VALUES FOR 3D MOMENTUM CORRECTION AND UPDATE
                IF(K == 3)THEN
                UARD = UARD + UA*D1
                VARD = VARD + VA*D1
                EGF  = EGF  + EL/ISPLIT
        !JQI<	 
                if (EQUI_TIDE==1)then
                EGF_EQI = EGF_EQI + EL_EQI/ISPLIT
                endif
     
        !JQI>                     
                UARDS = UARDS + UAS*D1
                VARDS = VARDS + VAS*D1
                END IF
                !CALCULATE VALUES USED FOR SALINITY/TEMP BOUNDARY CONDITIONS
                IF(K == 4.AND.IOBCN > 0) THEN
                DO I=1,IOBCN
                    J=I_OBC_N(I)
                    TMP=-(ELF(J)-ELRK(J))*ART1(J)/DTE-XFLUX_OBCN(I)
                    UARD_OBCN(I)=UARD_OBCN(I)+TMP/FLOAT(ISPLIT)
                END DO
                END IF
            endif !end defined (TWO_D_MODEL)

                    !UPDATE WET/DRY FACTORS
            IF(WET_DRY_ON)CALL WD_UPDATE(1)

            END DO     !! END RUNGE-KUTTA LOOP
        !  ************************   judge  umax  and vmax********************************** 
            uvmax=0.
	        do i=1,n
                if(sqrt(ua(i)**2+va(i)**2)>uvmax) uvmax=sqrt(ua(i)**2+va(i)**2)
	                if(uvmax>=10) then
                        write(*,*)' the velocity is too large , stop and check your flow field'
	                    write(*,*)  'vmax=', uvmax, 'imax=', i
	                    call tecplot
	                    stop
	            endif
            enddo
    END DO     !! EXTERNAL MODE LOOP

    !==============================================================================!
    !  END LOOP OVER EXTERN STEPS                                                  !
    !==============================================================================!
    if  (TWO_D_MODEL==0)then

    !==============================================================================!
    !    ADJUST INTERNAL VELOCITY FIELD TO CORRESPOND TO EXTERNAL                  !
    !==============================================================================!
    CALL ADJUST2D3D(2)
    ! New Open Boundary Condition ----9
    if  (MEAN_FLOW==1)then
    !     CALL BCOND_NG_3D       ! change BWI
        CALL BCOND_BKI_3D(2)      
        CALL FLUX_OBC3D
    endif

    !==============================================================================!
    !     CALCULATE INTERNAL VELOCITY FLUXES                                       |
    !==============================================================================!
                                !                                                    !
    CALL VERTVL_EDGE     ! Calculate/Update Sigma Vertical Velocity (Omega)   !
    IF(WET_DRY_ON) CALL WD_UPDATE(2)

    if  (GCN==1)then
        CALL ADV_UV_EDGE_GCN ! Horizontal Advect/Diff + Vertical Advection        !
    else
        CALL ADV_UV_EDGE_GCY ! Horizontal Advect/Diff + Vertical Advection        !
    endif
    CALL VDIF_UV         ! Implicit Integration of Vertical Diffusion of U/V  !

    DO I=1,N
        IF(H1(I) <= DJUST ) THEN
            DO K=1,KBM1
                UF(I,K)=UA(I)
                VF(I,K)=VA(I)
            END DO
        END IF
    END DO

    if  (GCN==1)then
        CALL BCOND_GCN(3)    ! Boundary Condition on U/V At River Input           !
    else
        CALL BCOND_GCY(3)    ! Boundary Condition on U/V At River Input           !
    endif
    CALL WREAL           ! Calculate True Vertical Velocity (W)               !
                                !                                                    !
			  
    CALL VISCOF_H        ! Calculate horizontal diffusion coefficient for     !
     			        ! the scalar                                         !

    !==============================================================================!
    !    TURBULENCE MODEL SECTION                                                  |
    !==============================================================================!


    !=================General Ocean Turbulence Model==========================!
    !===================Original FVCOM MY-2.5/Galperin 1988 Model=============!
    CALL ADV_Q(Q2,Q2F)       !!Advection of Q2 
    CALL ADV_Q(Q2L,Q2LF) 
    IF(TS_FCT) CALL FCT_Q2                             !Conservation Correction   !
    IF(TS_FCT) CALL FCT_Q2L                             !Conservation Correction   !
    CALL VDIF_Q                  !! Solve Q2,Q2*L eqns for KH/KM/KQ 
    DO I=1,M
        DO K=1,KBM1
            Q2(I,K) =Q2F(I,K)
            Q2L(I,K)=Q2LF(I,K)
        END DO
    END DO

    CALL N2E3D(KM,KM1)
    !==============================================================================!
    !    SEDIMENT MODEL SECTION                                                    |
    !==============================================================================!
    if (SEDIMENT==1)then
        ALLOCATE(FTEMP(0:MT),FTEMP2(0:NT))
        FTEMP2 = SQRT(WUBOT**2 + WVBOT**2)
        CALL E2N2D(FTEMP2,FTEMP)
        IF(SEDIMENT_ON)CALL ADVANCE_SED(DTI,THOUR*3600,FTEMP)
        DEALLOCATE(FTEMP,FTEMP2)
    !     call tecplot
    !	 stop
    endif
    !==============================================================================!
    !    UPDATE TEMPERATURE IN NON-BAROTROPIC CASE                                 !
    !==============================================================================!
    IF(TEMP_ON)THEN                                !                          !
    CALL ADV_T                                     !Advection                 !
    !#                                                   if !defined (DOUBLE_PRECISION)
    IF(TS_FCT) CALL FCT_T                             !Conservation Correction   !
    !#                                                   endif
    !qxu{
            IF(CASENAME == 'gom')THEN
            CALL VDIF_TS_GOM(1,TF1)
            ELSE  
            CALL VDIF_TS(1,TF1)                            !Vertical Diffusion        !
            END IF
    !qxu}     
            CALL BCOND_TS(1)                               !Boundary Conditions       !
    !#    if defined (DATA_ASSIM)
    !qxu{     CALL TEMP_NUDGING                              !Nudge Temperature         !
    !     IF(TS_ASSIM)THEN
    !      IF(ASSIM_FLAG==0 .AND. .NOT. SST_ASSIM)CALL TEMP_NUDGING
    !      IF(ASSIM_FLAG==1 .AND. SST_ASSIM)CALL TEMP_NUDGING
    !     END IF
    !qxu}     
    !#    endif
            T1 = TF1                                       !Update to new time level  !
            CALL N2E3D(T1,T)                               !Shift to Elements         !
            END IF                                         !                          !
    !==============================================================================!
    !    UPDATE SALINITY IN NON-BAROTROPIC CASE                                    !
    !==============================================================================!
            IF(SALINITY_ON)THEN                            !                          !
            CALL ADV_S                                     !Advection                 !

    !#                                                   if !defined (DOUBLE_PRECISION)
            IF(TS_FCT) CALL FCT_S                             !Conservation Correction   !
    !#                                                   endif
            CALL VDIF_TS(2,SF1)                            !Vertical Diffusion        !
            CALL BCOND_TS(2)                               !Boundary Conditions       !
    !#    if defined (DATA_ASSIM)
    !qxu{     CALL SALT_NUDGING                              !Nudge Salinity            !
    !     IF(TS_ASSIM)THEN
    !      IF(ASSIM_FLAG==0 .AND. .NOT. SST_ASSIM)CALL SALT_NUDGING
    !      IF(ASSIM_FLAG==1 .AND. SST_ASSIM)CALL SALT_NUDGING
    !     END IF
    !qxu}      
    !#    endif
            S1 = SF1                                       !Update to new time level  !
            CALL N2E3D(S1,S)                               !Shift to Elements         !
            END IF                                         !                          !

    !==============================================================================!
        if  (DYE_RELEASE==1)then

    !==============================================================================!
    !    UPDATE DYE IN NON-BAROTROPIC CASE                                         !
    !==============================================================================!
    !     IF(DYE_ON.AND.IINT.GE.IINT_SPE_DYE_B) THEN     !                          !                          !
            IF(DYE_ON) THEN     !                          !                          !
            CALL ADV_DYE                                   !Advection                 !





            CALL VDIF_DYE(DYEF)                            !Vertical Diffusion        !

    !check
    !!     DYE = DYEF                                     !Update to new time level  !
    !!     IF(IINT.GE.IINT_SPE_DYE_B) CALL ARCHIVE
    !check    
            CALL BCOND_DYE                                 !Boundary Conditions       !
            DYE = DYEF                                     !Update to new time level  !

    !!     IF(MSR) WRITE(IPT,*) 'CALL Dye_on--iint=',iint,IINT_SPE_DYE_B
            END IF                                         !                          !
    !==================================================================================!
        endif
    !QXU}

        if  (MPDATA==0)then
            IF(POINT_ST_TYPE == 'calculated')THEN
    !    ADJUST TEMPERATURE AND SALINITY AT RIVER MOUTHS
            CALL ADJUST_TS
            END IF  
        endif   

        if  (WATER_QUALITY==1)then
    !==============================================================================!
    !    CALCULATE WATER QUALITY VARIABLES CONCENTRATIONS                          |
    !==============================================================================!
                                                        !                          !
            IF(WQM_ON)THEN                                 !Water Quality Active?     !
            CALL ADV_WQM                                   !Advection                 !
            CALL VDIF_WQM(WQM_F)                           !Vertical Diffusion        !
            CALL EXCHANGE_WQM                              !Interprocessor Exchange   !
            CALL BCOND_WQM                                 !Boundary Conditions       !
            WQM(1:M,1:KBM2,1:NB) = WQM_F(1:M,1:KBM2,1:NB)  !Update                    !
            END IF                                         !                          !
    !==============================================================================!
        endif


    !#    if defined (DATA_ASSIM)
    !==============================================================================!
    !     ASSIMILATE SEA SURFACE TEMPERATURE DATA                                  |
    !==============================================================================!
    !     IF(SST_ASSIM .AND. ASSIM_FLAG == 1) CALL SST_NUDGING !NUDGE SST           !
    !     IF(SST_ASSIM .AND. ASSIM_FLAG == 0) CALL SST_INT     !STORE HOURLY SST    !
    !     IF(SST_ASSIM .AND. ASSIM_FLAG == 1) CALL N2E3D(T1,T) !RECALCULATE T       !
    !qxu{
    !
    !----Recalculate Element Based Temperatures Based on Nudged Temp Field T1------!
    !
    !     IF(SST_ASSIM .AND. ASSIM_FLAG == 1)CALL N2E3D(T1,T)
    !qxu}
    !==============================================================================!
    !#    endif


    !==============================================================================!
    !     UPDATE THE DENSITY IN NON-BAROTROPIC CASE                                |
    !==============================================================================!
            IF(.NOT.BAROTROPIC)THEN
            IF(CTRL_DEN == 'pdensity'   ) CALL DENS
            IF(CTRL_DEN == 'sigma-t'    ) CALL DENS2
            IF(CTRL_DEN == 'sigma-t_stp') CALL DENS3
            END IF
    !==============================================================================!
    !     MIMIC CONVECTIVE OVERTURNING TO STABILIZE VERTICAL DENSITY PROFILE       |
    !==============================================================================!

            IF(VERT_STAB)THEN
            CALL CONV_OVER
            IF(.NOT.BAROTROPIC)THEN
                IF(CTRL_DEN == 'pdensity'   ) CALL DENS
                IF(CTRL_DEN == 'sigma-t'    ) CALL DENS2
                IF(CTRL_DEN == 'sigma-t_stp') CALL DENS3
            END IF
            END IF  
    !#  endif     
        end if !defined (TWO_D_MODEL)

    !#    if defined (DATA_ASSIM)
    !==============================================================================!
    !     DATA ASSIMILATION FOR CURRENT FIELD                                      |
    !==============================================================================!
    !     IF(CURRENT_ASSIM)THEN
    !       IF(ASSIM_FLAG==0 .AND. .NOT. SST_ASSIM)CALL CURRENT_NUDGING
    !       IF(ASSIM_FLAG==1 .AND. SST_ASSIM)CALL CURRENT_NUDGING
    !# if defined (GCN)
    !       CALL BCOND_GCN(5)
    !# else
    !       CALL BCOND_GCY(5)
    !# endif
    !     END IF
    !==============================================================================!
    !#    endif

    !==============================================================================!
    !    LAGRANGIAN PARTICLE TRACKING                                              |
    !==============================================================================!
    !     IF(LAG_ON) CALL LAG_UPDATE
   
        if(TWO_D_MODEL==0)then
    !==============================================================================!
    !     UPDATE VELOCITY FIELD (NEEDED TO WAIT FOR SALINITY/TEMP/TURB/TRACER)     |
    !==============================================================================!
            U = UF
            V = VF
    !==============================================================================!
    !    PERFORM DATA EXCHANGE FOR ELEMENT BASED INFORMATION AT PROC BNDRIES       |
    !==============================================================================!


    !==============================================================================!
    !     PERFORM DATA EXCHANGE FOR WATER QUALITY VARIABLES                        |
    !==============================================================================!
        if (WATER_QUALITY==1)then
            CALL EXCHANGE_WQM
        endif
        endif  !defined (TWO_D_MODEL)    

    !
    !----SHIFT SEA SURFACE ELEVATION AND DEPTH TO CURRENT TIME LEVEL---------------!
    !
            ET  = EL  
            DT  = D 
            ET1 = EL1
            DT1 = D1
            IF(WET_DRY_ON) CALL WD_UPDATE(3)

    ! New Open Boundary Condition ----10
        if  (MEAN_FLOW==1)then
            ELTDT = ELT
        endif

    !==============================================================================!
    !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
    !==============================================================================!

            IF(MSR)CALL REPORT_TIME(IINT,ISTART,IEND,TIME*86400,IPT) 
            IF(MOD(IINT,IREPORT)==0) CALL REPORT("FLOW FIELD STATS")
            CALL DUMP_PROBE_DATA 

    !#    if !defined (AIX) && !defined(ABSOFT)
    !     CALL FLUSH(6)
    !#    endif
    !
    !-------------UPDATE THE VISUALIZATION SERVER----------------------------------!
    !
    !#    if defined (PV3)
    !     CALL PV_UPDATE(THOUR)
    !#    endif

    !#    if !defined (DATA_ASSIM)
    !     CALL ARCHIVE
    !#    endif

    ! New Open Boundary Condition ----11
    !     CALL PRINT_VALS

    !
    !---------------WRITE OUTPUT FILES---------------------------------------------!
    !
    !#    if defined (DATA_ASSIM)
    !     IF((SST_ASSIM .AND. ASSIM_FLAG==1) .OR. PURE_SIM) CALL  ARCHIVE
    !#    endif
  
            CALL SHUTDOWN_CHECK

            if(mod(IINT,360)==0)then
            call tecplot
            endif


    END DO !!MAIN LOOP

    !  close(112)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    END MAIN LOOP OVER PHYSICAL TIME STEPS                                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#  if defined (DATA_ASSIM)
!   IF(SST_ASSIM .AND. ASSIM_FLAG==1) CALL ARCRST
!#  endif

!#  if defined (DATA_ASSIM)
!   END DO
!   END DO
!#  endif
!==============================================================================!
!  END MAIN LOOP OVER DATA ASSIMILATION INTERVALS AND SWEEP #                  !
!==============================================================================!
!==============================================================================!
!  CLOSE UP COMPUTATION                                                        !
!==============================================================================!

    IINT = IEND
    CALL CLOSEFILES
    CALL ARCRST
    IF(MSR)THEN
        WRITE(IPT,*)'DUMPING RESTART'
    END IF
    IF(WET_DRY_ON) CALL WD_DUMP

    CALL REPORT('FINAL VALUES INFORMATION')

    IF(MSR)THEN
        WRITE(IPT,*) ; WRITE(IPT,*)'Computation completed, congratulations!'
        CALL GET_CLOCK 
        CALL GETTIME(TSTRING,INT(TCURRENT-TINIT))
    END IF
    CALL PSTOP
    pause
!
!----------------------FORMAT STATEMENTS---------------------------------------|
!
7000 FORMAT(1X,A28,A13)  
7001 FORMAT(1X,A28,I8)  

    END PROGRAM US_FVCOM




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!==============================================================================!

!==============================================================================!
!   IMPORT CASENAME FROM COMMAND LINE                                          !
!==============================================================================!

 