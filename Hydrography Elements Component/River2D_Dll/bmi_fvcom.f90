Module FVCOM
!这里的代码参考的us_fvcom，为了压缩代码行数，删除了部分注释，假如想要了解各个变量具体含义，请自行看us_fvcom里的说法
    USE ALL_VARS
    USE MOD_CLOCK
    USE MOD_SPHERICAL
    USE MOD_WD
    USE BCS
    USE PROBES  
    USE MOD_OBCS   
    USE MOD_MEANFLOW
    USE MOD_OBCS2
    USE MOD_OBCS3
    USE MOD_WQM
    USE MOD_NORTHPOLE
    USE MOD_EQUITIDE
    USE MOD_SED
    USE MOD_DYE
    USE MOD_TSOBC
    USE MOD_BALANCE_2D
    use wave
    use iso_c_binding, only: c_ptr, c_f_pointer, c_loc, c_null_char,c_int,c_char
    LOGICAL  :: FEXIST
    REAL(sp) :: TMP1,TMP,UTMP,VTMP,TTIME,  uvmax
    INTEGER  :: I,K,J,IERR,N1,J1,J2,JN
    REAL(sp), ALLOCATABLE :: FTEMP(:),FTEMP2(:)
    CHARACTER(LEN=13) :: TSTRING
    ! 王浩丞改：读文件
    contains
    
!初始化FVCOM
subroutine fvcom_init(cstring) bind(c,name='fvcom_init')
    !DEC$ ATTRIBUTES DLLEXPORT::fvcom_init
    IMPLICIT NONE
    type(c_ptr), dimension(3),target, intent(in) :: cstring
    character(len=80),dimension(3) ::data1
    character, pointer                            :: fstring(:)
    integer,dimension(3)                         :: slen
    do i = 1, 3
        call c_f_pointer(cstring(i), fstring,[80])
        slen(i) = 0
        do while(fstring(slen(i)+1) /= c_null_char)
            slen(i) = slen(i) + 1
        end do
        data1(i) = transfer(fstring(1:slen(i)), data1(i))
    end do
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
    INPDIR= data1(1)(1:slen(1))!"..\run\"
    OUTDIR= data1(2)(1:slen(2))!"..\output\"
    casename= data1(3)(1:slen(3))!'mdm'
    casename1= trim(INPDIR)//trim(casename)//"_run.dat" 
    INQUIRE(FILE=casename1,EXIST=FEXIST)
        IF(.NOT.FEXIST)THEN
            WRITE(IPT,*)'FILE ',TRIM(casename)//"_run.dat",' DOES NOT EXIST'
            WRITE(IPT,*)'STOPPING...'
            CALL PSTOP
        END IF
    spherical=0
    WET_DRY=1
    WATER_QUALITY=0
    DATA_ASSIM=0
    GCN=1
    GCY1=0
    GCY2=0
    EQUI_TIDE=0
    ATMO_TIDE=0
    NORTHPOLE=0
    MPDATA=0
    TWO_D_MODEL=1
    BALANCE_2D=0
    TS_OBC=0
    MEAN_FLOW=0
    DYE_RELEASE=0
    sediment=0
    iwave=0
    icyclone=0
    CALL DATA_RUN
    CALL TSOBC_TYPE
    CALL SET_WD_PARAM
    if (BALANCE_2D==1)then
        CALL SET_BALANCE_PARAM
        endif
    if (DYE_RELEASE==1)then
        CALL SET_DYE_PARAM
    endif   
    if  (TS_OBC==1)then
        CALL SET_TSOBC_PARAM
    endif   
    !  OPEN INPUT/OUTPUT FILES
    CALL IOFILES
    !  DETERMINE NUMBER OF ELEMENTS AND NODES IN THE MODEL
    CALL GETDIM
    !  READ PARAMETERS CONTROLLING WATER QUALITY MODEL 
    if (WATER_QUALITY==1)then
        CALL GET_WQMPAR
    endif
    !  MAP OPEN BOUNDARY CONDITION NODES TO LOCAL DOMAIN
    CALL BCMAP           !  开边界位置，类型
    !  INPUT AND SETUP BOUNDARY FORCING (HEAT/RIVERS/WIND/etc)
    CALL BCS_FORCE           ! 边界条件数据
    !  INPUT WATER QUALITY MODEL VARIABLES
    if (WATER_QUALITY==1)then
    IF(WQM_ON)CALL BCS_FORCE_WQM
    endif
    !  INPUT T/S OBC Series NUDGING VARIABLES
    IF(TSOBC_ON) CALL READ_TSOBC
    !  ALLOCATE FLOWFIELD VARIABLES
    CALL ALLOC_VARS
    !  ALLOCATE ELM1 AND ELM2 FOR ORLANSKI RADIATION OPEN BOUNDARY CONDITION
    CALL ALLOC_OBC_DATA
    !  ALLOCATE AND WET/DRY CONTROL ARRAYS
    CALL ALLOC_WD_DATA
    !  ALLOCATE WATER QUALITY MODEL VARIABLES
    if (WATER_QUALITY==1)then
    IF(WQM_ON)CALL ALLOC_WQM_VARS
    endif
        !  ALLOCATE SPHERICAL COORDINATE SYSTEM VARS
        if (SPHERICAL==1)then
            CALL ALLOC_SPHERE_VARS
        endif
        !  ALLOCATE MOMENTUM BALANCE CHEKING VARS
        if  (BALANCE_2D==1)then
            CALL ALLOC_BALANCE_VARS
        endif
        !  ALLOCATE DYE VARS
        if (DYE_RELEASE==1)then
            CALL ALLOC_VARS_DYE
        endif
        !  SHIFT GRID/CORIOLIS/BATHYMETRY TO LOCAL DOMAIN
        CALL PDOMDEC      ! 读网格
        !  ALLOCATE EQUILIBRIUM AND ATMOSPHERIC TIDE VARS
        if (EQUI_TIDE==1)then
            CALL ALLOCATE_EQUI
        endif  
        !  SET UP GRID METRICS (FLUX EDGES/CONTROL VOLUMES/ETC)
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
            call wavechazhi
            call wavepump
        endif
    !  SETUP OPEN BOUNDARY METRICS                                   
        CALL SETUP_OBC
        CALL SET_BNDRY               !Boundary Condition Metrics  inflow condtion
    !  INITIALIZE FLOWFIELD  --> [T,S,U,V,EL,D,Q2,Q2L]
        CALL STARTUP
    !  INITIALIZE GOTM
    !  INITIALIZE SEDIMENT MODEL
        if(sediment==1) then
            IF(SEDIMENT_ON)CALL SETUP_SED(RESTART_SED)
        endif   
    !  CALCULATE DEPTH HORIZONTAL DERIVATIVES
        CALL DEPTH_GRADIENT
        CALL REPORT('INITIAL VALUE INFORMATION')
    !  CALCULATE INTERNAL TIME STEP AND SET INTEGRATION LIMITS
        ISTART=IINT+1
        CALL START_CLOCK
    !  REPORT INTEGRATION INITIAL TIMES 
        IF(MSR)CALL REPORT_SIMTIME
    ! New Open Boundary Condition ----2
        if  (MEAN_FLOW==1)then
            CALL FIND_OBSIDE
            CALL ALLOC_OBC3_DATA
            CALL SETUP_OBC3
            CALL ALLOC_OBC2_DATA
        endif
        open(112, file=TRIM(OUTDIR)//"/"//trim(casename)//'_topo.dat')
        IINT=ISTART
    end subroutine fvcom_init
!更新FVCOM
subroutine update_fvcom() bind(c,name='update_fvcom')
       !DEC$ ATTRIBUTES DLLEXPORT::update_fvcom
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
            !IF(MSR)CALL REPORT_TIME(IINT,ISTART,IEND,TIME*86400,IPT) 
            !IF(MOD(IINT,IREPORT)==0) CALL REPORT("FLOW FIELD STATS")
            CALL DUMP_PROBE_DATA 
            CALL SHUTDOWN_CHECK

            if(mod(IINT,36)==0)then
            call tecplot
            endif
    IINT = IINT+1
end subroutine update_fvcom
!结束FVCOM
subroutine finalize_fvcom() bind(c,name='finalize_fvcom')
    !DEC$ ATTRIBUTES DLLEXPORT::finalize_fvcom
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
    !关闭两个监测点
    close(1121)
    close(1122)
    pause
    7000 FORMAT(1X,A28,A13)  
    7001 FORMAT(1X,A28,I8)  
end subroutine finalize_fvcom
!获取开始时间
subroutine Get_StartTime(START) bind(c,name='Get_StartTime')
    !DEC$ ATTRIBUTES DLLEXPORT::Get_StartTime
    real(4), intent(out) :: START
    START=DTI*FLOAT(ISTART)
end subroutine Get_StartTime
!获取结束时间
subroutine Get_EndTime(I_END) bind(c,name='Get_EndTime')
    !DEC$ ATTRIBUTES DLLEXPORT::Get_EndTime
    real(4), intent(out) :: I_END
    I_END = DTI*FLOAT(IEND)
end subroutine Get_EndTime
!获取当前时间
subroutine Get_current_time_fvcom(I_Int) bind(c,name='Get_current_time_fvcom')
    !DEC$ ATTRIBUTES DLLEXPORT::Get_current_time_fvcom
    I_Int=DTI*FLOAT(IINT)
end subroutine Get_current_time_fvcom
subroutine Set_Obc(shuiwei,I) bind(c,name='Set_Obc')
    !DEC$ ATTRIBUTES DLLEXPORT::Set_Obc
    real ,intent(in):: shuiwei
    integer , value:: I
    ELSBC(I,:)=shuiwei
end subroutine Set_Obc
!FVCOM抛出水位
!h 水位，i不规则格点索引
subroutine get_depth_fvcom(hhh,i) bind(c,name='get_depth_fvcom')
    !DEC$ ATTRIBUTES DLLEXPORT::get_depth_fvcom
    integer , value:: i
    hhh=D(i)
end subroutine get_depth_fvcom
!FVCOM抛出水位
!h 水位，i不规则格点索引
subroutine get_elevation_fvcom(hhh,i) bind(c,name='get_elevation_fvcom')
    !DEC$ ATTRIBUTES DLLEXPORT::get_elevation_fvcom
    integer , value:: i
    hhh=H(i)
end subroutine get_elevation_fvcom
subroutine get_water_depth_fvcom(hhh,i) bind(c,name='get_water_depth_fvcom')
    !DEC$ ATTRIBUTES DLLEXPORT::get_water_depth_fvcom
    integer , value:: i
    hhh=D(i)-H(i)
end subroutine get_water_depth_fvcom
!FVCOM接收流量
!Q 流量 i流量边界号
subroutine LISFLOOD_To_FVCOM(Q,i) bind(c,name='LISFLOOD_To_FVCOM')
    !DEC$ ATTRIBUTES DLLEXPORT::LISFLOOD_To_FVCOM
    integer ,value:: i
    real ,intent(in):: Q
    DO I = 1, QBC_TM%NTIMES
        DQDIS(1:NUMQBC_GL,I)=Q
    END DO
end subroutine LISFLOOD_To_FVCOM
    end Module FVCOM
    
    