!==============================================================================!
!   Open Input Files for Model Parameters and Output Files for Results         !
!==============================================================================!


   SUBROUTINE IOFILES
   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_BALANCE_2D
   USE MOD_MEANFLOW
   IMPLICIT NONE
   LOGICAL:: CHECK 
   CHARACTER(LEN=80)  :: TEMP,ISTR,OSTR
   CHARACTER(LEN=100) :: MKOUTDIR,MKSMSDIR,MKMDMDIR,MKRESDIR,MKTSODIR,MKNCDDIR
   INTEGER IERR,ISTAT

!==============================================================================!
!                Definitions of input files                                    !
!                                                                              !
! inrun: casename_run.dat: input data and parameter file controlling           !
!                          the model run                                       !
! indep: casename_dep.dat: input water depth at the rest. For estuary,         !
!                          water depth refers to water depth at the            !
!                          lowest water level. The adjustment (adjust)         !
!                          must be specified in "casename_run.dat"             !
! ingrd: casename_grd.dat: input triangular mesh files. This file is           !
!                          generated from the SMS mesh generation              !
! inobc: casename_obn.dat: input data for open boundary triangles              !
!                          nodes                                               !            
! incwh: casename_mc.dat: boundary input meteorlogical forcing:                !
!                         wind velocity,heat flux, and                         !
!                         precipitation/evaporation. This file only            !
!                         work for a uniform meteorological forcing            !
!                         case.                                                !
! inriv: casename_riv.dat: river discharge input data that include             !
!                          number of rivers, transport, etc.                   !
! inits: casename_its.dat: initial temperature and salinity fields             !
! inoel: casename_el.dat: tidal amplitudes and phases at the open              !
!                         boundary (with no julian time)                       !
! injul: casename_elj.dat: tidal amplitudes and phases at the open             !
!                          boundary generated using Foreman s program          !
!                          with a julian time.                                 !
! inhfx: casename_hfx.dat: real-time field of heat flux input. In our          !
!                          current experiments, it is an ouput file            !
!                          from MM5 results.                                   !
! inwnd: casename_wnd.dat: real-time field of wind velocity or wind            !
!                          stress used for external mode. In our               !
!                          current experiments, it is an output file           !
!                          from MM5 results.                                   !                
! inelf: casename_el_ini.dat: the initial field of surface elevation           !
!                            at each triangular node calculated from           !
!                            the Foreman s harmonic tidal analysis             !
!                            program. This file needs to be generated          !
!                            based on calibrated amplitudes and phases         !
!                            of major tidal constituents from the              !
!                            model tidal simulation.                           !
! inuvf: casename_uv_ini.dat: the initial field of currents at each            !
!                            triangular centriod calculated from the           !
!                            Foreman s harmonic tidal analysis program         !
!                            . This file needs to be generated based           !
!                            on calibrated tidal ellipses of major             !
!                            tidal constituents from the model tidal           !
!                            simulation.                                       !
! inlag: casename_lag_ini.dat: the initial positions of particles for          !
!                             the lagrangian tracking                          !
! inspo: casename_spg.dat:   the parameters for sponge layers at the           !
!                            open boundary                                     !
! incor: casename_cor.dat: the latitudes of triangular nodes which are         !
!                          used to calculate the Coriolis parameter            !
! inbfw: casename_bfw.dat: the inFORMATion about bottom freshw water           !
!                          input                                               !
! injmp: casename_jmpobc.dat number of nodes and list of nodes for frictional  !
!                          geostrophic inflow correction                       !
! inmf : casename_meanflow.dat: open bndy volume transport input file that     !
!                         includes number of open boundary cell, transport,etc !
! intcell : casename_tide_cell.dat: the cell index of open bndy tidal cells    !
!                                   used in MEAN_FLOW bndy condition           !
! intuv   : casename_tide_uv.dat: the time series of u and v velocity at       !
!                                   tide_cell points                           !
! intnode : casename_tide_node.dat: the node index of open bndy tidal nodes    !
!                                   used in MEAN_FLOW bndy condition           !
! intelel : casename_tide_el.dat: the time series of water elevation at        !
!                                   tide_node points                           !
!==============================================================================!
         
!==============================================================================!
!                  Definitions of output files                                 !
!                                                                              !
! ioprt: casename_prt.dat: the file printing all input parametes, data         !
!                          and forcings                                        !
! ioplt: casename_plt.dat: the archive file including averaged field           !
!                          of currents, temperature, salinity,etc for          !
!                          a specified time interval (such as atidal           !
!                          cycle)                                              !
! iotsr: casename_tsr.dat: the archive file including time series of           !
!                           selected variables at selected locations           !            
! iomob: casename_mob.dat: the archive file including time series of           !
!                           momentum balance variables at selected locations   !            
! iopar: casename_lag_out.dat: the output trajectories of particles            !
!                              during tracking periods.                        !
! iopuv: casename_lag_ouv.dat: the output velocity of particles                !
!                              during tracking periods.                        !
! iosmsd: casename_dep.xy: the output file for water depth used for            !
!                          sms graphics software.                              !
! iosmsv: casename_uvi_uva.xy: the output files for internal and               !
!                              external velocities used for sms                !
!                              graphics software.                              !
! iosmst: casename_elts.xy:    the output fiels for elevation,                 !
!                              temperature and salinity used for sms           !
!                              graphics software.                              !
!                          iosmsv and isomst output files are updated          !
!                          at a specified time interval, so that they          !
!                          are opened in subroutine called  out_sms_one        !
! item90-94: the temporary files used to storage the data, which are           !
!            cleaned up at the end of the model run                            !
!                                                                              !
!  there are also some output files in subroutines of "out_binary" and         !
!  "out_binary_residual" for standard model binary output. All these           !
!  files could be directly used for graphics system developed by the           !
!  ocean ecosystem modeling group at SMAST/UMASSD.                             !
!==============================================================================!

   IF(MSR)WRITE(IPT,*)'!                                                                !'
   IF(MSR)WRITE(IPT,*)'!                  OPENING FILES                                 !'
   IF(MSR)WRITE(IPT,*)'!                                                                !'

   INDEP =11
   INGRD =12
   INOBC =13
   INCWH =14
   INRIV =15
   INITS =16
   INOEL =17
   INJUL =18
   INHFX =19
   INWND =20
   INELF =22
   INUVF =23
   INSPO =25
   INCOR =26
   INBFW =27

   IOPRT =41
   IOPLT =42
   IOTSR =43
!QXU{   
 if  (BALANCE_2D==1)  then
   IOMOB =44
 endif   
!QXU}

 if (MEAN_FLOW==1)  then
   INMF    =45
   INTCELL =46
   INTNODE =47
   INTELEL =48
   INTUV   =49
 endif   

   IOSMSD=51
   IOSMSV=52
   IOSMST=53

   IGOTM =59

   INRES=54
   INJMP=55
   IREST=60


!
!-----------------CHECK FOR EXISTENCE/CREATE DIRECTORIES-----------------------!
!



   IF(MSR)WRITE(IPT,*)
!
!---------CHECK EXISTENCE OF STANDARD FILES-AND OPEN---------------------------!
!
   ISTR = TRIM(INPDIR)//"/"//trim(casename)
   OSTR = TRIM(OUTDIR)//"/out"//trim(casename)
   CALL FOPEN(IOPRT, TRIM(OSTR)//'_prt.dat',"ofr")
   CALL FOPEN(INDEP, TRIM(ISTR)//'_dep.dat',"cfr")
   CALL FOPEN(INGRD, TRIM(ISTR)//'_grd.dat',"cfr")
   CALL FOPEN(INOBC, TRIM(ISTR)//'_obc.dat',"cfr")
   CALL FOPEN(INRIV, TRIM(ISTR)//'_riv.dat',"cfr")
   CALL FOPEN(INSPO, TRIM(ISTR)//'_spg.dat',"cfr")
  if (SPHERICAL==0) CALL FOPEN(INCOR, TRIM(ISTR)//'_cor.dat',"cfr")

   CALL FOPEN(INBFW, TRIM(ISTR)//'_bfw.dat',"cfr")


!
!-----------------INITIAL TEMPERATURE AND SALINITY-----------------------------!
!
   IF(RESTART /= 'hot_start') CALL FOPEN(INITS, TRIM(ISTR)//'_its.dat',"cfr")
!
!-----------------OPEN METEOROLOGICAL FORCING FILES----------------------------!
!
   IF(M_TYPE == 'uniform')THEN
     CALL FOPEN(INCWH, TRIM(ISTR)//'_mc.dat',"cfr")
   ELSE 
     CALL FOPEN(INHFX, TRIM(ISTR)//'_hfx.dat' ,"cur")
     CALL FOPEN(INWND, TRIM(ISTR)//'_wnd.dat' ,"cur")
   END IF

!
!-----------------OPEN TIDAL FORCING AND INITIAL VELOCITY FIELD FILES----------!
!
   IF(S_TYPE=='non-julian')THEN
     CALL FOPEN(INOEL,TRIM(ISTR)//'_el_obc.dat',"cfr")
   ELSE 
     CALL FOPEN(INJUL,TRIM(ISTR)//'_sea.dat',"cfr")
     IF(RESTART=='cold_start') THEN
  !     CALL FOPEN(INELF,TRIM(ISTR)//'_el_ini.dat' ,"cfr")
  !     CALL FOPEN(INUVF,TRIM(ISTR)//'_uv_ini.dat' ,"cfr")
     END IF
   END IF


!
!-----------------------FILES FOR GEOSTROPHIC CORRECTION AT INFLOW-------------!
!
   JMPOBC = .FALSE.
   INQUIRE(FILE=TRIM(ISTR)//'_jmpobc.dat',EXIST=CHECK) 
   IF(CHECK)THEN
     JMPOBC = .TRUE.
     CALL FOPEN(INJMP,TRIM(ISTR)//'_jmpobc.dat',"cfr")
   END IF
!
!-----------------------FILES FOR RESTART--------------------------------------!
!
   IF(RESTART == 'hot_start') CALL FOPEN(INRES,TRIM(ISTR)//'_restart.dat',"cur")
   IF(RESTART == 'hot_cold_s') CALL FOPEN(INRES,TRIM(ISTR)//'_restart.dat',"cur")
!
!-----------------------FILES FOR ARCHIVING------------------------------------!
!
   CALL FOPEN(IOPLT,TRIM(OSTR)//'_plt.dat' ,"our")
   CALL FOPEN(IOTSR,TRIM(OSTR)//'_tsr.dat' ,"ofr")

!
!-----------------------DEPTH OUTPUT FOR SMS PLOT------------------------------!
!
!   CALL FOPEN(IOSMSD,TRIM(OUTDIR)//"/sms/"//trim(casename)//"_dep.xy" ,"ofr")
!QXU{
  if (BALANCE_2D==1)CALL FOPEN(IOMOB,TRIM(OSTR)//'_mob.dat' ,"ofr")
   
!QXU}

  if (MEAN_FLOW==1)then
   CALL FOPEN(INMF,   TRIM(ISTR)//'_meanflow.dat'  ,"cfr")
   CALL FOPEN(INTCELL,TRIM(ISTR)//'_tide_cell.dat' ,"cfr")
   CALL FOPEN(INTNODE,TRIM(ISTR)//'_tide_node.dat' ,"cfr")
   CALL FOPEN(INTELEL,TRIM(ISTR)//'_tide_el.dat'   ,"cfr")
   CALL FOPEN(INTUV,  TRIM(ISTR)//'_tide_uv.dat'   ,"cfr")
  endif 


   RETURN
   END SUBROUTINE IOFILES
!==============================================================================!




