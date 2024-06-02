!==============================================================================|
!     SET BOUNDARY CONDITIONS FOR ALMOST ALL VARIABLES                         |
!                                                                              |
!         idx: identifies which variables are considered                       |
!              1=tidal forcing                                                 |
!              2=solid bcs for external mode uaf and vaf                       |
!              3=solid bcs for internal mode uf and vf                         |
!              4=open bcs for s and t                                          |
!              5=solid bcs for internal mode u and v                           |
!              6=unused                                                        |
!              7=unused                                                        |
!              8=the surface forcings for internal mode                        |
!              9=the surface forcings for external mode                        |
!                                                                              |
!==============================================================================|

   SUBROUTINE BCOND_GCY(IDX)

!==============================================================================|
   USE ALL_VARS
   USE BCS

!JQI<modified for new obcs on 07/14/04
   USE MOD_OBCS
!>JQI 07/14/04

   USE MOD_WQM

   USE MOD_WD
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: IDX
   REAL(SP) :: ZREF(KBM1),ZREFJ(KBM1),TSIGMA(KBM1),SSIGMA(KBM1)
   REAL(SP) :: TTMP(KBM1),STMP(KBM1),TREF(KBM1),SREF(KBM1)
   REAL(SP) :: PHY_Z(KBM1),PHY_Z1(KBM1)
   REAL(SP) :: TT1(KBM1),TT2(KBM1),SS1(KBM1),SS2(KBM1)
   REAL(SP) :: TIME1,FACT,UFACT,FORCE,QPREC,QEVAP,UI,VI,UNTMP,VNTMP,TX,TY,HFLUX
   REAL(SP) :: DTXTMP,DTYTMP,QPREC2,QEVAP2,SPRO,WDS,CD,SPCP,ROSEA
   REAL(SP) :: PHAI_IJ,ALPHA1,DHFLUXTMP,DHSHORTTMP,HSHORT,TIMERK1
   REAL(SP) :: ETAXFER,ANGALONG,ANGWIND,WNDALONG,TXJMP2,TYJMP2,RHOILST,RHOINT,CUMEL
   REAL(SP) :: DXBC,ETATAN,CC,DELTAN,CP

   INTEGER  I,J,K,I1,I2,J1,J2,II,L1,L2,IERR
   INTEGER  NNOW,NLAST,NBCJMP,IGL
   INTEGER  N1

   SELECT CASE(IDX)


!==============================================================================|
   CASE(1) !Surface Elevation Boundary Conditions (Tidal Forcing)              !
!==============================================================================|

!JQI<modified on 07/14/04
!   TIME1=TIME*86400.0_SP
!   TIMERK1 = TIMERK*86400
!
!-Julian: Set Elevation Based on Linear Interpolation Between Two Data Times---|
!
!   IF(S_TYPE == 'julian') THEN
!!     CALL BRACKET(ELO_TM,TIMERK1,L1,L2,FACT,UFACT,IERR)
!     IF(IOBCN > 0)CALL BRACKET(ELO_TM,TIME1,L1,L2,FACT,UFACT,IERR)   !!GWC REMOVE
!     DO I = 1, IOBCN
!       IF(TYPE_OBC(I) == 1)CALL BCOND_ASL_JULIANDAY(I,L1,L2,FACT,UFACT)
!       IF(TYPE_OBC(I) == 2)CALL BCOND_ASL_CLP(I)
!       IF(TYPE_OBC(I) == 3)CALL BCOND_GWI(I)
!       IF(TYPE_OBC(I) == 4)CALL BCOND_GWI(I)
!       IF(TYPE_OBC(I) == 5)CALL BCOND_ORE(I)
!       IF(TYPE_OBC(I) == 6)CALL BCOND_ORE(I)
!    DO J = 1, ISBCN
!       II=I_SBC_N(J)
!       ELF(II) = UFACT*ELSBC(J,L1) + FACT*ELSBC(J,L2)
!     END DO
!   END IF


!
!-Non-Julian: Set Elevation Based on Input Amplitude and Phase of Tidal Comps-|
!

!   IF(S_TYPE == 'non-julian') THEN
!     DO I = 1, IOBCN
!       IF(TYPE_OBC(I) == 1)CALL BCOND_ASL(I,TIME1)
!       IF(TYPE_OBC(I) == 2)CALL BCOND_ASL_CLP(I)
!       IF(TYPE_OBC(I) == 3)CALL BCOND_GWI(I)
!       IF(TYPE_OBC(I) == 4)CALL BCOND_GWI(I)
!       IF(TYPE_OBC(I) == 5)CALL BCOND_ORE(I)
!       IF(TYPE_OBC(I) == 6)CALL BCOND_ORE(I)
!!       FORCE = 0.0_SP
!!       DO J = 1, 6
!!         PHAI_IJ= PHAI(I,J)*PI2/360.0_SP
!!!         FORCE  = APT(I,J) * COS(PI2/PERIOD(J)*TIMERK1 -PHAI_IJ) + FORCE !!GWC
!!         FORCE  = APT(I,J) * COS(PI2/PERIOD(J)*TIME1 -PHAI_IJ) + FORCE
!!       END DO
!!       FORCE= FORCE + EMEAN(I)
!!       ELF(I_SBC_N(I)) = FORCE * RAMP
!     END DO
!   END IF

!
!-Radiation Boundary Conditions (Sommerfeld) for Free Surface Elevation---------|
!
!   DO J = 1, IRBCN
!     I1 = I_RBC_N(J)
!     I2 = NEXT_RBC(J)
!     CC = SQRT(GRAV*H(I1))*DTE/DLTN_RBC(J)
!     CP = CC + 1.

!     ELF(I1) =  (CC*ELF(I2)+EL(I1))/CP
!   END DO

   CALL BCOND_ASL
   CALL BCOND_ASL_CLP
   CALL BCOND_GWI
   CALL BCOND_BKI
   CALL BCOND_ORE
!>JQI 07/14/04

!
!--Allow setup/down on north boundary in response to longshore wind
!--Corrects for Wind-Driven Barotropic Response (See Schwing 1989)
!--Implemented by Jamie Pringle
!

!  CASEUNIQUE
   IF(CASENAME == "gom")THEN
     ETAXFER = 0.5_SP
     ANGALONG = 68.0_SP/360.0_SP*6.283185_SP
     DO J=1,IOBCN
       II = I_OBC_N(J)
       IGL = II
       IF(IGL == 1 .OR. IGL == 2)THEN
         TXJMP2 = WUSURF2(2)*1000.0_SP
         TYJMP2 = WVSURF2(2)*1000.0_SP
         ANGWIND=ATAN2(TYJMP2,TXJMP2)-ANGALONG
         WNDALONG=COS(ANGWIND)*SQRT(TXJMP2**2+TYJMP2**2)
         IF (IGL == 1) THEN
           ELF(II)=ELF(II)-WNDALONG*ETAXFER
         ELSEIF (IGL == 2) THEN
           ELF(II)=ELF(II)-0.5_SP*WNDALONG*ETAXFER
         ENDIF
       ENDIF
     END DO
   END IF
!  END CASEUNIQUE

!  CASEUNIQUE
   IF(CASENAME == "gom")THEN

     CUMEL=0.0_SP
     DO NBCJMP=1,NOBCGEO
       NNOW=IBCGEO(NBCJMP)

       !INTEGRATE RHO IN DEPTH, CONVERT TO MKS UNITS DAMIT
       RHOILST=RHOINT
       RHOINT=0.0_SP
       DO K=1,KBM1
         RHOINT=RHOINT+(1.0_SP+RHO1(NNOW,K))*1.0E3_SP*DZ(K)
       END DO

       !FIND DENSITY GRADIENT, AND MODIFY BOUNDARY ELEVATION
       !NOTE THE FACTOR OF 1000 AND 2 TO COMPENSATE FOR THE
       !FACT THAT THE MODEL STORES RHO1 AS SIGMA, AND IN CGS.
       IF (NBCJMP /= 1) THEN
         NLAST=IBCGEO(NBCJMP-1)
         DXBC=SQRT((VX(NNOW)-VX(NLAST))**2+(VY(NNOW)-VY(NLAST))**2)
         ETATAN=-(1.0_SP/(0.5_SP*(RHOINT+RHOILST))) &
               *((H(NNOW)*RHOINT-H(NLAST)*RHOILST)/DXBC &
               -0.5_SP*1.0e3_SP*(2.0_SP+RHO1(NNOW,KBM1)+RHO1(NLAST,KBM1)) &
                *(H(NNOW)-H(NLAST))/DXBC)
         CUMEL=CUMEL+ETATAN*DXBC
         ELF(NNOW)=ELF(NNOW)+CUMEL

       ENDIF
     END DO

  END IF
! END CASEUNIQUE


!==============================================================================|
   CASE(2) !External Mode Velocity Boundary Conditions                         |
!==============================================================================|

   DO I=1,N

!
!--2 SOLID BOUNDARY EDGES------------------------------------------------------|
!
!Q&C< commented on 06/22/04
!   IF(ISBCE(I) == 3) THEN
!     UAF(I)=0.0_SP
!     VAF(I)=0.0_SP
!   END IF

!  GWC SPEED REPLACE ABOVE 4 LINES
!   UAF(LISBCE_3(1:NISBCE_3)) = 0.
!   VAF(LISBCE_3(1:NISBCE_3)) = 0.

!
!--1 SOLID BOUNDARY EDGE-------------------------------------------------------|
!
   IF(ISBCE(I) == 1) THEN
     ALPHA1=ALPHA(I)
     IF(NUMQBC > 0) THEN
       IF(INFLOW_TYPE == 'node') THEN
         DO J=1,NUMQBC
           I1=INODEQ(J)
           J1=NBVE(I1,1)
           J2=NBVE(I1,NTVE(I1))
           IF((I == J1).OR.(I == J2)) THEN
             UNTMP=UAF(I)*COS(ANGLEQ(J))+VAF(I)*SIN(ANGLEQ(J))
             VNTMP=-UAF(I)*SIN(ANGLEQ(J))+VAF(I)*COS(ANGLEQ(J))
             UNTMP=MAX(UNTMP,0.0_SP)
             UAF(I)=UNTMP*COS(ANGLEQ(J))-VNTMP*SIN(ANGLEQ(J))
             VAF(I)=UNTMP*SIN(ANGLEQ(J))+VNTMP*COS(ANGLEQ(J))
               GOTO 21
               END IF
             END DO
            ELSE IF(INFLOW_TYPE == 'edge') THEN
             DO J=1,NUMQBC
               J1=ICELLQ(J)
               IF(I == J1) THEN
                UNTMP=UAF(I)*COS(ANGLEQ(J))+VAF(I)*SIN(ANGLEQ(J))
                VNTMP=-UAF(I)*SIN(ANGLEQ(J))+VAF(I)*COS(ANGLEQ(J))
                UNTMP=MAX(UNTMP,0.0_SP)
                UAF(I)=UNTMP*COS(ANGLEQ(J))-VNTMP*SIN(ANGLEQ(J))
                VAF(I)=UNTMP*SIN(ANGLEQ(J))+VNTMP*COS(ANGLEQ(J))
               GOTO 21
               END IF
             END DO
           END IF

           END IF

!Q&C< commented on 06/22/04
!           UI= UAF(I)*(SIN(ALPHA1))**2-VAF(I)*SIN(ALPHA1)*COS(ALPHA1)
!           VI=-UAF(I)*SIN(ALPHA1)*COS(ALPHA1)+VAF(I)*(COS(ALPHA1))**2
!           UAF(I)=UI
!           VAF(I)=VI

!           UI= UAS(I)*(SIN(ALPHA1))**2-VAS(I)*SIN(ALPHA1)*COS(ALPHA1)
!           VI=-UAS(I)*SIN(ALPHA1)*COS(ALPHA1)+VAS(I)*(COS(ALPHA1))**2
!           UAS(I)=UI
!           VAS(I)=VI
!           UAS(I)=UAF(I)
!           VAS(I)=VAF(I)
!Q&C> 06/22/04

21        CONTINUE
          END IF
       END DO


!==============================================================================|
   CASE(3) !3-D Velocity Boundary Conditions                                   !
!==============================================================================|

!  GWC SPEED REPLACE NEXT 4 LINES
!   UF(LISBCE_3(1:NISBCE_3),1:KBM1) = 0.
!   VF(LISBCE_3(1:NISBCE_3),1:KBM1) = 0.

       do i= 1, n
         do k =1, kbm1
!Q&C< commented on 06/22/04
!          if(isbce(i).eq.3) then
!           uf(i,k)=0.0_SP
!           vf(i,k)=0.0_SP
!          end if
!Q&C> 06/22/04

          if(isbce(i).eq.1) then
          alpha1=alpha(i)
          if(numqbc.ge.1) then
          if(inflow_type.eq.'node') then
            do j=1,numqbc
              i1=inodeq(j)
              j1=nbve(i1,1)
              j2=nbve(i1,ntve(i1))
              if((i.eq.j1).or.(i.eq.j2)) then
               untmp=uf(i,k)*cos(angleq(j))+vf(i,k)*sin(angleq(j))
               vntmp=-uf(i,k)*sin(angleq(j))+vf(i,k)*cos(angleq(j))
               untmp=max(untmp,0.0_SP)
               uf(i,k)=untmp*cos(angleq(j))-vntmp*sin(angleq(j))
               vf(i,k)=untmp*sin(angleq(j))+vntmp*cos(angleq(j))
              goto 31
              end if
            end do
           else if(inflow_type.eq.'edge') then
            do j=1,numqbc
              j1=icellq(j)
              if(i.eq.j1) then
               untmp=uf(i,k)*cos(angleq(j))+vf(i,k)*sin(angleq(j))
               vntmp=-uf(i,k)*sin(angleq(j))+vf(i,k)*cos(angleq(j))
               untmp=max(untmp,0.0_SP)
               uf(i,k)=untmp*cos(angleq(j))-vntmp*sin(angleq(j))
               vf(i,k)=untmp*sin(angleq(j))+vntmp*cos(angleq(j))
              goto 31
              end if
            end do
           else
             print*, 'inflow_type not correct'
             call pstop
          end if
          end if

!Q&C< commented on 06/22/04
!          ui= uf(i,k)*(sin(alpha1))**2-vf(i,k)*sin(alpha1)*cos(alpha1)
!          vi=-uf(i,k)*sin(alpha1)*cos(alpha1)+vf(i,k)*(cos(alpha1))**2
!          uf(i,k)=ui
!          vf(i,k)=vi

!          ui= us(i,k)*(sin(alpha1))**2-vs(i,k)*sin(alpha1)*cos(alpha1)
!          vi=-us(i,k)*sin(alpha1)*cos(alpha1)+vs(i,k)*(cos(alpha1))**2
!          us(i,k)=ui
!          vs(i,k)=vi
!          us(i,k)=uf(i,k)
!           vs(i,k)=vf(i,k)
!Q&C> 06/22/04

31        continue
          end if
         end do
       end do
!==============================================================================|
   CASE(4) !Blank                                                              !
!==============================================================================|
!==============================================================================|
   CASE(5) !!SOLID BOUNDARY CONDITIONS ON U AND V                              !
!==============================================================================|
!  GWC SPEED REPLACE NEXT 4 LINES
!   U(LISBCE_3(1:NISBCE_3),1:KBM1) = 0.
!   V(LISBCE_3(1:NISBCE_3),1:KBM1) = 0.

       do i= 1, n
         do k =1, kbm1
!Q&C< commented on 06/22/04
!          if(isbce(i).eq.3) then
!          u(i,k)=0.0_SP
!          v(i,k)=0.0_SP
!          end if
!Q&C> 06/22/04
          if(isbce(i).eq.1) then
          alpha1=alpha(i)
          if(numqbc.ge.1) then
          if(inflow_type.eq.'node') then
            do j=1,numqbc
              i1=inodeq(j)
              j1=nbve(i1,1)
              j2=nbve(i1,ntve(i1))
              if((i.eq.j1).or.(i.eq.j2)) then
               untmp=u(i,k)*cos(angleq(j))+v(i,k)*sin(angleq(j))
               vntmp=-u(i,k)*sin(angleq(j))+v(i,k)*cos(angleq(j))
               untmp=max(untmp,0.0_SP)
               u(i,k)=untmp*cos(angleq(j))-vntmp*sin(angleq(j))
               v(i,k)=untmp*sin(angleq(j))+vntmp*cos(angleq(j))
              goto 51
              end if
            end do
           else if(inflow_type.eq.'edge') then
            do j=1,numqbc
              j1=icellq(j)
              if(i.eq.j1) then
               untmp=u(i,k)*cos(angleq(j))+v(i,k)*sin(angleq(j))
               vntmp=-u(i,k)*sin(angleq(j))+v(i,k)*cos(angleq(j))
               untmp=max(untmp,0.0_SP)
               u(i,k)=untmp*cos(angleq(j))-vntmp*sin(angleq(j))
               v(i,k)=untmp*sin(angleq(j))+vntmp*cos(angleq(j))
              goto 51
              end if
            end do
           else
             print*, 'inflow_type not correct'
             call pstop
          end if
          end if

!Q&C< commented on 06/22/04
!          ui= u(i,k)*(sin(alpha1))**2-v(i,k)*sin(alpha1)*cos(alpha1)
!          vi=-u(i,k)*sin(alpha1)*cos(alpha1)+v(i,k)*(cos(alpha1))**2
!          u(i,k)=ui
!          v(i,k)=vi

!          ui= us(i,k)*(sin(alpha1))**2-vs(i,k)*sin(alpha1)*cos(alpha1)
!          vi=-us(i,k)*sin(alpha1)*cos(alpha1)+vs(i,k)*(cos(alpha1))**2
!          us(i,k)=ui
!          vs(i,k)=vi
!          us(i,k)=u(i,k)
!          vs(i,k)=v(i,k)
!Q&C> 06/22/04

51        continue
          end if
         end do
       end do
!==============================================================================|
   CASE(6) !Blank                                                              !
!==============================================================================|
!==============================================================================|
   CASE(7) !Blank                                                              !
!==============================================================================|
!==============================================================================|
   CASE(8) !!SURFACE FORCING FOR INTERNAL MODE                                 !
!==============================================================================|
!
!--Fresh Water Discharge-------------------------------------------------------|
!
   IF(NUMQBC > 0)THEN
     CALL BRACKET(QBC_TM,THOUR,L1,L2,FACT,UFACT,IERR)
     QDIS(:) = UFACT*DQDIS(:,L1) + FACT*DQDIS(:,L2)
     TDIS(:) = UFACT*DTDIS(:,L1) + FACT*DTDIS(:,L2)
     SDIS(:) = UFACT*DSDIS(:,L1) + FACT*DSDIS(:,L2)
     QDIS    = QDIS*RAMP

  if  (WATER_QUALITY==1)then
     DO N1 = 1, NB
       WDIS(:,N1) = UFACT*DWDIS(:,N1,L1) + FACT*DWDIS(:,N1,L2)
     END DO
  endif
   END IF
   IF(M_TYPE == 'uniform') THEN

     CALL BRACKET(UMF_TM,THOUR,L1,L2,FACT,UFACT,IERR)
!
!----Surface Evaporation and Precipitation---------------------------------------|
!
     QPREC = UFACT*DQPREC(L1) + FACT*DQPREC(L2)
     QEVAP = UFACT*DQEVAP(L1) + FACT*DQEVAP(L2)

!
!--- Heat Flux and Short Wave Radiation----------------------------------------!
!

     SPCP=4.2174E3_SP
     ROSEA = 1.023E3_SP
     SPRO = SPCP*ROSEA
     WTSURF(:) = UFACT*UHFLUX(L1)  + FACT*UHFLUX(L2)
     SWRAD(:)  = UFACT*UHSHORT(L1) + FACT*UHSHORT(L2)
     WTSURF    = -WTSURF/SPRO*RAMP
     SWRAD     = -SWRAD/SPRO*RAMP

!
!--- Wind Stress for the Internal Mode-----------------------------------------!
!

     TX = UFACT*UWIND(L1) + FACT*UWIND(L2)
     TY = UFACT*VWIND(L1) + FACT*VWIND(L2)

     IF(WINDTYPE == 'speed')THEN
       WDS=SQRT(TX*TX+TY*TY)
       CD=1.2E-3
       IF (WDS >= 11.0_SP) CD=(0.49_SP+0.065_SP*WDS)*1.E-3_SP
       IF (WDS >= 25.0_SP) CD=(0.49_SP+0.065_SP*25.0_SP)*1.E-3_SP
       TX=1.2_SP*CD*TX*WDS
       TY=1.2_SP*CD*TY*WDS
       UUWIND(:)=-1.0E-3_SP*TX
       VVWIND(:)=-1.0E-3_SP*TY
       WUSURF(:)=-1.0E-3_SP*TX*RAMP
       WVSURF(:)=-1.0E-3_SP*TY*RAMP
     ELSE IF(WINDTYPE == 'stress') THEN
       TX=0.001_SP*TX
       TY=0.001_SP*TY
       UUWIND(:)=-TX
       VVWIND(:)=-TY
       WUSURF(:)=-TX*RAMP
       WVSURF(:)=-TY*RAMP
     END IF

   END IF !! M_TYPE='uniform'


   IF(M_TYPE == 'non-uniform')THEN

!
!--- Heat flux and short wave radiation----------------------------------------!
!
     CALL BRACKET(HFX_TM,THOUR,L1,L2,FACT,UFACT,IERR)

     QPREC = 0.0_SP
     QEVAP = 0.0_SP

     SPCP  = 4.2174E3_SP
     ROSEA = 1.023E3_SP
     SPRO  = SPCP*ROSEA

     IF(IERR==-1)THEN
       WTSURF = 0.0_SP
       SWRAD  = 0.0_SP
     ELSE
       WTSURF(1:M) = UFACT*DHFLUX(1:M,L1)  + FACT*DHFLUX(1:M,L2)
       SWRAD(1:M)  = UFACT*DHSHORT(1:M,L1) + FACT*DHSHORT(1:M,L2)
       WTSURF = -WTSURF/SPRO*RAMP
       SWRAD  = -SWRAD/SPRO*RAMP
     END IF


!
!--- Wind Stress for the Internal Mode-----------------------------------------!
!

     CALL BRACKET(WND_TM,THOUR,L1,L2,FACT,UFACT,IERR)

     IF(IERR == -1)THEN
       WUSURF = 0.0_SP
       WVSURF = 0.0_SP
     ELSE

     IF(WINDTYPE == 'speed')THEN
       DO I=1,N
         TX = UFACT*DTX(I,L1) + FACT*DTX(I,L2)
         TY = UFACT*DTY(I,L1) + FACT*DTY(I,L2)
         WDS=SQRT(TX*TX+TY*TY)
         CD=1.2E-3_SP
         IF (WDS >= 11.0_SP) CD=(0.49_SP+0.065_SP*WDS)*1.E-3_SP
         IF (WDS >= 25.0_SP) CD=(0.49_SP+0.065_SP*25.0_SP)*1.E-3_SP
         TX = 1.2_SP*CD*TX*WDS
         TY = 1.2_SP*CD*TY*WDS
         UUWIND(I)=-1.0E-3_SP*TX
         VVWIND(I)=-1.0E-3_SP*TY
         WUSURF(I)=-1.0E-3_SP*TX*RAMP
         WVSURF(I)=-1.0E-3_SP*TY*RAMP
       END DO
     ELSE IF(WINDTYPE == 'stress') THEN
       DO I=1,N
         TX = UFACT*DTX(I,L1) + FACT*DTX(I,L2)
         TY = UFACT*DTY(I,L1) + FACT*DTY(I,L2)
         TX = 0.001_SP*TX
         TY = 0.001_SP*TY
         UUWIND(I) = -TX
         VVWIND(I) = -TY
         WUSURF(I) = -TX*RAMP
         WVSURF(I) = -TY*RAMP
       END DO
     END IF
     END IF

   END IF !! MTYPE='non-uniform'

!==============================================================================|
   CASE(9) !External Mode Surface BCs (River Flux/Wind Stress/Heat/Moist)      !
!==============================================================================|

!
!-Freshwater Flux: Set  Based on Linear Interpolation Between Two Data Times---|
!

   IF (NUMQBC /= 0) THEN
     CALL BRACKET(QBC_TM,THOUR1,L1,L2,FACT,UFACT,IERR)
     QDIS2(:) = UFACT*DQDIS(:,L1) + FACT*DQDIS(:,L2)
     QDIS2    = QDIS2*RAMP
   END IF

!
!--Uniform Meteorology -> Set Precipitation/Evaporation/Surface Wind-----------|
!

   IF(M_TYPE == 'uniform') THEN

     CALL BRACKET(UMF_TM,THOUR1,L1,L2,FACT,UFACT,IERR)

     IF(IERR == -1)THEN
       QPREC2=0.0_SP
       QEVAP2=0.0_SP
       WUSURF2 = 0.0_SP
       WVSURF2 = 0.0_SP
     ELSE

       QPREC2 = UFACT*DQPREC(L1) + FACT*DQPREC(L2)
       QEVAP2 = UFACT*DQEVAP(L1) + FACT*DQEVAP(L2)
       TX = UFACT*UWIND(L1) + FACT*UWIND(L2)
       TY = UFACT*VWIND(L1) + FACT*VWIND(L2)

       IF(WINDTYPE == 'speed')THEN
         WDS=SQRT(TX*TX+TY*TY)
         CD=1.2E-3_SP
         IF (WDS >= 11.0_SP) CD=(0.49_SP+0.065_SP*WDS)*1.E-3_SP
         IF (WDS >= 25.0_SP) CD=(0.49_SP+0.065_SP*25.0_SP)*1.E-3_SP
         TX      = 1.2_SP*CD*TX*WDS
         TY      = 1.2_SP*CD*TY*WDS
         UUWIND  = 1.0E-3_SP*TX
         VVWIND  = 1.0E-3_SP*TY
         WUSURF2 = 1.0E-3_SP*TX*RAMP
         WVSURF2 = 1.0E-3_SP*TY*RAMP
       ELSE IF(WINDTYPE == 'stress')THEN
         TX      = 0.001_SP*TX
         TY      = 0.001_SP*TY
         UUWIND  = TX
         VVWIND  = TY
         WUSURF2 = TX*RAMP
         WVSURF2 = TY*RAMP
       END IF
       END IF

   END IF !!M_TYPE = 'uniform'

!
!--Non-Uniform Meteorology -> Set Precipitation/Evaporation/Surface Wind-------|
!

   IF(M_TYPE == 'non-uniform') THEN

     QPREC2=0.0_SP
     QEVAP2=0.0_SP

     CALL BRACKET(WND_TM,THOUR1,L1,L2,FACT,UFACT,IERR)

     IF(IERR == -1)THEN
       WUSURF2 = 0.0_SP
       WVSURF2 = 0.0_SP
     ELSE
     IF(WINDTYPE == 'speed') THEN
!       CD=1.2E-3_SP
       DO I=1,N
         TX = UFACT*DTX(I,L1) + FACT*DTX(I,L2)
         TY = UFACT*DTY(I,L1) + FACT*DTY(I,L2)
         WDS=SQRT(TX*TX+TY*TY)
         CD=1.2E-3_SP
         IF (WDS >= 11.0_SP) CD=(0.49_SP+0.065_SP*WDS)*1.E-3_SP
         IF (WDS >= 25.0_SP) CD=(0.49_SP+0.065_SP*25.0_SP)*1.E-3_SP
         TX = 1.2_SP*CD*TX*WDS
         TY = 1.2_SP*CD*TY*WDS
         UUWIND(I)  = 1.0E-3_SP*TX
         VVWIND(I)  = 1.0E-3_SP*TY
         WUSURF2(I) = 1.0E-3_SP*TX*RAMP
         WVSURF2(I) = 1.0E-3_SP*TY*RAMP
       END DO

     ELSE IF(WINDTYPE == 'stress') THEN
       DO I=1,N
         TX = UFACT*DTX(I,L1) + FACT*DTX(I,L2)
         TY = UFACT*DTY(I,L1) + FACT*DTY(I,L2)
         TX = 0.001_SP*TX
         TY = 0.001_SP*TY
         UUWIND(I) = TX
         VVWIND(I) = TY
         WUSURF2(I) = TX*RAMP
         WVSURF2(I) = TY*RAMP
       END DO
     END IF
     END IF


   END IF !!MTYPE='non-uniform'

   END SELECT

   RETURN
   END SUBROUTINE BCOND_GCY
!==============================================================================|


