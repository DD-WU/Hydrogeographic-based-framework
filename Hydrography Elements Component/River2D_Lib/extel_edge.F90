!==============================================================================|
!  CALCULATE FLUXES OF FREE SURFACE ELEVATION (CONTINUITY) EQUATION            |
!==============================================================================|
    SUBROUTINE EXTEL_EDGE(K)       
!==============================================================================|
    USE ALL_VARS
    USE BCS
    USE MOD_OBCS
    USE MOD_NORTHPOLE
!QXU{
    USE MOD_BALANCE_2D
!QXU}
    IMPLICIT NONE
    REAL(SP) :: XFLUX(0:MT)
    REAL(SP) :: DIJ,UIJ,VIJ,DTK,UN,EXFLUX
    INTEGER  :: I,J,K,I1,IA,IB,JJ,J1,J2
!QXU{
    REAL(SP), DIMENSION(0:MT) :: XFLUXU,XFLUXV
    REAL(SP), DIMENSION(0:NT) :: XFLUXU1,XFLUXV1
    REAL(SP) :: EXFLUX_U,EXFLUX_V
!QXU}
!==============================================================================|

!----------INITIALIZE FLUX ARRAY ----------------------------------------------!

    XFLUX = 0.0_SP
!QXU{
    if  (BALANCE_2D==1)then

        XFLUXU= 0.0_SP
        XFLUXV= 0.0_SP
        XFLUXU1= 0.0_SP
        XFLUXV1= 0.0_SP
    endif
!QXU}
!---------ACCUMULATE FLUX BY LOOPING OVER CONTROL VOLUME HALF EDGES------------!

   DO I=1,NCV
     I1  = NTRG(I)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)
     DIJ = D1(I1)

     UIJ = UA(I1)
     VIJ = VA(I1)
     EXFLUX = DIJ*(-UIJ*DLTYE(I) + VIJ*DLTXE(I))  
     XFLUX(IA) = XFLUX(IA)-EXFLUX
     XFLUX(IB) = XFLUX(IB)+EXFLUX
!QXU{
    if  (BALANCE_2D==1)  then
         EXFLUX_U = -DIJ*UIJ*DLTYE(I)
         EXFLUX_V =  DIJ*VIJ*DLTXE(I)
         XFLUXU(IA) = XFLUXU(IA)-EXFLUX_U
         XFLUXU(IB) = XFLUXU(IB)+EXFLUX_U
         XFLUXV(IA) = XFLUXV(IA)-EXFLUX_V
         XFLUXV(IB) = XFLUXV(IB)+EXFLUX_V
    endif
!QXU}
    END DO

    if  (SPHERICAL==1)then
    if  (NORTHPOLE==1)then
        CALL EXTEL_EDGE_XY(K,XFLUX)
    endif
    endif
   
!--SAVE ACCUMULATED FLUX ON OPEN BOUNDARY NODES AND ZERO OUT OPEN BOUNDARY FLUX!

    IF(IOBCN > 0) THEN  
        DO I=1,IOBCN
            XFLUX_OBCN(I)=XFLUX(I_OBC_N(I))
            XFLUX(I_OBC_N(I)) = 0.0_SP
        END DO
    END IF

!---------ADJUST FLUX FOR FRESH WATER DISCHARGE--------------------------------!

    IF(NUMQBC >= 1) THEN   
        IF(INFLOW_TYPE == 'node') THEN
        DO J=1,NUMQBC
            JJ=INODEQ(J)
            XFLUX(JJ)=XFLUX(JJ)-QDIS(J)
        END DO
    ELSE IF(INFLOW_TYPE == 'edge') THEN
       DO J=1,NUMQBC
            J1=N_ICELLQ(J,1)
            J2=N_ICELLQ(J,2)
            XFLUX(J1)=XFLUX(J1)-QDIS(J)*RDISQ(J,1)
            XFLUX(J2)=XFLUX(J2)-QDIS(J)*RDISQ(J,2)
        END DO
    END IF
    END IF


!----------PERFORM UPDATE ON ELF-----------------------------------------------!
    DTK = ALPHA_RK(K)*DTE
    ELF = ELRK - DTK*XFLUX/ART1

!QXU{
!
!--STORE VARIABLES FOR MOMENTUM BALANCE CHECK----------------------------------|
!
    if(BALANCE_2D==1)then

    DO I=1,N
        XFLUXU1(I)=ONE_THIRD*(XFLUXU(NV(I,1))+ XFLUXU(NV(I,2))+ XFLUXU(NV(I,3)))
        XFLUXV1(I)=ONE_THIRD*(XFLUXV(NV(I,1))+ XFLUXV(NV(I,2))+ XFLUXV(NV(I,3)))
    END DO
   IF(K == 4) THEN
     DIVX2D2 = DIVX2D2 + XFLUXU1/ART/FLOAT(ISPLIT)            !dUD/dx
     DIVY2D2 = DIVY2D2 + XFLUXV1/ART/FLOAT(ISPLIT)            !dVD/dy
     DEDT2   = DEDT2 + (ELF1-ELRK1*(H1+ELRK1)/(H1+ELF1))/DTE/FLOAT(ISPLIT)    
    END IF     
    endif
!QXU}
    RETURN
    END SUBROUTINE EXTEL_EDGE
!==============================================================================|

