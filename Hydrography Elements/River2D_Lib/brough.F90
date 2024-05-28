!==============================================================================|
!  Calculate Bottom Drag Coefficient based on Bottom Roughness                 !
!   note:                                                                      !
!   when the log function derived from the constant stress log-viscous         !
!   layer is applied to an estuary, if the value of z0 is close to             !
!   (zz(kbm1)-z(kb)*dt1, drag coefficient "cbc" could become a huge            !
!   number due to near-zero value of alog function. In our application         !
!   we simply cutoff at cbc=0.005. One could adjust this cutoff value          !
!   based on observations or his or her experiences.                           !   
!   CALCULATES:   WUBOT(N), WVBOT(N) : BOTTOM SHEAR STRESSES                   !
!==============================================================================|

   SUBROUTINE BOTTOM_ROUGHNESS

!==============================================================================!
   USE ALL_VARS
   USE MOD_WD
   use wave
   IMPLICIT NONE
   INTEGER :: I,II, n1,n2,n3
   REAL(SP), PARAMETER  :: KAPPA = .40_SP   !!VON KARMAN LENGTH SCALE
   REAL(SP), PARAMETER  :: VK2   = .160_SP  !!KAPPA SQUARED
   REAL(SP)             :: CBCMIN,Z0,ZTEMP,BTPS,RR,U_TAUB,Z0B_GOTM,Z0B_TEMP
   real(sp)             :: Twave, Hwave, theta, Ab,uw,vw, wavek,wavekn,ee,domga,fw
!==============================================================================!


!
!  SET CONSTANTS
!
   CBCMIN = BFRIC
   Z0     = Z0B

!==============================================================================|
   IF(BROUGH_TYPE == 'orig')THEN !USE ORIGINAL FVCOM FORM FOR BOTTOM FRICTION  |
!==============================================================================|
     DO I=1,N
       IF(DT1(I) > 3.0)THEN
        ZTEMP=(ZZ(KBM1)-Z(KB))*DT1(I)/Z0
        CBC(I) = MAX(CBCMIN,VK2/(LOG(ZTEMP))**2)
       ELSE 
        ZTEMP=(ZZ(KBM1)-Z(KB))*3.0/Z0
        CBC(I) = MAX(CBCMIN,VK2/(LOG(ZTEMP))**2)
       END IF

     END DO
!==============================================================================|
   ELSE IF(BROUGH_TYPE == 'gotm')THEN !GOTM FORMULATION FOR BOTTOM FRICTION    |
!==============================================================================|

!----Convert Input Z0B to GOTMS H0B
     Z0B_TEMP = Z0B/.03  

     DO I=1,N
     U_TAUB = 0.
     DO II=1,40       
        IF (UMOL <= 0.) THEN
           Z0B_GOTM=0.03*Z0B  
        ELSE
           Z0B_GOTM=0.1*UMOL/MAX(UMOL,U_TAUB)+0.03*Z0B_TEMP
        END IF
          ztemp=(zz(kbm1)-z(kb))*dt1(i)
        RR=KAPPA/(LOG((Z0B_GOTM+ZTEMP)/Z0B_GOTM))
      U_TAUB = RR*SQRT( U(I,KBM1)*U(I,KBM1) + V(I,KBM1)*V(I,KBM1) )
     END DO
     CBC(I) =   RR*RR
     END DO

!==============================================================================|
   ELSE IF(BROUGH_TYPE == 'user_defined')THEN !Use User Defined broud_ud.F     | 
!==============================================================================|
   
     CALL BROUGH_UD

   END IF


!==============================================================================|
!  CALCULATE SHEAR STRESS ON BOTTOM  --> WUBOT/WVBOT                           |
!==============================================================================|
  DO  I = 1, N
     IF(D1(I) > 0.0_SP) THEN
      if  (TWO_D_MODEL==0)then     
       BTPS= CBC(I)*SQRT(U(I,KBM1)**2+V(I,KBM1)**2)
       WUBOT(I) = -BTPS * U(I,KBM1)
       WVBOT(I) = -BTPS * V(I,KBM1)


!	  n1=nv(i,1)
!	  n2=nv(i,2)
!	  n3=nv(i,3)
!     if(iwave==1.and.wave_h(n1)*wave_h(n2)*wave_h(n3)>0.) then
    
!	  Twave=0.3333*(wave_t(n1)+ wave_t(n2) + wave_t(n3))
!	  Hwave=0.3333*(wave_H(n1)+ wave_H(n2) + wave_H(n3))
!	  theta=0.3333*(wave_theta(n1)+ wave_theta(n2)+ wave_theta(n3))
!      domga=6.28/Twave

!	  wavekn=domga**2/9.81
!      ee=1.
!		   do	while (ee>=1e-5)
!	       wavek=wavekn+(domga**2-9.81*wavekn*tanh(wavekn*d1(i)))/(9.81*tanh(wavekn*d1(i))+9.81*wavekn/(cosh(wavekn*d1(i)))**2)
!           ee=abs(wavek-wavekn)
!		   wavekn=wavek
!           if(ee>100.)stop'wave number not convergence'	
!	       enddo
!      Ab=Hwave/sinh(wavek*d1(i))
!      uw=Ab*domga*cos((270-theta)*3.1415926/180)
!	  vw=Ab*domga*sin((270-theta)*3.1415926/180)
!	if(30*Z0/Ab<0.08)then
!	fw=0.13*(30*Z0/Ab)**0.4
!	else if(30*Z0/Ab<1)then
!	fw=0.23*(30*Z0/Ab)**0.62
!	else
!	fw=0.23
!	endif


!    WUBOT(I)=  -(sqrt(abs(WUBOT(I)))+0.026*Ab*domga)**2*sign(1.,u(i,kbm1))
!    WVBOT(I)=  -(sqrt(abs(WVBOT(I)))+0.026*Ab*domga)**2*sign(1.,v(i,kbm1))
     
	



      else
       BTPS= CBC(I)*SQRT(UA(I)**2+VA(I)**2)
       WUBOT(I) = -BTPS * UA(I)
       WVBOT(I) = -BTPS * VA(I)
      endif       
     ELSE
       WUBOT(I) = 0.0_SP
       WVBOT(I) = 0.0_SP
     END IF
   END DO
   
   

   RETURN
   END SUBROUTINE BOTTOM_ROUGHNESS
!==============================================================================|
