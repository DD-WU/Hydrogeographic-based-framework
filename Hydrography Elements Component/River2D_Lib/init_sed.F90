!=======================================================================
! FVCOM Sediment Initialization Routine (User Defined)
!   Configuration:  default    
!=======================================================================
  Subroutine Init_Sed  
  use all_vars
  Use Mod_Prec 
  Use Mod_Sed
  use wave
  Use Lims, only: m,kbm1
  implicit none 
  integer :: i,k,ised
  real(sp) :: bed_thickness, chezy, dd,uvm, k0,Ab,fw
  real(sp),allocatable :: concd( :,: )

  !--------------------------------------------------
  !Initialize Bed Properties
  !--------------------------------------------------
  bed_thickness = 10.00 !meters

  Do k=1,Nbed
    Do i=1,m
       bed(i,k,iaged) = float(iint+1)*dti                         !  床龄初始化
       bed(i,k,ithck) = bed_thickness/float(Nbed)
       bed(i,k,iporo) = 0.40
     End Do
   End do

  !--------------------------------------------------
  !Initialize Bed_Frac properties
  !--------------------------------------------------

  Do k=1,Nbed
    Do i=1,m
       Do ised=1,Nsed
         sed(ised)%frac(i,k)=1.0/float(nsed)   
       End Do
     End Do
  End Do

  !--------------------------------------------------
  !Initialize Sediment Concentrations
  !--------------------------------------------------
  
  allocate (concd(0:n,1:kb))

 call  BOTTOM_ROUGHNESS    !  obtain  CBC
!if(iwave==1) then
!do i=1, 4*12
!call readwave(2)                       
!enddo
!call wavechazhi
!call wavepump
!endif
do ised=1, nsed

  do k=1,kbm1
    do i=1,n
!    采用挟沙力公式进行初始浓度给定
    uvm=sqrt(ua(i)**2+va(i)**2)
    Chezy=sqrt(9.8/cbc(i))
	dd=d1(i)
	if(d1(i)<=0.1) dd=0.1
	concd(i,k)=0.023*1000*2.65/1.65*uvm**3/Chezy**2/dd/(sed(ised)%Wset*0.001) *0.33333
   
    enddo
  enddo
  

         
        call  E2N3D(concd,sed(ised)%conc)
     do k=1,kbm1
	 do i=1,m
	if(iwave==1)then
	if(wave_T(i)/=0) then
	dd=d(i)
	if(dd<=0.5) dd=0.5

     k0=0.001
	 Ab=0.5*wave_h(i)/sinh(wave_k(i)*dd)
    if(k0/Ab<0.08)then
	fw=0.13*(k0/Ab)**0.4
	else if(k0/Ab<1)then
	fw=0.23*(k0/Ab)**0.62
	else
	fw=0.23
	endif


!	sed(ised)%conc(i,k)=Sed(ised)%force(i)              &
!	  +0.004*10000*2.65/1.65*fw*wave_h(i)**3/(wave_T(i)**3*9.81*dd*sed(ised)%Wset*0.001*sinh(dd*wave_k(i))**3)
    endif
	endif
	enddo
	enddo


	  sed(ised)%cnew = min(sed(ised)%conc,  50.) 

      End Do
!   call tecplot
!   stop
  deallocate (concd)
  Return

  End Subroutine Init_Sed
!=======================================================================
