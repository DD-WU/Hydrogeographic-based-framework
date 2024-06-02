            MODULE wave


            REAL(4), ALLOCATABLE :: wave_h(:)               !!GLOBAL X-COORD AT NODE 
            REAL(4), ALLOCATABLE :: wave_T(:)      
            REAL(4), ALLOCATABLE :: wave_theta(:)
			real(4), allocatable :: wave_k(:)

            real(4), allocatable :: radxxij(:)
            real(4), allocatable :: radyyij(:)
		    real(4), allocatable :: radxyij(:)

			REAL(4), ALLOCATABLE :: RTEMP_H(:),RTEMP_T(:),RTEMP_theta(:), Rtempxx(:), Rtempyy(:)
            integer   num

		    CONTAINS
            
			SUBROUTINE ALLOC_VARS_WAVE
			USE ALL_VARS
            ALLOCATE(wave_h(0:mgl))           ;wave_h    = 0
            ALLOCATE(wave_T(0:mgl))            ;wave_T    = 0
            ALLOCATE(wave_theta(0:mgl))       ;wave_theta    = 0
			allocate(wave_k(0:mgl))           ; wave_k  =0 
			allocate(radxxij(1:ne))           ; radxxij = 0
			allocate(radyyij(1:ne))           ; radyyij = 0
			allocate(radxyij(1:ne))           ; radxyij = 0   
			end SUBROUTINE ALLOC_VARS_WAVE
!-------------------------------------------------------------------------------------------------------------			
			subroutine  readwave(idx)
            USE ALL_VARS
           integer idx
           select case (idx)
           case  (1) 
           num=30351
           allocate(Rtempxx(num))
           allocate(Rtempyy(num))
           allocate(RTEMP_H(num))
		   allocate(RTEMP_T(num))
		   allocate(RTEMP_theta(num))
           open(221,file='..\wave\xycor.txt')
		   do i=1,num  
		   read(221,*)Rtempxx(i),Rtempyy(i)   !, RTEMP_H(i), RTEMP_T(i), RTEMP_theta(i)
		   
          Rtempxx(i)=Rtempxx(i)-VXMIN
          Rtempyy(i)=Rtempyy(i)-vymin
		   enddo


           close(221)
           open(223,file='..\wave\tm01.txt')
!		   read(223,*) (Rtemp_t(i),i=1,num)
           open(224,file='..\wave\hs.txt')
!           read(224,*) (Rtemp_h(i), i=1,num)
		   open(225,file='..\wave\dir.txt') 
!		   read(225,*)(Rtemp_theta(i),i=1,num)
           do j=1, 24*4 
           read(223,*) (Rtemp_t(i),i=1,num)

           read(224,*) (Rtemp_h(i), i=1,num)

            read(225,*)(Rtemp_theta(i),i=1,num)
           enddo
  !         open(110,file='..\output\wave.dat')
  !         do i=1,num
  !		  if(Rtemp_h(i)>0) write(150,*) Rtempxx(i), Rtempyy(i), Rtemp_h(i)
  !		  enddo
     
		   case(2)
		   read(223,*) (Rtemp_t(i),i=1,num)
           read(224,*) (Rtemp_h(i), i=1,num)
           read(225,*)(Rtemp_theta(i),i=1,num)



          end select

	       end subroutine readwave

!-------------------------------------------------------------------------------------------------		  
		   subroutine  wavechazhi
           use all_vars
		   g=9.81
		   do i=1,m     !  每一个结点
		   Rtotal=0.
		   do j=1,num
		   if(Rtemp_t(j)>0.and.Rtemp_h(j)>0.and.Rtemp_theta(j)>0)then
           rij=(xg(i)-Rtempxx(j))**2+(yg(i)-Rtempyy(j))**2      !  xijc
           Rtotal=Rtotal+1./rij 
		   endif
	       enddo
		   wave_h(i)=0.
		   wave_t(i)=0.
		   wave_theta(i)=0.
		   do j=1,num
		   if(Rtemp_t(j)>0.and.Rtemp_h(j)>0.and.Rtemp_theta(j)>0)then
		   rij=(xg(i)-Rtempxx(j))**2+(yg(i)-Rtempyy(j))**2
           alpha1=1./rij/Rtotal
           wave_h(i)= alpha1*Rtemp_H(j) +wave_h(i)
           wave_T(i)=alpha1*Rtemp_T(j)+wave_T(i)
           wave_theta(i)=alpha1*Rtemp_theta(j)+wave_theta(i)
		   if(yg(i)>=60000) wave_h(i)=0
		   endif
		   enddo       !  j  =1, num
          enddo   ! i
           do i=1,m
		   domga=2*3.14/wave_t(i)
           wavekn=domga**2/9.81 
	       ee=1.
           hij=d(i)
		   if(hij<=0.1) hij=0.1
		   itt=0
		   do	while (ee>=1e-5)
	       wave_k(i)=wavekn+(domga**2-9.81*wavekn*tanh(wavekn*hij))/(9.81*tanh(wavekn*hij)+9.81*wavekn/(cosh(wavekn*hij))**2)
           ee=abs(wave_k(i)-wavekn)
		   wavekn = wave_k(i)
           itt=itt+1
           if(ee>100..or.itt>1000)stop'wave number not convergence'	
	       enddo


          enddo

		  end subroutine wavechazhi
!-----------------------------------------------------输出波浪要素检查一下
 !        call tecplot
          subroutine wavepump 
         USE ALL_VARS

      open(110,file='..\output\wave.dat')
	    write(110,*)'VARIABLES = "X","Y","waveh","waveT","wavetheta"'
	 write(110,*)'Zone n=', MGL,',e=',ngl,', f=feblock,et=triangle'    !  quadrilateral
!	 write(112,*)'VARLOCATION = ([6-9]=CELLCENTERED)'
  
	 write(110,"(200f18.4)") (xg(i),i=1,mgl)
     write(110,"(200f18.4)") (yg(i),i=1,mgl)
	 write(110,"(200f18.4)") (wave_h(i),i=1,mgl)
	 write(110,"(200f18.4)") (wave_T(i),i=1,mgl)
     write(110,"(200f18.4)") (wave_theta(i),i=1,mgl) 


	 do i=1,ngl
	 write(110,*) (nv(i,j),j=1,3)
	 enddo
      close(110)

          end subroutine wavepump  

!-----------------------------------------------------
           subroutine  wavestress
           use all_vars
		   integer itt
         
		   do i=1,ne
		

   		   n1=ienode(i,1)
   		   n2=ienode(i,2)
           Tij=0.5*(wave_T(n1)+wave_T(n2))
           domga=2*3.1415926/Tij
          
		   thetaij=(270.-0.5*(wave_theta(n1)+wave_theta(n2)))*3.1415926/180.
		   wavehij=0.5*(wave_h(n1)+wave_h(n2))

           hij=0.5*(d(n1)+d(n2))
		   wavekn=domga**2/9.81 
	       ee=1.
		   if(hij<=0.1) hij=0.1
		   itt=0
		   do	while (ee>=1e-5)
	       wavek=wavekn+(domga**2-9.81*wavekn*tanh(wavekn*hij))/(9.81*tanh(wavekn*hij)+9.81*wavekn/(cosh(wavekn*hij))**2)
           ee=abs(wavek-wavekn)
		   wavekn = wavek
           itt=itt+1
           if(ee>100..or.itt>1000)stop'wave number not convergence'	
	       enddo
          
           if(hij<=0.5)then
		    cgbc=1
		   else
		    cgbc=0.5+wavek*hij/sinh(wavek*hij)
		   endif                                                            !  计算辐射应力
		   E=0.125*9.81*wavehij**2 
	       radxxij(i)=E*(cgbc*cos(thetaij)**2+(cgbc-0.5))
		   radyyij(i)=E*(cgbc*sin(thetaij)**2+(cgbc-0.5))
		   radxyij(i)=E*cgbc*sin(thetaij)*cos(thetaij)

           enddo

		   end SUBROUTINE wavestress





					END MODULE wave
