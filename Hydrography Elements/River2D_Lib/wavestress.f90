            subroutine  wavestress
            USE ALL_VARS

            REAL(4), ALLOCATABLE :: RTEMP_H(:),RTEMP_T(:),RTEMP_theta(:), Rtempxx(:), Rtempyy(:)
            g=9.81
           open(121,file='..\run\wave.dat')
           read(121,*)num
           allocate(RTEMP_H(num))
		   allocate(RTEMP_T(num))
		   allocate(RTEMP_theta(num))

		   do i=1,  num
		   read(121,*)Rtempxx(i),Rtempyy(i), RTEMP_H(i), RTEMP_T(i), RTEMP_theta(i)
		   enddo
           close(121)



          
		  
		  
		  
		  
		  
		   do i=1,ne     !  每一条边
		   do j=1,num
           rij=(xijc(i)-Rtempxx(j))**2+(yijc(i)-Rtempyy(j))**2
           Rtotal=Rtotal+1./rij 
	       enddo
		   do j=1,num
		   rij=(xijc(i)-Rtempxx(j))**2+(yijc(i)-Rtempyy(j))**2
           alpha=1./rij/Rtotal
           wave_h(i)=alpha*Rtemp_H(j)+wave_h(i)
           wave_T(i)=alpha*Rtemp_T(j)+wave_T(i)
           wave_theta(i)=alpha*Rtemp_theta(j)+wave_theta(i)
		   enddo
!          enddo    i
!-----------------------------------------------------输出波浪要素检查一下


!-----------------------------------------------------
           domga=2*3.1415926/wave_T(i)
           wavekn=domga**2/9.81 
	       ee=1.
		   do	while (ee>=1e-7)
	       wavek=wavekn+(domga**2-g*kn*tanh(kn*h))/(g*tanh(kn*h)+g*kn/(cosh(kn*h))**2)
           ee=abs(k-kn)
           if(ee>100.)stop'wave number not convergence'	
	       enddo
           n1=ienode(i,1)
		   n2=ienode(i,2)
           hij=0.5*(d1(n1)+d1(n2))
           if(hij<=0.5)then
		    cgbc=1
		   else
		    cgbc=0.5+k*hij/sinh(k*hij)
		   endif                                                            !  计算辐射应力
		   E=0.125*g*wave_h(i)**2 
	       radxxij(i)=E*(cgbc*cos(wave_theta(i))**2+(cgbc-0.5))
		   radyyij(i)=E*(cgbc*sin(wave_theta(i))**2+(cgbc-0.5))
		   radxyij(i)=E*cgbc*sin(wave_theta(i))*cos(wave_theta(i))

           enddo









			end