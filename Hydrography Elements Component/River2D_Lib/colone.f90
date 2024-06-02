                     SUBROUTINE   cyclone
                      USE ALL_VARS
					  USE BCS

                      REAL(SP) ::   pre0,pre1,  r,R0, beta, theta,  c1, c2, dpdr,w1
					  real(sp) ::    uvc0, vc0, uc0, xc0, yc0, alpha1
					  real(sp)::    xx(20), yy(20), pp(20),lat(20),ddt, dts(20)
					  INTEGER :: I,K,J,NT1
                      logical  check
					 

                     INQUIRE(FILE='../cyclone/cyclone.dat',EXIST=CHECK) 
                     if(.not.check) then
					 stop 'cyclone.dat does not exist' 
                     else

                    open(121, file= '../cyclone/cyclone.dat',status='old')
                    read(121,*)  nt1     ,ddt 
					read(121,*) ihour
					ihour=ihour-1
					do i=1,nt1
					read(121,*) xx(i), yy(i),  pp(i) , lat(i)       ! 台风路径，中心压强,   纬度
					xx(i)=xx(i)-vxmin
					yy(i)=yy(i)-vymin
					enddo
				    close(121)
					endif

					  
					 beta=30*3.1415926/180.                !  偏角
                     c1=0.8                                 !  订正系数
                     c2=0.5
                     R0=40000.
                     pre1=1013.3   ! hpa                   !  外围气压

                     WND_TM%NTIMES=nsteps*dti/3600            ! 每小时计算一次
					
					 ALLOCATE(WND_TM%TIMES(WND_TM%NTIMES))
                     do i=1,WND_TM%NTIMES
                     WND_TM%TIMES(i)=1.*(i-1)
                     enddo  

					 ALLOCATE(DTX(N,WND_TM%NTIMES),DTY(N,WND_TM%NTIMES),pre(n,wnd_tm%ntimes))


					  Dtotal=0.  
                      do i=1,nt1-1 
					  Dts(i)=sqrt((xx(i+1)-xx(i))**2+(yy(i+1)-yy(i))**2)
					  enddo

					  do k=1, WND_TM%NTIMES
					    if(k>ihour) then 
					    j=(k-ihour)*3600./ddt+1          !  j ,  j+1
						if(j+1>nt1)  goto 320
						alpha1=mod((k-ihour)*3600., ddt)/ddt
		!				if(alpha1==0) alpha1=1.
						uvc0=Dts(j)/ddt
						theta=atan2(yy(j+1)-yy(j), xx(j+1)-xx(j))
						uc0=uvc0*cos(theta)
						vc0=uvc0*sin(theta)
						xc0=  (xx(j+1)-xx(j))*alpha1 +xx(j)
					    yc0= (yy(j+1)-yy(j))*alpha1 +yy(j)
					    pre0= (pp(j+1)-pp(j))*alpha1 +pp(j)
						R0=28.52*tanh(0.0873*(lat(j)-28))+12.22/exp((pre1-pp(j))*100/33.86)+0.2*uvc0+37.22
					    R0=R0*100. 
					  do i=1, n
                        
					    r=sqrt((xc(i)-xc0)**2+(yc(i)-yc0)**2)
						pre(i,k)=pre1*100-(pre1-pre0)*100./sqrt(1.+(r/R0)**2)
                        dpdr=100.*(pre1- pre0)*(1+(r/R0)**2)**(-1.5)*(r/R0)/r0
                        w1=sqrt(cor(i)**2*r**2/4.+r/1.025*dpdr)-0.5*cor(i)*r
                        theta=atan2(yc(i)-yc0, xc(i)-xc0)
                        if(theta<0) theta=theta+6.28
                        dtx(i, k)=-w1*sin(theta+beta)*c1+ uc0*exp(-3.14*r/500000.)*c2
                        dty(i, k)=w1*cos(theta+beta)*c1+vc0*exp(-3.14*r/500000.)*c2

                     enddo   !  i 
                          else      ! if(k<ihour)
						  do i=1,n
                        dtx(i, k)=0.
                        dty(i, k)=0.
						pre(i,k)=0.
						  enddo


					 endif  !k 

                     enddo   !k  

!  *********************** write to tecplot for check wind  field*****************************************
320                  continue 
 !                    open(121,file='../output/cyclone.dat')
 !			   		  do k=1,  WND_TM%NTIMES
                   
 !                    write(121,*)'VARIABLES = "X","Y","P","Wu","Wv" '
 !	                 write(121,*)'Zone n=', MGL,',e=',ngl,', f=feblock,et=triangle'    !  quadrilateral
 !	                 write(121,*)'VARLOCATION = ([3-5]=CELLCENTERED)'

 !	                 write(121,"(200f18.4)") (xg(i),i=1,mgl)
 !                    write(121,"(200f18.4)") (yg(i),i=1,mgl)
 !	                 write(121,"(200f18.4)") (pre(i,k),i=1,ngl)
 !	                 write(121,"(200f18.4)") (dtx(i,k),i=1,ngl)
 !	                 write(121,"(200f18.4)") (dty(i,k),i=1,ngl)



 !                    	 do i=1,ngl
 !	                     write(121,*) (nv(i,j),j=1,3)
 !	                     enddo


 !                     enddo
 !                     close(121)

					 END SUBROUTINE cyclone