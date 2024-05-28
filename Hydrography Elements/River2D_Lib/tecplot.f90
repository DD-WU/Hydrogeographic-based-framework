     subroutine tecplot
	  USE ALL_VARS
	  use mod_sed
 
    ! open(112, file=TRIM(OUTDIR)//"/"//trim(casename)//'_topo.dat')
	 if(sediment==1)then
	 write(112,*)'VARIABLES = "X","Y","h","d","sed","dep","u","v","Wu","Wv"'
	 else
	 write(112,*)'VARIABLES = "X","Y","h","el","u","v"'
     endif
	 write(112,*)'Zone n=', MGL,',e=',ngl,', f=feblock,et=triangle'    !  quadrilateral
	 if(sediment==1)then
	 write(112,*)'VARLOCATION = ([7-10]=CELLCENTERED)'
     else
	 write(112,*)'VARLOCATION = ([5-6]=CELLCENTERED)'
	 endif
	 write(112,"(200f18.4)") (xg(i),i=1,mgl)
     write(112,"(200f18.4)") (yg(i),i=1,mgl)
	 write(112,"(200f18.4)") (hg(i),i=1,mgl)
	 write(112,"(200f18.4)") (d(i),i=1,mgl)
	 if(sediment==1)then
	 write(112,"(200f18.4)") (sed(1)%conc(i,1),i=1,mgl)
	 write(112,"(200f18.4)") (bottom(i,dthck), i=1,mgl)
     endif
	 write(112,"(200f18.4)") (ua(i),i=1,ngl)
	 write(112,"(200f18.4)") (va(i),i=1,ngl)
     if(sediment==1)then
	 write(112,"(200f18.4)") (1000*wubot(i),i=1,ngl)
	 write(112,"(200f18.4)") (1000*wvbot(i),i=1,ngl)
     endif





	 do i=1,ngl
	 write(112,*) (nv(i,j),j=1,3)
	 enddo




!	 close(112)

!    write(121,"(10f8.3)")  el(14009), el(13612), el(11798), el(13775), el(13933),el(17613) !, 竹银， 黄金， 灯笼山，官冲， 西炮台， 石咀
	!  灯笼山   马骝洲   大横琴   三灶   莲花大桥  大三洲 
	 !write(121,"(10f8.3)")  el(11798), el(4582),el(10248),el(8798), el(6371),el(6162)

	
	
		
	!write(122,"(100f8.3)")   ua(22360), va(22360),  ua(20414), va(20414),  ua(16519), va(16519),  ua(13488), va(13488),   &
	!                        ua(16539), va(16539),  ua(13556), va(13556),  ua(10889), va(10889),  ua(15911), va(15911)
   !   angle1=atan2(ua(22360),va(22360))*180/3.1415926
	  !if(angle1<0)  angle1=angle1+360.
	  !angle2=atan2(ua(20414),va(20414))*180/3.1415926
	  !if(angle2<0)  angle2=angle2+360.
	  !angle3=atan2(ua(16519),va(16519))*180/3.1415926
	  !if(angle3<0)  angle3=angle3+360.
	  !angle4=atan2(ua(13488),va(13488))*180/3.1415926
	  !if(angle4<0)  angle4=angle4+360.
	  !angle5=atan2(ua(16539),va(16539))*180/3.1415926
	  !if(angle5<0)  angle5=angle5+360.
	  !angle6=atan2(ua(13556),va(13556))*180/3.1415926
	  !if(angle6<0)  angle6=angle6+360.
	  !angle7=atan2(ua(10889),va(10889))*180/3.1415926
	  !if(angle7<0)  angle7=angle7+360.
	  !angle8=atan2(ua(15911),va(15911))*180/3.1415926
	  !if(angle8<0)  angle8=angle8+360.
!	  angle9=atan2(ua(26128),va(26128))*180/3.1415926
!	  if(angle9<0)  angle9=angle9+360.

   !   write(123,"(100f8.3)")   angle1, angle2 , angle3, angle4, angle5, angle6, angle7, angle8
   !
   !   if(sediment==1)  &
	  !write(124,"(100f8.3)") sed(1)%conc(11798,1),sed(1)%conc(11431,1),sed(1)%conc(8331,1), &
	  !                       sed(1)%conc(7450,1), sed(1)%conc(9010,1), sed(1)%conc(7300,1),&
			!				 sed(1)%conc(5068,1),  sed(1)%conc(8175,1)
	                        
	 return
	 end

	 
	
