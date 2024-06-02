          subroutine  chazhi
          USE ALL_VARS
          USE MOD_OBCS
          use bcs

		 REAL(SP),dimension(19,2000) :: cw

        do i=1, ELO_TM%NTIMES
    	READ(INJUL,*) (cw(J,I),J=1,19)
        
		
		do j=1,10
	!	ELSBC(J,I)=cw(1,i)
		enddo
		
		do j=11,40
		a=float(j-11)/29. 
!		ELSBC(J,I)=cw(3,i)*a  +  cw(1,i)*(1.-a)
		enddo

		do j=1,43

		ELSBC(J,I)=cw(3,i)
		enddo



		do j=44,82
		a=float(j-44)/38.
		ELSBC(J,I)=cw(6,i)*a  +  cw(3,i)*(1.-a)
		enddo
        do j=83,124
        a= float(j-83)/41.
		ELSBC(J,I)=cw(9,i)*a  +  cw(6,i)*(1.-a)
        enddo 
        do j=125,166
		a=float(j-125)/41.
		ELSBC(J,I)=cw(12,i)*a  +  cw(9,i)*(1.-a)
        enddo
        do j=167,206 
        a=float(j-167)/39.
		ELSBC(J,I)=cw(15,i)*a  +  cw(12,i)*(1.-a)
        enddo 


		do j=207,252    
		
		a=float(j-207)/45.
	!	ELSBC(J,I)=cw(18,i)*a  +  cw(15,i)*(1.-a)
	    elsbc(j,i)=cw(15,i) 
		enddo
        do j=253,267

	!	ELSBC(J,I)=cw(18,i)
	    elsbc(j,i)=cw(15,i)

		enddo
!#########################################
        do j=1,267

       elsbc(j,i)=cw(4,i)*0.7
        enddo
!############################################

		do j=268,269

         elsbc(j,i)=cw(19,i)
		
		enddo

       
      
           enddo    ! i
     

		  end