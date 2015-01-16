*----------------------------------------------------------------------*
* COMPILE USING
*  f2py -m seg_impact -c update_seg_impact.f	

	subroutine update_seg_impact(N_mp,i_seg_mp,wei_mp, N_seg, hist)
Cf2py intent(in)  N_mp   
Cf2py intent(in)  i_seg_mp                                      
Cf2py intent(in)  wei_mp
Cf2py intent(in)  N_seg
Cf2py intent(inout) hist
        implicit none
	integer  N_mp
	integer  i_seg_mp(N_mp)
	real*8   wei_mp(N_mp)
	integer  N_seg
	real*8   hist(N_seg)
	integer  p
	real*8   wei_mp_curr
	integer  i
	
        
	!write(*,*) 'ciao fortran'
	
	do p=1,N_mp
		!loop over particles
        	
		i=i_seg_mp(p)+1  	
		wei_mp_curr=wei_mp(p);
        
        	!interpolate charge to nodes
        	if (i<=N_seg .and. i>=1) then
        	hist(i) = hist(i) + wei_mp_curr;
		end if

        	
        
	end do
        


	end subroutine
	
	
