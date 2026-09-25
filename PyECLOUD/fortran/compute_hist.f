*----------------------------------------------------------------------*
* COMPILE USING
*  f2py -m rhocompute -c compute_rho.f	

	subroutine compute_hist(N_mp,x_mp,wei_mp, bias_x, 
     +	Dx, Nxg, hist)
Cf2py intent(in)  N_mp   
Cf2py intent(in)  x_mp                                      
Cf2py intent(in)  wei_mp
Cf2py intent(in)  bias_x
Cf2py intent(in)  Dx
Cf2py intent(in)  Nxg
Cf2py intent(inout) hist
        implicit none
	integer  N_mp
	real*8   x_mp(N_mp), wei_mp(N_mp)
	real*8   bias_x, Dx
	integer  Nxg
	real*8   hist(Nxg)
	integer  p
	real*8   wei_mp_curr, fi, hx
	integer  i
	
        
	!write(*,*) 'ciao fortran'
	
	do p=1,N_mp
		!loop over particles
        	
        	fi = 1+(x_mp(p)-bias_x)/Dx;   !real i index of particle's cell 
        	i = int(fi);                  !integral part
        	hx = fi-dble(i);              !the remainder
        
        	
         	wei_mp_curr=wei_mp(p);
        
        	!interpolate charge to nodes
		!write (*,*) 'Nxg=', Nxg
		!write (*,*) 'i=', i

		if (i<Nxg) then
        	hist(i) = hist(i) + wei_mp_curr*(1-hx);
        	hist(i+1) = hist(i+1) + wei_mp_curr*hx;
        	else
        	hist(Nxg)=hist(Nxg)+ wei_mp_curr
		end if
       
	end do
        


	end subroutine
	
	
