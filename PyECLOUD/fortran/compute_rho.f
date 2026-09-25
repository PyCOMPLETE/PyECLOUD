*----------------------------------------------------------------------*
* COMPILE USING
*  f2py -m rhocompute -c compute_rho.f	

	subroutine compute_sc_rho(N_mp,x_mp, y_mp,nel_mp, bias_x, 
     +	bias_y, Dh, Nxg, Nyg,rho)
Cf2py intent(in)  N_mp   
Cf2py intent(in)  x_mp                                      
Cf2py intent(in)  y_mp   
Cf2py intent(in)  nel_mp
Cf2py intent(in)  bias_x
Cf2py intent(in)  bias_y 
Cf2py intent(in)  Dh
Cf2py intent(in)  Nxg
Cf2py intent(in)  Nyg
Cf2py intent(out) rho
        implicit none
	integer  N_mp
	real*8   x_mp(N_mp), y_mp(N_mp), nel_mp(N_mp)
	real*8   bias_x, bias_y, Dh
	integer  Nxg, Nyg
	real*8   rho(Nxg, Nyg)
	integer  p
	real*8   nel_mp_curr, fi, fj, hx, hy
	integer  i, j
	
        !rho=0d0;
	do j=1,Nyg
		do i=1,Nxg
		rho(i,j)=0.0
                end do
	end do
	!write(*,*) 'ciao fortran'
	
	do p=1,N_mp
		!loop over particles
        	
        	fi = 1+(x_mp(p)-bias_x)/Dh;   !real i index of particle's cell 
        	i = int(fi);                !integral part
        	hx = fi-dble(i);                    !the remainder
        
        	
        	fj = 1+(y_mp(p)-bias_y)/Dh;   !real j index of particle's cell (C-like indexing)
        	j = int(fj);                !integral part
        	hy = fj-dble(j);                    !the remainder
        
        	nel_mp_curr=nel_mp(p);
            
            if (i>0 .and. j>0 .and. i<Nxg .and. j<Nyg) then
        
        	!interpolate charge to nodes
        	rho(i,j) = rho(i,j) + nel_mp_curr*(1-hx)*(1-hy);
        	rho(i+1,j) = rho(i+1,j) + nel_mp_curr*hx*(1-hy);
        	rho(i,j+1) = rho(i,j+1) + nel_mp_curr*(1-hx)*hy;
        	rho(i+1,j+1) = rho(i+1,j+1) + nel_mp_curr*hx*hy;
            end if 
	
	end do
        


	end subroutine
	
	
