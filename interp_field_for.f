*----------------------------------------------------------------------*
* COMPILE USING
*  f2py -m int_field_for -c interp_field_for.f	

        subroutine int_field(N_mp,xn,yn, bias_x,bias_y, dx,dy,efx, efy,
     +	Nxg, Nyg, Ex_n, Ey_n)
Cf2py intent(in)  N_mp   
Cf2py intent(in)  xn                                      
Cf2py intent(in)  yn   
Cf2py intent(in)  bias_x
Cf2py intent(in)  bias_y 
Cf2py intent(in)  dx
Cf2py intent(in)  dy
Cf2py intent(in)  efx
Cf2py intent(in)  efy
Cf2py intent(in)  Nxg
Cf2py intent(in)  Nyg
Cf2py intent(out) Ex_n
Cf2py intent(out) Ey_n


        implicit none
	integer  N_mp
	real*8   xn(N_mp), yn(N_mp)
	real*8   bias_x, bias_y, dx, dy
	integer  Nxg,Nyg 
	real*8   efx(Nxg, Nyg), efy(Nxg, Nyg)
	integer  p
	real*8   fi, fj, hx, hy
	integer  i, j
        real*8   Ex_n(N_mp), Ey_n(N_mp)
        
        
        
        do p=1,N_mp
        fi = 1+(xn(p)-bias_x)/dx;             !i index of particle's cell 
     	i  = int(fi);
     	hx = fi-dble(i);                      !fractional x position in cell
    

     	fj = 1+(yn(p)-bias_y)/dy;             !j index of particle' cell(C-like!!!!)
     	j = int(fj);
     	hy = fj-dble(j);                      !fractional y position in cell

     	
     	!gather electric field
        if (i>0 .and. j>0 .and. i<Nxg .and. j<Nyg) then
            Ex_n(p)=efx((i),(j))*(1-hx)*(1-hy);   
            Ex_n(p) = Ex_n(p) + efx((i+1),(j))*hx*(1-hy);
            Ex_n(p) = Ex_n(p) + efx((i),(j+1))*(1-hx)*hy;
            Ex_n(p) = Ex_n(p) + efx((i+1),(j+1))*hx*hy;
        
            Ey_n(p)=efy((i),(j))*(1-hx)*(1-hy);   
            Ey_n(p) = Ey_n(p) + efy((i+1),(j))*hx*(1-hy);
            Ey_n(p) = Ey_n(p) + efy((i),(j+1))*(1-hx)*hy;
            Ey_n(p) = Ey_n(p) + efy((i+1),(j+1))*hx*hy;
        end if
     	end do
     	
     	end subroutine

