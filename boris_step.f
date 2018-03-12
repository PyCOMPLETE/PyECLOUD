*----------------------------------------------------------------------*
* COMPILE USING
*  f2py -m boris_step -c boris_step.f        

        subroutine boris_step(N_mp, Dtt, xn1, yn1,  zn1, 
     +                                vxn1, vyn1,  vzn1,
     +                                Ex_n, Ey_n, Ez_n,
     +                                Bx_n, By_n, Bz_n,
     +                                mass, charge)
Cf2py intent(in) N_mp
Cf2py intent(in) Dtt
Cf2py intent(inout) xn1
Cf2py intent(inout) yn1
Cf2py intent(inout) zn1
 
Cf2py intent(inout) vxn1
Cf2py intent(inout) vyn1
Cf2py intent(inout) vzn1

Cf2py intent(in) Ex_n 
Cf2py intent(in) Ey_n
Cf2py intent(in) Ez_n

Cf2py intent(in) Bx_n
Cf2py intent(in) By_n
Cf2py intent(in) Bz_n

Cf2py intent(in) mass
Cf2py intent(in) charge


        implicit none
        integer  N_mp
        real*8   Dtt
        real*8   xn1(N_mp), yn1(N_mp),  zn1(N_mp)
        real*8   vxn1(N_mp), vyn1(N_mp),  vzn1(N_mp)
        real*8   Ex_n(N_mp), Ey_n(N_mp),  Ez_n(N_mp)
        real*8   Bx_n(N_mp), By_n(N_mp),  Bz_n(N_mp)

        real*8   mass, charge, qm
        integer  p
        real*8   tBx, tBy, tBz, tBsq
        real*8   sBx, sBy, sBz
        real*8   vx_prime, vy_prime, vz_prime
        real*8   vx_min, vy_min, vz_min
        real*8   vx_plus, vy_plus, vz_plus
        
        real*8  Ex_np, Ey_np, Ez_np 
        real*8  vxn1p, vyn1p, vzn1p
        
        qm=charge/mass
        
        do p=1,N_mp
        
        Ex_np = Ex_n(p)
        Ey_np = Ey_n(p)
        Ez_np = Ez_n(p)
        
        vxn1p = vxn1(p)
        vyn1p = vyn1(p)
        vzn1p = vzn1(p)

        tBx = 0.5*qm*Dtt*Bx_n(p)
        tBy = 0.5*qm*Dtt*By_n(p)
        tBz = 0.5*qm*Dtt*Bz_n(p)

        tBsq = tBx*tBx + tBy*tBy + tBz*tBz
        
        sBx = 2.*tBx/(1.+tBsq)
        sBy = 2.*tBy/(1.+tBsq)
        sBz = 2.*tBz/(1.+tBsq)
        
        vx_min = vxn1p + 0.5*qm*Ex_np*Dtt
        vy_min = vyn1p + 0.5*qm*Ey_np*Dtt
        vz_min = vzn1p + 0.5*qm*Ez_np*Dtt
        
        !v_prime = v_min + cross(v_min, tB)
        vx_prime = vy_min*tBz-vz_min*tBy + vx_min 
        vy_prime = vz_min*tBx-vx_min*tBz + vy_min
        vz_prime = vx_min*tBy-vy_min*tBx + vz_min                
        
        !v_plus = v_min + cross(v_prime, sB)
        vx_plus = vy_prime*sBz-vz_prime*sBy + vx_min 
        vy_plus = vz_prime*sBx-vx_prime*sBz + vy_min
        vz_plus = vx_prime*sBy-vy_prime*sBx + vz_min
        
        vxn1p = vx_plus + 0.5*qm*Ex_np*Dtt
        vyn1p = vy_plus + 0.5*qm*Ey_np*Dtt
        vzn1p = vz_plus + 0.5*qm*Ez_np*Dtt
        
        xn1(p) = xn1(p) + vxn1p * Dtt
        yn1(p) = yn1(p) + vyn1p * Dtt
        zn1(p) = zn1(p) + vzn1p * Dtt
        
        vxn1(p) = vxn1p
        vyn1(p) = vyn1p
        vzn1(p) = vzn1p
        
        
        end do
             
        end subroutine

