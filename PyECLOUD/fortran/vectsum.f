*----------------------------------------------------------------------*
* COMPILE USING
*  f2py -m vectsum -c vectsum.f       

        subroutine vectsum(N_mp, x, res)
Cf2py intent(in) N_mp
Cf2py intent(in) x
Cf2py intent(out) res
        implicit none
        integer N_mp
        integer p
        real*8  x(N_mp)
        real*8  res

        res = 0

        do p=1,N_mp        
        res = res+x(p)        
        end do
             
        end subroutine

