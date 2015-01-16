from cmath import *
from numpy import *
from errffor import errf

def wfun(z):
    x=z.real
    y=z.imag
    wx,wy=errf(x,y)
    return wx+1j*wy

def BassErsk(xin,yin,sigmax,sigmay):
        
    x=abs(xin);
    y=abs(yin);
    
    eps0=8.854187817620e-12;
    
    
    if sigmax>sigmay:
    
        S=sqrt(2*(sigmax*sigmax-sigmay*sigmay));
        factBE=1/(2*eps0*sqrt(pi)*S);
        etaBE=sigmay/sigmax*x+1j*sigmax/sigmay*y;
        zetaBE=x+1j*y;
        
        val=factBE*(wfun(zetaBE/S)-exp( -x*x/(2*sigmax*sigmax)-y*y/(2*sigmay*sigmay))*wfun(etaBE/S) );
           
        Ex=abs(val.imag)*sign(xin);
        Ey=abs(val.real)*sign(yin);
    
    else:
    
        S=sqrt(2*(sigmay*sigmay-sigmax*sigmax));
        factBE=1/(2*eps0*sqrt(pi)*S);
        etaBE=sigmax/sigmay*y+1j*sigmay/sigmax*x;
        yetaBE=y+1j*x;
        
        val=factBE*(wfun(yetaBE/S)-exp( -y*y/(2*sigmay*sigmay)-x*x/(2*sigmax*sigmax))*wfun(etaBE/S) );
           
        Ey=abs(val.imag)*sign(yin);
        Ex=abs(val.real)*sign(xin);
         
    return Ex, Ey

def ImageTerms(x,y,a,b,x0,y0, nimag):
    
        
    eps0=8.854187817620e-12;    
    
    if abs((a-b)/a)>1e-3:    
        g=sqrt(a*a-b*b)
        
        z=x+1j*y
        q=acosh(z/g)
        mu=q.real
        phi=q.imag
        
        z0=x0+1j*y0
        q0=acosh(z0/g)
        mu0=q0.real
        phi0=q0.imag
        
        mu1=0.5*log((a+b)/(a-b))
        
        Ecpx=0+0j
        
        
        q=conj(q)
        for nn in range(1,nimag+1):
            Ecpx=Ecpx+exp(-nn*mu1) * ( (cosh(nn*mu0)*cos(nn*phi0)) / (cosh(nn*mu1)) + 1j * (sinh(nn*mu0)*sin(nn*phi0)) / (sinh(nn*mu1))   )* (sinh(nn*q))/(sinh(q))
            
        
        Ecpx=Ecpx/(4*pi*eps0)*4/g
        Ex=Ecpx.real
        Ey=Ecpx.imag
    else:
        if (x0==0) and (y0==0):
            Ex=0.
            Ey=0.
        else:
            print('This case has not been implemented yet')
    
    return Ex, Ey
    
#    eps0=8.854187817620e-12;
#    
#    
#    if sigmax>sigmay:
#    
#        S=sqrt(2*(sigmax*sigmax-sigmay*sigmay));
#        factBE=1/(2*eps0*sqrt(pi)*S);
#        etaBE=sigmay/sigmax*x+1j*sigmax/sigmay*y;
#        zetaBE=x+1j*y;
#        
#        val=factBE*(wfun(zetaBE/S)-exp( -x*x/(2*sigmax*sigmax)-y*y/(2*sigmay*sigmay))*wfun(etaBE/S) );
#           
#        Ex=abs(val.imag)*sign(xin);
#        Ey=abs(val.real)*sign(yin);
#    
#    else:
#    
#        S=sqrt(2*(sigmay*sigmay-sigmax*sigmax));
#        factBE=1/(2*eps0*sqrt(pi)*S);
#        etaBE=sigmax/sigmay*y+1j*sigmay/sigmax*x;
#        yetaBE=y+1j*x;
#        
#        val=factBE*(wfun(yetaBE/S)-exp( -y*y/(2*sigmay*sigmay)-x*x/(2*sigmax*sigmax))*wfun(etaBE/S) );
#           
#        Ey=abs(val.imag)*sign(yin);
#        Ex=abs(val.real)*sign(xin);
#    
#    
#    
#    return Ex, Ey
