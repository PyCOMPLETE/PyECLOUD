#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 6.1.1
#
#
#     Author and contact:   Giovanni IADAROLA
#                           BE-ABP Group
#                           CERN
#                           CH-1211 GENEVA 23
#                           SWITZERLAND
#                           giovanni.iadarola@cern.ch
#
#                contact:   Giovanni RUMOLO
#                           BE-ABP Group
#                           CERN
#                           CH-1211 GENEVA 23
#                           SWITZERLAND
#                           giovanni.rumolo@cern.ch
#
#
#
#     Copyright  CERN,  Geneva  2011  -  Copyright  and  any   other
#     appropriate  legal  protection  of  this  computer program and
#     associated documentation reserved  in  all  countries  of  the
#     world.
#
#     Organizations collaborating with CERN may receive this program
#     and documentation freely and without charge.
#
#     CERN undertakes no obligation  for  the  maintenance  of  this
#     program,  nor responsibility for its correctness,  and accepts
#     no liability whatsoever resulting from its use.
#
#     Program  and documentation are provided solely for the use  of
#     the organization to which they are distributed.
#
#     This program  may  not  be  copied  or  otherwise  distributed
#     without  permission. This message must be retained on this and
#     any other authorized copies.
#
#     The material cannot be sold. CERN should be  given  credit  in
#     all references.
#----------------------------------------------------------------------

import dynamics_dipole as dyndip
import dynamics_Boris as dynB
import dynamics_strong_B_generalized as dyngen
from numpy import array

Dt=25e-12;
N_steps=100
B=0.1
N_sub_steps=100

dynamicsd=dyndip.pusher_dipole_magnet(Dt,B)
dynamicsB=dynB.pusher_Boris(Dt, 0., B, 0., \
                 None, None, None,N_sub_steps=N_sub_steps)
dynamicsGen=dyngen.pusher_strong_B_generalized(Dt, 0., B, \
                 None, None)

N_mp=2
Ex_n=array([0.001,0.001])
Ey_n=array([1.,0.])
x_mpB=array([0.,0.])
y_mpB=array([0.,0.])
z_mpB=array([0.,0.])

x_mpd=x_mpB.copy()
y_mpd=y_mpB.copy()
z_mpd=z_mpB.copy()

x_mpg=x_mpB.copy()
y_mpg=y_mpB.copy()
z_mpg=z_mpB.copy()

vx_mpB=array([1.,0.])
vy_mpB=array([0.,0.])
vz_mpB=array([0.,1.])

vx_mpd=vx_mpB.copy()
vy_mpd=vy_mpB.copy()
vz_mpd=vz_mpB.copy()

vx_mpg=vx_mpB.copy()
vy_mpg=vy_mpB.copy()
vz_mpg=vz_mpB.copy()


x_lisB=[];
y_lisB=[];
z_lisB=[];

x_lisd=[];
y_lisd=[];
z_lisd=[];

x_lisg=[];
y_lisg=[];
z_lisg=[];

vx_lisB=[];
vy_lisB=[];
vz_lisB=[];

vx_lisd=[];
vy_lisd=[];
vz_lisd=[];

vx_lisg=[];
vy_lisg=[];
vz_lisg=[];


for ii in range(N_steps):
    #print ii
    x_lisB.append(x_mpB.copy())
    y_lisB.append(y_mpB.copy())
    z_lisB.append(z_mpB.copy())

    x_lisd.append(x_mpd.copy())
    y_lisd.append(y_mpd.copy())
    z_lisd.append(z_mpd.copy())

    x_lisg.append(x_mpg.copy())
    y_lisg.append(y_mpg.copy())
    z_lisg.append(z_mpg.copy())

    vx_lisB.append(vx_mpB.copy())
    vy_lisB.append(vy_mpB.copy())
    vz_lisB.append(vz_mpB.copy())

    vx_lisd.append(vx_mpd.copy())
    vy_lisd.append(vy_mpd.copy())
    vz_lisd.append(vz_mpd.copy())

    vx_lisg.append(vx_mpg.copy())
    vy_lisg.append(vy_mpg.copy())
    vz_lisg.append(vz_mpg.copy())


    x_mpB[0:N_mp], y_mpB[0:N_mp], z_mpB[0:N_mp], vx_mpB[0:N_mp], vy_mpB[0:N_mp], vz_mpB[0:N_mp]=\
                   dynamicsB.step(x_mpB[0:N_mp], y_mpB[0:N_mp], z_mpB[0:N_mp],\
                   vx_mpB[0:N_mp], vy_mpB[0:N_mp],vz_mpB[0:N_mp],Ex_n[0:N_mp],Ey_n[0:N_mp]);

    x_mpd[0:N_mp], y_mpd[0:N_mp], z_mpd[0:N_mp], vx_mpd[0:N_mp], vy_mpd[0:N_mp], vz_mpd[0:N_mp]=\
                   dynamicsd.step(x_mpd[0:N_mp], y_mpd[0:N_mp], z_mpd[0:N_mp],\
                   vx_mpd[0:N_mp], vy_mpd[0:N_mp],vz_mpd[0:N_mp],Ex_n[0:N_mp],Ey_n[0:N_mp]);

    x_mpg[0:N_mp], y_mpg[0:N_mp], z_mpg[0:N_mp], vx_mpg[0:N_mp], vy_mpg[0:N_mp], vz_mpg[0:N_mp]=\
                   dynamicsGen.step(x_mpg[0:N_mp], y_mpg[0:N_mp], z_mpg[0:N_mp],\
                   vx_mpg[0:N_mp], vy_mpg[0:N_mp],vz_mpg[0:N_mp],Ex_n[0:N_mp],Ey_n[0:N_mp]);



x_lisB=array(x_lisB)
y_lisB=array(y_lisB)
z_lisB=array(z_lisB)

x_lisd=array(x_lisd)
y_lisd=array(y_lisd)
z_lisd=array(z_lisd)

x_lisg=array(x_lisg)
y_lisg=array(y_lisg)
z_lisg=array(z_lisg)

vx_lisB=array(vx_lisB)
vy_lisB=array(vy_lisB)
vz_lisB=array(vz_lisB)

vx_lisd=array(vx_lisd)
vy_lisd=array(vy_lisd)
vz_lisd=array(vz_lisd)

vx_lisg=array(vx_lisg)
vy_lisg=array(vy_lisg)
vz_lisg=array(vz_lisg)

import pylab as pl
pl.close('all')

for ii in range(len(x_lisB[1])):
    t_vect = Dt*pl.arange(len(x_lisB[:,ii]))*1e9

    pl.figure(2*ii)
    fst=pl.subplot(3,1,1)
    pl.plot(t_vect, x_lisd[:,ii])
    pl.hold('on')
    pl.plot(t_vect, x_lisB[:,ii],'.r')
    pl.plot(t_vect, x_lisg[:,ii],'xg')

    pl.subplot(3,1,2, sharex=fst)
    pl.plot(t_vect, y_lisd[:,ii])
    pl.hold('on')
    pl.plot(t_vect, y_lisB[:,ii],'.r')
    pl.plot(t_vect, y_lisg[:,ii],'xg')


    pl.subplot(3,1,3, sharex=fst)
    pl.plot(t_vect, z_lisd[:,ii])
    pl.hold('on')
    pl.plot(t_vect, z_lisB[:,ii],'.r')
    pl.plot(t_vect, z_lisg[:,ii],'xg')
    pl.xlabel('Time [ns]')

    pl.figure(2*ii+1)
    fst=pl.subplot(3,1,1)
    pl.plot(t_vect, vx_lisd[:,ii])
    pl.hold('on')
    pl.plot(t_vect, vx_lisB[:,ii],'.r')
    pl.plot(t_vect, vx_lisg[:,ii],'xg')

    pl.subplot(3,1,2, sharex=fst)
    pl.plot(t_vect, vy_lisd[:,ii])
    pl.hold('on')
    pl.plot(t_vect, vy_lisB[:,ii],'.r')
    pl.plot(t_vect, vy_lisg[:,ii],'xg')


    pl.subplot(3,1,3, sharex=fst)
    pl.plot(t_vect, z_lisd[:,ii])
    pl.hold('on')
    pl.plot(t_vect, z_lisB[:,ii],'.r')
    pl.plot(t_vect, z_lisg[:,ii],'xg')
    pl.xlabel('Time [ns]')

pl.show()
