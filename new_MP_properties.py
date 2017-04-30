from __future__ import division
import numpy as np
import numpy.random as random
from scipy.constants import e as qe, m_e

qm = qe/m_e

def velocities_angle_cosine(N_new_MP, En_gen, Norm_x, Norm_y):
    v_gen_mod=np.sqrt(2.*qm*En_gen)

    phi_p = random.rand(N_new_MP)*2*np.pi
    sin_phi_p = np.sin(phi_p)
    cos_phi_p = np.cos(phi_p)

    sin_theta_sq = random.rand(N_new_MP)
    cos_theta_p = np.sqrt(1-sin_theta_sq)
    sin_theta_p = np.sqrt(sin_theta_sq)
    vx_gen = v_gen_mod*\
        (cos_theta_p*Norm_x+sin_theta_p*sin_phi_p*Norm_y)
    vy_gen = v_gen_mod*\
        (cos_theta_p*Norm_y-sin_theta_p*sin_phi_p*Norm_x)
    vz_gen = v_gen_mod*(sin_theta_p*cos_phi_p)

    return vx_gen, vy_gen, vz_gen

def velocities_angle_cosine_old(N_new_MP, En_gen, Norm_x, Norm_y):
    v_gen_mod=np.sqrt(2.*qm*En_gen)

    phi_p = random.rand(N_new_MP)*2*np.pi
    sin_phi_p = np.sin(phi_p)
    cos_phi_p = np.cos(phi_p)

    sin_theta_p = random.rand(N_new_MP)
    cos_theta_p = np.sqrt(1-sin_theta_p**2)
    vx_gen = v_gen_mod*\
        (cos_theta_p*Norm_x+sin_theta_p*sin_phi_p*Norm_y)
    vy_gen = v_gen_mod*\
        (cos_theta_p*Norm_y-sin_theta_p*sin_phi_p*Norm_x)
    vz_gen = v_gen_mod*(sin_theta_p*cos_phi_p)

    return vx_gen, vy_gen, vz_gen

