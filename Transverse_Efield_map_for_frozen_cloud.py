import numpy as np
from scipy.constants import c

from PyHEADTAIL.particles.slicing import UniformBinSlicer
from PyPIC.PyPIC_Scatter_Gather import PyPIC_Scatter_Gather


class Transverse_Efield_map(object):
    def __init__(self, xg, yg, Ex, Ey, L_interaction, slicer,
                 flag_clean_slices=False, wrt_slice_centroid=False,
                 x_beam_offset=0., y_beam_offset=0., slice_by_slice_mode=False):

        self.slicer = slicer
        self.L_interaction = L_interaction
        self.flag_clean_slices = flag_clean_slices
        self.wrt_slice_centroid = wrt_slice_centroid

        self.Ex = Ex
        self.Ey = Ey
        self.pic = PyPIC_Scatter_Gather(xg=xg, yg=yg)

        self.x_beam_offset = x_beam_offset
        self.y_beam_offset = y_beam_offset

        self.slice_by_slice_mode = slice_by_slice_mode
        if self.slice_by_slice_mode:
            self.sid = 0
            self.track = self._track_in_single_slice_mode
            self.finalize_and_reinitialize = self._finalize_and_reinitialize
        else:
            self.finalize_and_reinitialize = None

    def get_beam_x(self, beam):
        return beam.x

    def get_beam_y(self, beam):
        return beam.y

    def track(self, beam):
        if self.flag_clean_slices:
            beam.clean_slices()

        slices = beam.get_slices(self.slicer)

        sid = 0
        for _ in range(slices.n_slices):

            sid -= 1

            # select particles in the slice
            pid = slices.particle_indices_of_slice(sid)

            # slice size
            dz = (slices.z_bins[sid + 1] - slices.z_bins[sid])

            self._track_single_slice(beam, sid, pid)

    def _track_single_slice(self, beam, sid, pid):
        x = self.get_beam_x(beam)[pid]
        y = self.get_beam_y(beam)[pid]

        self.pic.efx = np.squeeze(self.Ex[sid, :, :])
        self.pic.efy = np.squeeze(self.Ey[sid, :, :])

        centroid_x = 0
        centroid_y = 0
        if self.wrt_slice_centroid:
            centroid_x = np.mean(x)
            centroid_y = np.mean(y)

        Ex_sc_p, Ey_sc_p = self.pic.gather(
            x - centroid_x + self.x_beam_offset,
            y - centroid_y + self.y_beam_offset)

        # kick beam particles
        fact_kick = beam.charge / (beam.p0 * beam.beta * c) * self.L_interaction
        beam.xp[pid] += fact_kick * Ex_sc_p
        beam.yp[pid] += fact_kick * Ey_sc_p

    def _track_in_single_slice_mode(self, beam):

        if beam.slice_info is not 'unsliced':
            self.sid -= 1
            self._track_single_slice(beam=beam, sid=self.sid, pid=np.arange(beam.macroparticlenumber))

    def _finalize_and_reinitialize(self):
        self.sid = 0
