# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

from numpy import *


def sincc(x):
    th = 1e-4

    if abs(x) < th:
        val = 1.
    else:
        val = sin(x) / x

    return val


def cosincc(x):
    th = 1e-4
    if abs(x) < th:
        val = 0.5
    else:
        val = (1 - cos(x)) / (x * x)

    return val


def cubsincc(x):
    th = 1e-4
    if abs(x) < th:
        val = 1. / 6.
    else:
        val = (x - sin(x)) / (x * x * x)

    return val
