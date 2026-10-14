# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

#for compatibility with old PyECLOUD versions

try:
    from .buildup_simulation import BuildupSimulation
    print('Imported from local folder')
except ImportError:
    from PyECLOUD.buildup_simulation import BuildupSimulation
    print('Imported from PyECLOUD folder')


sim = BuildupSimulation()
sim.run()
