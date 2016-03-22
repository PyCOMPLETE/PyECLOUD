# secondary emission model
import sys,os
BIN = os.path.expanduser('../../../')
sys.path.append(BIN)
BIN = os.path.expanduser('../../')
sys.path.append(BIN)
BIN = os.path.expanduser('../')
sys.path.append(BIN)
from PyECLOUD.buildup_simulation import BuildupSimulation
print 'Imported from PyECLOUD folder'

import datetime
Start=datetime.datetime.now()
sim = BuildupSimulation()
sim.run()
print datetime.datetime.now()-Start

