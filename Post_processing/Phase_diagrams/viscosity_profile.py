import numpy as np
import pandas as pd


mineral_data = pd.read_csv('mineral_dislocation_creep.ods')
###########################
## load mineral volume fractions
olivine=np.loadtxt('Olivine.dat')
garnet=np.loadtxt('Garnet.dat')
