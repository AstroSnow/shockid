# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pipreadmods as PIPpy
import matplotlib.pyplot as plt
import numpy as np
from shockid import shockid


fname='/home/bs428/Downloads/thermal_tearing2d/adiabatic'

ds=PIPpy.pipread('/home/bs428/Documents/simdata/MHD_ref/',10)

#Convert the data into a 2D array
ro=np.tile(ds['ro_p'],(20,1))
pr=np.tile(ds['pr_p'],(20,1))
vx=np.tile(ds['vx_p'],(20,1))
vy=np.tile(ds['vy_p'],(20,1))
vz=np.tile(ds['vz_p'],(20,1))
bx=np.tile(ds['bx'],(20,1))
by=np.tile(ds['by'],(20,1))
bz=np.tile(ds['bz'],(20,1))

#plt.contourf(divv)
plt.plot(ro[1,:])