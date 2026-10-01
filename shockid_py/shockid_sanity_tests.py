from shockid import getWaveSpeeds
import numpy as np 


ro, pr = np.array([1.0]), np.array([0.6])       # cs^2 = 1 for gamma=5/3
# perpendicular propagation: Bn=0 -> slow=0, fast=cs2+va2
s = getWaveSpeeds(ro, pr, np.array([0.0]), np.array([2.0]), 0.0, np.array([0.0]))
print(s['vslow2'], s['vfast2'], 'should be [0.], [5.]')                 # [0.], [5.]
# parallel propagation, va2=4 > cs2=1 -> slow=cs2, fast=va2
s = getWaveSpeeds(ro, pr, np.array([2.0]), np.array([0.0]), 0.0, np.array([2.0]))
print(s['vslow2'], s['vfast2'], 'should be [1.], [4.]')                 # [1.], [4.]
# hydro limit
s = getWaveSpeeds(ro, pr, np.array([0.0]), np.array([0.0]), 0.0, np.array([0.0]))
print(s['vslow2'], s['vfast2'], 'should be [0.], [1.]')                 # [0.], [1.]
