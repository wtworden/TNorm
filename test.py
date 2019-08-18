from sage.all import *

import os

tmp = str(sage.misc.misc.SAGE_TMP)
print(tmp)
os.listdir(tmp)

S=sage.plot.plot3d.shapes.Sphere(.5)
p=S.plot()
p.show(viewer='threejs',online=True)

os.listdir(tmp)