from tnorm.utilities.temp_dir import make_temp_directory
from tnorm.utilities.sage_types import save,load

from zipfile import ZipFile
import os

import regina
import snappy


def save_file(TN_wrapper, filename):
    W = TN_wrapper
    with make_temp_directory() as temp_dir:
        T = W._triangulation
        T_file = os.path.join(temp_dir, 'T_file.rga')
        T.save(T_file)
        W._manifold = None
        W._triangulation = None
        W_file = os.path.join(temp_dir, 'W_file')
        save(W,W_file)
        with ZipFile(filename, 'w') as myzip:
            myzip.write(T_file, arcname='T_file.rga')
            myzip.write(W_file+'.sobj', arcname='W_file.sobj')

def open_file(filename):
    with make_temp_directory() as temp_dir:
        with ZipFile(filename) as myzip:
            myzip.extract('W_file.sobj', path=temp_dir)
            W_file = os.path.join(temp_dir, 'W_file.sobj')
            W = load(W_file)
            myzip.extract('T_file.rga', path=temp_dir)
            T_file = os.path.join(temp_dir, 'T_file.rga')
            T = regina.open(T_file)
    W._triangulation = T
    W._manifold = snappy.Manifold(T.snapPea())
    return W

