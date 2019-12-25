name = "tnorm"

import sys
try:
    sys.path.remove('/Applications/Regina.app/Contents/MacOS/python')
except ValueError:
    pass

from tnorm.version import __version__
version = __version__

import tnorm.kernel
import tnorm.GUI

from tnorm.TN_wrapper import TN_wrapper as load

#IN_SAGE=False