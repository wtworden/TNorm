name = "tnorm"

load_attach_path().remove('/Applications/Regina.app/Contents/MacOS/python')

from tnorm.version import __version__

import tnorm.kernel
import tnorm.GUI

from tnorm.TN_wrapper import TN_wrapper as load

