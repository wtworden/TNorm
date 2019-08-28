
from tnorm.kernel.matrices import sign_matrix

QUIET = None
SIGN_MATRIX = sign_matrix()
ABS_SIGN_MATRIX = sign_matrix().apply_map(abs)