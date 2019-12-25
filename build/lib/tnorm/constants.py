
from tnorm.kernel.matrices import sign_matrix

QUIET = None
SIGN_MATRIX = sign_matrix()
ABS_SIGN_MATRIX = sign_matrix().apply_map(abs)

##### for sympy implementation
#ABS_SIGN_MATRIX = abs(sign_matrix())

#IN_SAGE = False

#try:
#    from sage.all import RR
#    IN_SAGE = True
#except:
#    IN_SAGE = False