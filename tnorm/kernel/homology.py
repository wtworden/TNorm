
from tnorm.utilities.sage_types import vector


def _map_to_H1bdy(oriented_spun_surface, TN_wrapper):
    s = oriented_spun_surface
    pos_slopes, neg_slopes = TN_wrapper.boundary_slopes(s).values()

    slopes = []
    for i in range(len(pos_slopes)):
        slope = (pos_slopes[i][0]+neg_slopes[i][0], pos_slopes[i][1]+neg_slopes[i][1])
        slopes.append(slope)

    return slopes


### this maps an oriented spun normal surface to H_2(M,\partial M; \ZZ)
def _map_to_H2(oriented_spun_surface, TN_wrapper):
    s = oriented_spun_surface
    slopes = TN_wrapper.H1bdy_slopes(s)
    coords = vector([slopes[i][1] for i in range(len(slopes))])
    return coords


#################--------------------------------------------------------------------



