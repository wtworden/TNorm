
import contextlib
import tempfile
import shutil
from math import gcd as GCD

class cached_property(object):
    """
    Descriptor (non-data) for building an attribute on-demand on first use.
    """
    def __init__(self, factory):
        """
        <factory> is called such: factory(instance) to build the attribute.
        """
        self._attr_name = factory.__name__
        self._factory = factory

    def __get__(self, instance, owner):
        # Build the attribute.
        attr = self._factory(instance)

        # Cache the value; hide ourselves.
        setattr(instance, self._attr_name, attr)

        return attr


@contextlib.contextmanager
def make_temp_directory():
    temp_dir = tempfile.mkdtemp()
    try:
        yield temp_dir
    finally:
        shutil.rmtree(temp_dir)


def lcm(list_of_ints):
    L = list_of_ints
    LCM = 1
    for k in L:
        LCM = LCM*k//gcd(LCM,k)
    return LCM

def gcd(a,b):
    return abs(GCD(a,b))







