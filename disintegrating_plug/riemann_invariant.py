import sympy
from probe import show

from acoustic import calc_amplitude_ratio
from rhd import (
    p,
    psi)


def calc_planar_riemann_invariant():

    q = calc_amplitude_ratio()
    return sympy.log(p) + q*psi

if __name__ == '__main__':

    show(locals())
