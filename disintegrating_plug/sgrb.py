import sympy
from probe import show
from rhd import (r)

sympy.var('h_e M_e', positive=True)

def make_density_profile():

    return (M_e/r**3)*sympy.exp(-r**2/h_e**2)

if __name__ == '__main__':

    show(locals())

