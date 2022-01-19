import sympy
from probe import show

def calc_hydrostatic_density():

    sympy.var('M_e R_e x omega', positive=True)

    return (M_e/R_e**3)*(x/R_e)**omega

if __name__ == '__main__':

    show(locals())
