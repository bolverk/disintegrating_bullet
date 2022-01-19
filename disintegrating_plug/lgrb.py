import sympy
from probe import show
from equation_of_motion import alpha
from rhd import c

sympy.var('x rho_a', positive=True)

def calc_hydrostatic_density():

    sympy.var('M_e R_e omega', positive=True)

    return (M_e/R_e**3)*(x/R_e)**omega

def eruption_criterion():

    sympy.var('L')

    return L/(alpha**2*x**2*c) - rho_a*c**2/alpha**4

def calc_eruption_depth():

    _ = eruption_criterion()
    _ = _.subs(rho_a, calc_hydrostatic_density())
    return sympy.solve(_, x)[0]

if __name__ == '__main__':

    show(locals())
