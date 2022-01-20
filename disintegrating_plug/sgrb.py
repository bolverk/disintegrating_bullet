import sympy
from probe import show
from rhd import (r, c)
from equation_of_motion import (alpha, M_1, gamma_i, t_i)

sympy.var('h_e M_e L', positive=True)

def make_density_profile():

    return (M_e/r**3)*sympy.exp(-r**2/h_e**2)

def make_eruption_radius():

    return h_e*sympy.sqrt(sympy.log(M_e*c**3/L/alpha**2/h_e))

def make_eruption_density():

    return alpha**2*L/(h_e**2*c**3)

def make_eruption_mass():

    return alpha**4*L*h_e/c**3

if __name__ == '__main__':

    show(locals())

