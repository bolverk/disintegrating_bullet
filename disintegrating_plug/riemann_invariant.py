import sympy
from probe import show
from caching import luggage
from my_logging import logger

from acoustic import (
    calc_amplitude_ratio,
    calc_sound_speed)
from rhd import (
    eta,
    c,
    p,
    t,
    r,
    psi,
    calc_time_derivs)


def calc_planar_riemann_invariant():

    q = calc_amplitude_ratio()
    return sympy.log(p) + q*psi

@luggage.memory
def retrieve_planar_time_derivatives():

    return calc_time_derivs(0)

@luggage.memory
def retrieve_spherical_time_derivatives():

    return calc_time_derivs(2)

def calc_spherical_riemann_invariant():

    beta = sympy.Symbol('beta', positive=True)
    xi = sympy.Symbol('xi', positive=True)

    logger.debug('begin planar riemann invariant calculation')
    j = calc_planar_riemann_invariant()
    logger.debug('end planar riemann invariant calculation')
    def rel_vel_add(v1, v2):

        return (v1+v2)/(1+v1*v2/c**2)

    def characteristic_derivative(f):
        b_a = calc_sound_speed()
        b = rel_vel_add(b_a,c*sympy.tanh(psi))
        return f.diff(t)+b*f.diff(r)

    def massage(expr, tder):
        _ = expr
        _ = characteristic_derivative(_)
        _ = _.subs(tder)
        _ = _.doit()
        _ = [_.diff(var.diff(r)) for var in [psi,p]]
        _ = sympy.Matrix(_)
        _ = _.subs(sympy.tanh(psi), beta)
        _ = _.subs(sympy.sinh(psi), beta/sympy.sqrt(1-beta**2))
        _ = _.subs(sympy.sinh(2*psi), 2*beta/(1-beta**2))
        _ = _.subs(eta, xi**2+1)
        _.simplify()
        return _
    assert all([itm == 0 for itm in
                massage(j, retrieve_planar_time_derivatives())])
    _ = massage(j, retrieve_spherical_time_derivatives())
    return _
    #return massage(j, retrieve_planar_time_derivatives())
    
    """
    _ = j
    _ = characteristic_derivative(_)
    _ = _.doit()
    _ = _.subs(retrieve_planar_time_derivatives())
    _ = _.doit()
    
    _ = [_.diff(var.diff(r)) for var in [psi,p]]
    _ = sympy.Matrix(_)
    _ = _.subs(sympy.tanh(psi), beta)
    _ = _.subs(sympy.sinh(psi), beta/sympy.sqrt(1-beta**2))
    _ = _.subs(sympy.sinh(2*psi), 2*beta/(1-beta**2))
    _ = _.subs(eta, xi**2+1)
    _.simplify()
    return _
    """

if __name__ == '__main__':

    show(locals())
