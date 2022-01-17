import sympy
from probe import show
from my_logging import logger
from caching import luggage

from acoustic import (
    calc_amplitude_ratio,
    calc_sound_speed)
from rhd import (
    nu,
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
def calc_spherical_riemann_invariant():

    beta = sympy.Symbol('beta', positive=True)
    gamma = sympy.Symbol('gamma', positive=True)
    xi = sympy.Symbol('xi', positive=True)
    zeta = sympy.Symbol('zeta', real=True)

    logger.debug('begin planar riemann invariant calculation')
    j_p = calc_planar_riemann_invariant()
    logger.debug('end planar riemann invariant calculation')
    def rel_vel_add(v1, v2):

        return (v1+v2)/(1+v1*v2/c**2)

    def characteristic_derivative(f):
        b_a = calc_sound_speed()
        b = rel_vel_add(b_a,c*sympy.tanh(psi))
        return f.diff(t)+b*f.diff(r)

    def massage(expr):
        logger.debug('inside massage')
        _ = expr
        logger.debug('load expression')
        _ = characteristic_derivative(_)
        logger.debug('char. derivative')
        _ = _.subs(calc_time_derivs())
        logger.debug('time derivative')
        _ = _.doit()
        logger.debug('doit')
        _ = [_.diff(var.diff(r)) for var in [psi,p]]+\
            [_.subs({var.diff(r):0 for var in [psi,p]})]
        _ = sympy.Matrix(_)
        logger.debug('make matrix')
        _ = _.subs(sympy.tanh(psi), beta)
        _ = _.subs(sympy.sinh(psi), beta/sympy.sqrt(1-beta**2))
        _ = _.subs(sympy.sinh(2*psi), 2*beta/(1-beta**2))
        _ = _.subs(sympy.cosh(2*psi),2*sympy.cosh(psi)**2-1)
        _ = _.subs(sympy.cosh(psi),1/sympy.sqrt(1-beta**2))
        #_ = _.subs(beta,sympy.sqrt(1-1/gamma**2))
        logger.debug('get rid of hyp. trig.')
        _ = _.subs(eta, xi**2+1)
        logger.debug('before simplify')
        _.simplify()
        _ = _.subs(xi, sympy.sqrt(eta-1))
        logger.debug('after simplify')
        return _
    j_s = j_p + zeta*sympy.log(r)
    _ = massage(j_s)
    assert _[0] == 0 and _[1]==0
    _ = _[2]
    _ = sympy.solve(_, zeta)[0].subs(beta, 1)
    return j_s.subs(zeta, _)

if __name__ == '__main__':

    show(locals())
