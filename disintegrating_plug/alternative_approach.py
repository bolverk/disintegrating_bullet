from box import Box
import sympy
from probe import show
from rhd import (psi, p, eta, c, t, r, nu)
from riemann_invariant import (
    calc_planar_riemann_invariant,
    calc_spherical_riemann_invariant
    )
from my_logging import logger
from caching import memory
from equation_of_motion import alpha

def acceleration_transformation():

    sympy.var('a t c gamma beta', positive=True)
    rf = sympy.Matrix([a*t**2/2,c*t])
    boost = gamma*sympy.Matrix([[1,sympy.sqrt(1-1/gamma**2)],
                                [sympy.sqrt(1-1/gamma**2),1]])
    lf = boost*rf
    lf.simplify()
    _ = (lf[0].diff(t)/lf[1].diff(t)).diff(t)/lf[1].diff(t)
    _ = _.series(t,0,1)
    _ = _.removeO()
    _ = _.simplify()
    return _

@memory.cache
def traditional_planar_breakup():

    sympy.var('M gamma t A p_i gamma_i alpha w c n n_1 w_1 eta a mu',
              positive=True)

    def prepare_riemann_invariant():

        _ = calc_planar_riemann_invariant()
        _ = _.subs(psi, sympy.log(gamma))
        _ = _ - _.subs({gamma:gamma_i, p:p_i})
        return _

    def make_eqns():

        return Box(
            {'eom':M*gamma*c/t - p*A,
             'ri':prepare_riemann_invariant(),
             'baryon':A*(n_1*w_1-n*w),
             'adiabatic':p/p_i-(n/n_1)**eta}
            )

    def make_gamma_vs_t():

        xi = sympy.Symbol('xi', positive=True)

        _ = u.eom
        _ = _.subs(M,A*w*p/c**2)
        _ = _/p/A
        _ = _.simplify()
        _ = _.subs(sympy.solve(u.baryon, w, dict=True)[0])
        _ = _.subs(sympy.solve(u.adiabatic, n, dict=True)[0])
        _ = _.subs(sympy.solve(u.ri, p, dict=True)[0])
        _ = _.subs(eta ,xi**2+1)
        _ = sympy.expand_power_base(_, force=True)
        _ = _.simplify()
        _ = sympy.solve(_, gamma)[0]
        _ = sympy.expand_power_base(_, force=True)
        _ = _.simplify()
        _ = _.subs(xi, sympy.sqrt(eta-1))
        return _

    @memory.cache
    def calc_t_breakup():

        xi = sympy.Symbol('xi', positive=True)

        _ = (t/gamma)**2*(a/w)*(mu*n*c**2)/p
        _ = _.subs(sympy.solve(u.baryon, w, dict=True)[0])
        _ = _.subs(sympy.solve(u.adiabatic, n, dict=True)[0])
        _ = _.subs(sympy.solve(u.ri, p, dict=True)[0])
        _ = _.subs(a, c*gamma/t)
        _ = _.subs(gamma, u.gamma_vs_t)
        _ = _.subs(eta, xi**2+1)
        _ = sympy.expand_power_base(_, force=True)
        _ = _.simplify()
        _ = sympy.solve(_-1, t)[0]
        _ = sympy.expand_power_base(_, force=True)
        _ = _.simplify()
        _ = _.subs(xi, sympy.sqrt(eta-1))
        return _

    def calc_gamma_breakup():

        xi = sympy.Symbol('xi', positive=True)

        _ = u.gamma_vs_t.subs(t, u.t_breakup)
        _ = _.subs(p_i, mu*gamma_i*n_1*c**2)
        _ = _.subs(eta, xi**2+1)
        _ = sympy.expand_power_base(_, force=True)
        _ = _.simplify()
        _ = _.subs(xi, sympy.sqrt(eta-1))
        return _

    logger.debug('begin alternative_approach')
    u = make_eqns()
    logger.debug('finished make_eqns')
    u['gamma_vs_t'] = make_gamma_vs_t()
    logger.debug('finished make_gamma_vs_t')
    u['t_breakup'] = calc_t_breakup()
    logger.debug('finished calc_t_breakup')
    u['gamma_breakup'] = calc_gamma_breakup()
    logger.debug('finished calc_gamma_breakup')
    return u

def traditional_spherical_breakup():

    sympy.var('gamma w n w_1 n_1 p_i M a mu gamma_i t_i',
              positive=True)

    def make_ri():

        xi = sympy.Symbol('xi', positive=True)

        _ = calc_spherical_riemann_invariant()
        _ = _.subs(psi, sympy.log(gamma))
        _ = _.subs(p,xi)
        _ = _.subs(r, t*c)
        _ = _ - _.subs({gamma:gamma_i,
                        t:t_i,
                        xi:p_i})
        _ = _.subs(xi, p)
        return _

    def make_eom():

        rhs = M*gamma*c/t
        area = (alpha*t*c)**2
        lhs = p*area
        _ = sympy.log(rhs) - sympy.log(lhs)
        _ = sympy.expand_power_base(_, force=True)
        _ = _.simplify()
        _ = sympy.expand(_)
        return _

    def make_baryon():

        rhs = alpha**2*c**2*t**2*w*n
        lhs = rhs.subs({t:t_i,w:w_1,n:n_1})
        _ = sympy.log(lhs) - sympy.log(rhs)
        _ = sympy.expand(_)
        return _

    def make_adiabatic():

        rhs = sympy.log(p) - eta*sympy.log(n)
        lhs = rhs.subs({p:p_i,n:n_1})
        _ = sympy.log(rhs) - sympy.log(lhs)
        _ = sympy.expand(_)
        return _

    def make_eqns():

        return Box(
            {'eom':make_eom(),
             'ri':make_ri(),
             'baryon':make_baryon(),
             'adiabatic':make_adiabatic()})

    def make_gamma_vs_t():

        xi = sympy.Symbol('xi', positive=True)

        _ = u.eom
        _ = _.subs(M, alpha**2*t**2*c**2*w*p/c**2)
        _ = _.subs(sympy.solve(u.baryon, w, dict=True)[0])
        _ = sympy.expand(_)
        _ = _.subs(sympy.solve(u.adiabatic, sympy.log(n), dict=True)[0])
        _ = sympy.expand(_)
        _ = _.subs(sympy.solve(u.ri, sympy.log(p), dict=True)[0])
        _ = sympy.expand(_)
        _ = _.subs(nu, 2)
        _ = sympy.expand(_)
        _ = sympy.solve(_, sympy.log(gamma))[0]
        _ = sympy.expand(_)
        return _

    def calc_t_breakup():

        xi = sympy.Symbol('xi', positive=True)

        _ = (t/gamma)**2*(a/w)*(mu*n*c**2)/p
        _ = _.subs(a, c*gamma/t)
        _ = _.subs(p, xi)
        _ = sympy.log(_)
        _ = sympy.expand(_)
        _ = _.subs(xi, p)
        _ = _.subs(sympy.solve(u.baryon, w, dict=True)[0])
        _ = sympy.expand(_)
        _ = _.subs(sympy.solve(u.adiabatic, sympy.log(n), dict=True)[0])
        _ = sympy.expand(_)
        _ = _.subs(sympy.solve(u.ri, sympy.log(p), dict=True)[0])
        _ = _.subs(nu, 2)
        _ = sympy.expand(_)
        _ = _.subs(sympy.log(gamma), u.gamma_vs_t)
        _ = sympy.expand(_)
        _ = _.subs(eta, sympy.Rational(4,3))
        _ = _.n()
        _ = sympy.solve(_, sympy.log(t))[0]
        #_ = -_.subs(sympy.log(t),0)/_.subs(sympy.log(t),xi).diff(xi)
        _ = sympy.expand(_)
        return _

    def calc_gamma_breakup():

        _ = u.gamma_vs_t
        _ = sympy.expand(_)
        _ = _.subs(sympy.log(t), u.t_breakup)
        _ = _.subs(p_i, n_1*mu*c**2*gamma_i)
        _ = _.subs(eta, sympy.Rational(4,3))
        _ = _.subs(w_1, t_i*c/gamma_i)
        _ = _.n()
        _ = sympy.expand(_)
        _ = _.simplify()
        return _

    logger.debug('begin spherical breakup calculation')
    u = make_eqns()
    logger.debug('finished make eqns')
    u['gamma_vs_t'] = make_gamma_vs_t()
    logger.debug('finished gamma_vs_t')
    u['t_breakup'] = calc_t_breakup()
    logger.debug('finished t breakup')
    u['gamma_breakup'] = calc_gamma_breakup()
    logger.debug('finished gamma breakup')
    return [u.gamma_breakup,
            u.t_breakup.subs({p_i:gamma_i*mu*n_1*c**2,
                              w_1:t_i*c/gamma_i}).simplify()]

def calc_planar_breakup():

    sympy.var('p_i gamma gamma_i p_t gamma_t rho_t rho_i xi',
              positive=True)

    initial = {p:p_i,
               gamma:gamma_i}
    final = {p:p_t,
             gamma:gamma_t}

    ri = calc_planar_riemann_invariant().subs(psi, sympy.log(gamma))
    eqn1 = ri.subs(initial) - ri.subs(final)
    eqn2 = sympy.log(p_t) - sympy.log(p_i) - eta*(sympy.log(rho_t)-sympy.log(rho_i))
    _ = [eqn1, eqn2]
    _ = sympy.Matrix(_)
    _ = _.subs(p_i, gamma_i*rho_i)
    _ = _.subs(rho_t, p_t)
    _ = _.subs(eta, sympy.Rational(4,3))
    _ = sympy.solve(_,[gamma_t, p_t])[0]
    sol = _
    return sol

if __name__ == '__main__':

    show(locals())
