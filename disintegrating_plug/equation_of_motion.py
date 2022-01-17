import sympy
from probe import show
from caching import luggage

from rhd import (
    nu, c, t, r, p, psi, eta
    )

from riemann_invariant import calc_spherical_riemann_invariant

M = sympy.Symbol('M', positive=True) # Bullet inertia
M_1 = sympy.Symbol('M_1', positive=True)
gamma = sympy.Function('gamma')
A = sympy.Symbol('A')
p_i, gamma_i, t_i = sympy.symbols('p_i gamma_i t_i', positive=True)
alpha = sympy.Symbol('alpha', positive=True)

def prepare_ri():

    y1, y2, y3 = sympy.symbols('y{1:4}')
    yvec = (y1,y2,y3)

    aliases = {t:y1,
               p:y2,
               gamma(t):y3}

    #factors = [t,p,gamma(t)]

    _ = calc_spherical_riemann_invariant()
    _ = _.subs(psi, sympy.log(gamma(t)))
    _ = _.subs(aliases).subs(r, t*c).subs(aliases)
    _ = [(_.diff(var)*var).simplify() for var in yvec]
    res = 1
    for pli, yn in zip(_, yvec):
        res *= yn**pli
    res = res.subs(y1,t/t_i)
    res = res.subs(y3,gamma(t)/gamma_i)
    res = res.subs(y2,p/p_i)
    return res

def make_equation_of_motion():

    return c*M*gamma(t).diff() - A*p

@luggage.memory
def calc_asymptotic_gamma():

    xi = sympy.Symbol('xi', positive=True)
    Gamma = sympy.Symbol('Gamma', positive=True)

    ri = prepare_ri()
    _ = make_equation_of_motion()
    _ = _.subs(A, (alpha*t*c)**2)
    _ = _.subs(M, M_1*(p/p_i)**(1-1/eta))
    _ = _.subs(sympy.solve(ri-1, p,dict=True)[0])
    _ = _.subs(gamma(t).diff(t), gamma(t)/t)
    _ = _.subs(gamma(t), Gamma)
    _ = sympy.expand_power_base(_,force=True)
    _ = _.subs(eta, xi**2+1)
    _ = sympy.solve(_, Gamma)[0]
    _ = sympy.expand_power_base(_, force=True).simplify()
    _ = _.subs(xi,sympy.sqrt(eta-1))
    return _

def eval_lf_pli():

    _ = calc_asymptotic_gamma()
    _ = sympy.log(_)
    _ = _.diff(t)
    _ *= t
    _ = _.subs(eta, sympy.Rational(4,3))
    _ = _.subs(nu, 2)
    _ = _.n()
    return _

if __name__ == '__main__':

    show(locals())
