import sympy
from probe import show

from rhd import (
    derive_hydro_eqns,
    r,
    t,
    e,
    p,
    psi,
    eta,
    eos)

df = sympy.Symbol('delta f', real=True)
dpsi = sympy.Symbol('delta psi', real=True)
omega = sympy.Symbol('omega', real=True)
k = sympy.Symbol('k', positive=True)

def calc_acoustic_equations():

    epsilon = sympy.Symbol('epsilon', positive=True)
    p_0 = sympy.Symbol('p_0', positive=True)
    xi = sympy.Symbol('xi', positive=True)

    mode = epsilon*sympy.exp(sympy.I*(k*r-omega*t))
    ansatz = {p:p_0*(1+df*mode),
              psi:dpsi*mode}

    _ = derive_hydro_eqns()
    _ = _.subs(sympy.solve(eos,e,dict=True)[0])
    _ = _.subs(ansatz)
    _ = _.doit()
    assert all([itm==0 for itm in _.subs(epsilon,0)])
    _ = _.diff(epsilon).subs(epsilon, 0)
    _ /= p_0*mode/epsilon/(eta-1)/sympy.I
    _.simplify()
    return _

def calc_dispersion_matrix():

    _ = calc_acoustic_equations()
    return sympy.Matrix(
        [[itm.diff(var)
          for var in [df, dpsi]]
         for itm in _])

def calc_dispersion_equation():

    _ = calc_dispersion_matrix()
    return _.det()

def calc_sound_speed():

    _ = calc_dispersion_equation()
    return sympy.solve(_,omega)[1]/k

def calc_amplitude_ratio():

    xi = sympy.Symbol('xi', positive=True)

    M = calc_dispersion_matrix()
    eqn = M.det()
    omega_sol = sympy.solve(eqn, omega)[1]
    _ = M
    _ = _.subs(omega, omega_sol)
    _ = _.subs(eta, xi**2+1)
    _.simplify()
    _ = _.nullspace()
    _ = _[0]
    _ = _[0]/_[1]
    _ = sympy.together(_)
    _ = _.subs(xi, sympy.sqrt(eta-1))
    return _

if __name__ == '__main__':

    show(locals())

