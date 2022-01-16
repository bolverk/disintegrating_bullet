import sympy
import argparse
from probe import show

t = sympy.Symbol('t', real=True) # Time
r = sympy.Symbol('r', real=True) # Radius
args = (t,r)
# Energy density, pressure and rapidity
e, p, psi, gamma = [f(*args) for f
                    in sympy.symbols('e p psi gamma',
                                     real=True,
                                     cls=sympy.Function)]
def make_lorentz_boost():

    gen = psi*sympy.Matrix([[0,1],
                            [1,0]])
    _ = sympy.exp(gen)
    _ = _.applyfunc(sympy.simplify)
    return _

def make_rest_frame_stress_energy_tensor():

    return sympy.Matrix([[e,0],
                         [0,p]])

def make_stress_energy_tensor():

    T = make_rest_frame_stress_energy_tensor()
    L = make_lorentz_boost()
    _ = L*T*L
    _.simplify()
    return _

def div(f, n):

    return (f*r**n).diff(r)/r**n

def derive_hydro_eqns(n=0):

    T = make_stress_energy_tensor()
    energy_conservation = T[0].diff(t) + div(T[1],n)
    momentum_conservation = T[2].diff(t) + div(T[3],n) - n*p/r
    _ = [energy_conservation, momentum_conservation]
    _ = sympy.Matrix(_)
    return _

eta = sympy.Symbol('eta', positive=True)
eos = e - p/(eta-1)

def calc_time_derivs(n=0):

    _ = derive_hydro_eqns(n)
    _ = _.subs(sympy.solve(eos,e,dict=True)[0])
    _ = _.doit()
    _ = sympy.solve(_,[var.diff(t) for var in (psi,p)])
    _ = {k:v.simplify() for k,v in _.items()}
    return _

def acoustic_amplitude_ratio():

    epsilon = sympy.Symbol('epsilon', positive=True)
    dpsi = sympy.Symbol('delta psi', real=True)
    df = sympy.Symbol('delta f', real=True)
    p_0 = sympy.Symbol('p_0', positive=True)
    omega = sympy.Symbol('omega', real=True)
    k = sympy.Symbol('k', positive=True)
    xi = sympy.Symbol('xi', positive=True)

    mode = epsilon*sympy.exp(sympy.I*(k*r-omega*t))
    ansatz = {p:p_0*(1+df*mode),
              psi:dpsi*mode}

    _ = derive_hydro_eqns(0)
    _ = _.subs(sympy.solve(eos,e,dict=True)[0])
    _ = _.subs(ansatz)
    _ = _.doit()
    assert all([itm==0 for itm in _.subs(epsilon,0)])
    _ = _.diff(epsilon).subs(epsilon, 0)
    _ /= p_0*mode/epsilon/(eta-1)/sympy.I
    _.simplify()
    dispersion_matrix = sympy.Matrix(
        [[itm.diff(var) for itm in _]
         for var in [df,dpsi]]).T
    dispersion_equation = dispersion_matrix.det()
    frequency_solutions = sympy.solve(dispersion_equation, omega)
    _ = dispersion_matrix.subs(omega, frequency_solutions[1])
    _ = _.subs(eta, xi**2+1)
    _.simplify()
    _ = _.nullspace()[0]
    _ = _[0]/_[1]
    _ = sympy.together(_)
    _ = _.subs(xi, sympy.sqrt(eta-1))
    return _

def calc_planar_riemann_invariant():

    q = acoustic_amplitude_ratio()
    return p*sympy.exp(q*psi)

def calc_planar_riemann_invariant_ur():

    q = acoustic_amplitude_ratio()
    return p*gamma**q

if __name__ == '__main__':

    show(locals())
