import sympy
import argparse
from probe import show
from caching import memory, luggage

nu = sympy.Symbol('nu', negative=False)
t = sympy.Symbol('t', real=True) # Time
r = sympy.Symbol('r', real=True) # Radius
c = sympy.Symbol('c', positive=True) # Speed of light
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

def div(f):

    return (f*r**nu).diff(r)/r**nu

def derive_hydro_eqns(n=0):

    T = make_stress_energy_tensor()
    energy_conservation = T[0].diff(t)/c + div(T[1])
    momentum_conservation = T[2].diff(t)/c + div(T[3]) - nu*p/r
    _ = [energy_conservation, momentum_conservation]
    _ = sympy.Matrix(_)
    return _

eta = sympy.Symbol('eta', positive=True)
eos = e - p/(eta-1)

@luggage.memory
def calc_time_derivs():

    _ = derive_hydro_eqns()
    _ = _.subs(sympy.solve(eos,e,dict=True)[0])
    _ = _.doit()
    _ = sympy.solve(_,[var.diff(t) for var in (psi,p)])
    _ = {k:v.simplify() for k,v in _.items()}
    return _

if __name__ == '__main__':

    show(locals())
