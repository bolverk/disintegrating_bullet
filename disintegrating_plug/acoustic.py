import sympy
import argparse

from rhd import (
    derive_hydro_eqns,
    r,
    t,
    e,
    p,
    psi,
    eta,
    eos)

def calc_acoustic_equations():

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

    _ = derive_hydro_eqns()
    _ = _.subs(sympy.solve(eos,e,dict=True)[0])
    _ = _.subs(ansatz)
    _ = _.doit()
    assert all([itm==0 for itm in _.subs(epsilon,0)])
    _ = _.diff(epsilon).subs(epsilon, 0)
    _ /= p_0*mode/epsilon/(eta-1)/sympy.I
    _.simplify()
    return _

def parse_input():

    parser = argparse.ArgumentParser(
        description='try internal functions')
    parser.add_argument('func_name',
                        type=str,
                        help='name of internal function to be tested')
    args = parser.parse_args()
    return args.func_name

if __name__ == '__main__':

    func_name = parse_input()

    sympy.pprint(eval(func_name+'()'))

