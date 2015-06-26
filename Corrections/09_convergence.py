import numpy as np
import sympy as sp
import pylab as plt
import pyLBM

u, X = sp.symbols('u, X')

def solution(x, t):
    a, b = 0., 1.
    milieu = 0.5*(a+b)
    largeur = 0.1*(b-a)
    milieu -= 0.5*c*Tf
    y = x - c*t
    return 1.0/largeur**8 * (y-milieu-largeur)**4 * (milieu-y-largeur)**4 * (abs(y-milieu)<largeur)
    #return abs(y-milieu)<largeur

def test(c, s, N):
    dico = {
        'box':{'x':[0., 1.], 'label':-1},
        'scheme_velocity':1.,
        'space_step':1./N,
        'schemes':[
            {
                'velocities':[1,2],
                'conserved_moments':u,
                'polynomials':[1, X],
                'equilibrium':[u, c*u],
                'relaxation_parameters':[0., s],
                'init':{u:(solution, (0.,))},
            },
        ],
        'generator':pyLBM.CythonGenerator,
    }
    sol = pyLBM.Simulation(dico)
    while sol.t < Tf:
        sol.one_time_step()
    sol.f2m()
    x = sol.domain.x[0][1:-1]
    y = sol.m[0][0][1:-1]
    plt.clf()
    plt.plot(x, y, 'k', x, solution(x, sol.t), 'r')
    plt.pause(1.e-5)
    return sol.domain.dx * np.linalg.norm(y - solution(x, sol.t), 1)

Tf = .4
c = 0.75
s = 1.

vdx = []
vE = []
for k in xrange(6,17):
    N = 2**k
    erreur = test(c, s, N)
    print "k={0:2d}, N={1:5d}, E={2:8.2e}".format(k, N, erreur)
    vdx.append(1./N)
    vE.append(erreur)
slope = (np.log2(vE[-1]) - np.log2(vE[-2])) / (np.log2(vdx[-1]) - np.log2(vdx[-2]))
print "order of convergence: {0:8.6f}".format(slope)
