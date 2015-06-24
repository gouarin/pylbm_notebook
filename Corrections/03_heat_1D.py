import numpy as np
import sympy as sp
import pylab as plt
import pyLBM

u, X = sp.symbols('u, X')

def solution(x, t, k):
    return np.sin(k*np.pi*x)*np.exp(-k**2*np.pi**2*mu*t)

# parameters
xmin, xmax = 0., 1.
N = 128
mu = 1.
Tf = .1
dx = (xmax-xmin)/N # spatial step
la = 1./dx
s1 = 2./(1+2*mu)
s2 = 1.
k = 1 # number of the wave

dico = {
    'box':{'x':[xmin, xmax], 'label':0},
    'space_step':dx,
    'scheme_velocity':la,
    'schemes':[
        {
            'velocities':range(3),
            'conserved_moments':u,
            'polynomials':[1, X, X**2/2],
            'equilibrium':[u, 0., .5*u],
            'relaxation_parameters':[0., s1, s2],
            'init':{u:(solution,(0.,k))},
        }
    ],
    'boundary_conditions':{
        0:{'method':{0:pyLBM.bc.anti_bounce_back,}, 'value':None},
    },
}

sol = pyLBM.Simulation(dico)
x = sol.domain.x[0][1:-1]
y = sol.m[0][0][1:-1]

plt.figure(1)
plt.clf()
plt.plot(x, y,'k', label='initial')

while sol.t < Tf:
    sol.one_time_step()

plt.plot(x, y,'b', label=r'$D_1Q_3$')
plt.plot(x, solution(x, sol.t, k),'r', label='exacte')
plt.title('Heat equation t={0:5.3f}'.format(sol.t))
plt.legend()
plt.show()
