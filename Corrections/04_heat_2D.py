import numpy as np
import sympy as sp
import pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pyLBM

u, X, Y = sp.symbols('u, X, Y')

def solution(x, y, t, k, l):
    return np.sin(k*np.pi*x)*np.sin(l*np.pi*y)*np.exp(-(k**2+l**2)*np.pi**2*mu*t)

def plot(i, j, z, title):
    im = axarr[i,j].imshow(z)
    divider = make_axes_locatable(axarr[i, j])
    cax = divider.append_axes("right", size="20%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax, format='%6.0e')
    axarr[i, j].xaxis.set_visible(False)
    axarr[i, j].yaxis.set_visible(False)
    axarr[i, j].set_title(title)


# parameters
xmin, xmax, ymin, ymax = 0., 1., 0., 1.
N = 128
mu = 1.
Tf = .1
dx = (xmax-xmin)/N # spatial step
la = 1./dx
s1 = 2./(1+4*mu)
s2 = 1.
k, l = 2, 1 # number of the wave

dico = {
    'box':{'x':[xmin, xmax], 'y':[ymin, ymax], 'label':0},
    'space_step':dx,
    'scheme_velocity':la,
    'schemes':[
        {
            'velocities':range(5),
            'conserved_moments':u,
            'polynomials':[1, X, Y, (X**2+Y**2)/2, (X**2-Y**2)/2],
            'equilibrium':[u, 0., 0., .5*u, 0.],
            'relaxation_parameters':[0., s1, s1, s2, s2],
            'init':{u:(solution,(0.,k,l))},
        }
    ],
    'boundary_conditions':{
        0:{'method':{0:pyLBM.bc.anti_bounce_back,}, 'value':None},
    },
    'generator':pyLBM.generator.CythonGenerator,
}

sol = pyLBM.Simulation(dico)
x = sol.domain.x[0][1:-1]
y = sol.domain.x[1][1:-1]

f, axarr = plt.subplots(2, 2)
f.suptitle('Heat equation', fontsize=20)

plot(0, 0, sol.m[0][0][1:-1,1:-1].copy(), 'initial')

while sol.t < Tf:
    sol.one_time_step()

sol.f2m()
z = sol.m[0][0][1:-1,1:-1]
ze = solution(x[:,np.newaxis], y[np.newaxis,:], sol.t, k, l)
plot(1, 0, z, 'final')
plot(0, 1, ze, 'exact')
plot(1, 1, z-ze, 'error')

plt.show()
