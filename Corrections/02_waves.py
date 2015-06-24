import numpy as np
import pylab as plt
import sympy as sp
import pyLBM
rho, q, X, LA = sp.symbols('rho, q, X, LA')

# parameters
c = 1.        # velocity for the transport equation
N = 128       # number of points in space
la = 1.       # scheme velocity
s = 1.9       # relaxation parameter

def rho0(x):
    return np.sin(x)

dico = {
    'box': {'x': [0., 2.*np.pi], 'label':0},
    'space_step': 2.*np.pi/N,
    'scheme_velocity': la,
    'schemes':[{
        'velocities':range(3),
        'conserved_moments':[rho,q],
        'polynomials':[1, LA*X, LA**2/2*X**2],
        'relaxation_parameters':[0., 0., s],
        'equilibrium':[rho, q, c**2/2*rho],
        'init':{rho:(rho0,), q:0.},
    },],
    'boundary_conditions':{
        0:{'method':{0:pyLBM.bc.bouzidi_anti_bounce_back,}, 'value':None},
    },
    'parameters':{LA:la},
}

sol = pyLBM.Simulation(dico)

# create the viewer to plot the solution
viewer = pyLBM.viewer.matplotlibViewer
fig = viewer.Fig()
ax = fig[0]
ymin, ymax = -1.2, 1.2
ax.axis(0., 2.*np.pi, ymin, ymax)

x = sol.domain.x[0][1:-1]
l1 = ax.plot(x, sol.m[0][0][1:-1], width=2, color='r', label=r'$\rho$')[0]
l2 = ax.plot(x, sol.m[0][1][1:-1], width=2, color='k', label=r'$q$')[0]

def update(iframe):
    sol.one_time_step()      # increment the solution of one time step
    l1.set_data(x, sol.m[0][0][1:-1])
    l2.set_data(x, sol.m[0][1][1:-1])
    ax.title = 'solution at t = {0:f}'.format(sol.t)
    ax.legend()

fig.animate(update)
fig.show()
