import numpy as np
import sympy as sp
import pyLBM
u, X, LA = sp.symbols('u,X,LA')

# parameters
N = 128 # number of points in space
la = 1. # scheme velocity
s = 1.8 # relaxation parameter

def u0(x):
    ug, ud = 0.25, -0.15
    xmin, xmax = .5*np.sum(x[:2]), .5*np.sum(x[-2:])
    xc = xmin + .5*(xmax-xmin)
    return ug*(x<xc) + ud*(x>xc) + .5*(ug+ud)*(x==xc)

dico = {
    'box': {'x': [0., 1.], 'label':0},
    'space_step': 1./N,
    'scheme_velocity': la,
    'schemes':[{
        'velocities':[2,1],
        'conserved_moments':u,
        'polynomials':[1,LA*X],
        'relaxation_parameters':[0.,s],
        'equilibrium':[u, u**2/2],
        'init':{u:(u0,)},
    },],
    'boundary_conditions':{
        0:{'method':{0:pyLBM.bc.neumann,}, 'value':None},
    },
    'parameters':{LA:la},
}

sol = pyLBM.Simulation(dico)

# create the viewer to plot the solution
viewer = pyLBM.viewer.matplotlibViewer
fig = viewer.Fig()
ax = fig[0]
ymin, ymax = -.2, .3
ax.axis(0., 1., ymin, ymax)

x = sol.domain.x[0][1:-1]
l1 = ax.plot(x, sol.m[0][0][1:-1], width=2, color='b', label='D1Q2')[0]

def update(iframe):
    for k in range(16):
        sol.one_time_step()      # increment the solution of one time step
    sol.f2m()
    l1.set_data(x, sol.m[0][0][1:-1])
    ax.title = 'solution at t = {0:f}'.format(sol.t)
    ax.legend()

fig.animate(update)
fig.show()
