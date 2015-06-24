import numpy as np
import sympy as sp
import pyLBM
u, X, LA = sp.symbols('u,X,LA')

# parameters
c = .5      # velocity for the transport equation
N = 128  # number of points in space
la = 1.     # scheme velocity
s = 1.8  # relaxation parameter

def solution(x, t):
    xmin, xmax = .5*np.sum(x[:2]), .5*np.sum(x[-2:])
    xm = .5*(xmin+xmax)
    L = .125*(xmax-xmin)
    y = (x-c*t)%1
    return np.maximum(-(y-xm-L)*(y-xm+L)/L**2, 0.)

dico = {
    'box': {'x': [0., 1.], 'label':-1},
    'space_step': 1./N,
    'scheme_velocity': la,
    'schemes':[{
        'velocities':[2,1],
        'conserved_moments':u,
        'polynomials':[1,LA*X],
        'relaxation_parameters':[0.,s],
        'equilibrium':[u, c*u],
        'init':{u:(solution,(0,))},
    },],
    'parameters':{LA:la},
}

# create the solution
sol = pyLBM.Simulation(dico)

# create the viewer to plot the solution
viewer = pyLBM.viewer.matplotlibViewer
fig = viewer.Fig()
ax = fig[0]
ymin, ymax = -.2, 1.2
ax.axis(0., 1., ymin, ymax)

x = sol.domain.x[0][1:-1]
l1 = ax.plot(x, sol.m[0][0][1:-1], width=2, color='b', label='D1Q2')[0]
l2 = ax.plot(x, solution(x,sol.t), width=2, color='k', label='exact')[0]

def update(iframe):
    sol.one_time_step()      # increment the solution of one time step
    l1.set_data(x, sol.m[0][0][1:-1])
    l2.set_data(x, solution(x,sol.t))
    ax.title = 'solution at t = {0:f}'.format(sol.t)
    ax.legend()

fig.animate(update)
fig.show()
