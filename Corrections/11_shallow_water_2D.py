import numpy as np
import sympy as sp
import pyLBM

X, Y, LA = sp.symbols('X, Y, LA')
h, qx, qy = sp.symbols('h, qx, qy')

def h0(x, y):
    return h_l * np.ones((x.size, y.size), dtype='float64') \
        + (h_h-h_l) * ((x-0.5*(xmin+xmax))**2+(y-0.5*(ymin+ymax))**2 < 0.25**2)

# parameters
dx = 1./128  # spatial step
la = 4.      # velocity of the scheme
h_l = 1.     # low value of the water height
h_h = 2.     # high value of the water height
L = 2        # size of the domain
g = 1.       # gravity
s_h1 = 2.
s_h2 = 1.5
s_q1 = 1.5
s_q2 = 1.2
# initialization
xmin, xmax, ymin, ymax = -.5*L, .5*L, -.5*L, .5*L
s_h = [0., s_h1, s_h1, s_h2]
s_q = [0., s_q1, s_q1, s_q2]

vitesse = range(1,5)
polynomes = [1, LA*X, LA*Y, X**2-Y**2]

dico = {
    'box':{'x':[xmin, xmax], 'y':[ymin, ymax], 'label':-1},
    'space_step':dx,
    'scheme_velocity':la,
    'parameters':{LA:la},
    'schemes':[
        {
            'velocities':vitesse,
            'conserved_moments':h,
            'polynomials':polynomes,
            'relaxation_parameters':s_h,
            'equilibrium':[h, qx, qy, 0.],
            'init':{h:(h0,)},
        },
        {
            'velocities':vitesse,
            'conserved_moments':qx,
            'polynomials':polynomes,
            'relaxation_parameters':s_q,
            'equilibrium':[qx, qx**2/h + 0.5*g*h**2, qx*qy/h, 0.],
            'init':{qx:0.},
        },
        {
            'velocities':vitesse,
            'conserved_moments':qy,
            'polynomials':polynomes,
            'relaxation_parameters':s_q,
            'equilibrium':[qy, qx*qy/h, qy**2/h + 0.5*g*h**2, 0.],
            'init':{qy:0.},
        },
    ],
    'generator': pyLBM.generator.CythonGenerator,
}

sol = pyLBM.Simulation(dico)
viewer = pyLBM.viewer.matplotlibViewer
fig = viewer.Fig()
ax = fig[0]
im = ax.image(sol.m[0][0][1:-1,1:-1].transpose(), clim = [.5*h_l, 1.1*h_h])
ax.title = 'water height at t = {0:f}'.format(sol.t)

def update(iframe):
    for k in range(16):
        sol.one_time_step()
    sol.f2m()
    im.set_data(sol.m[0][0][1:-1,1:-1].transpose())
    ax.title = 'water height at t = {0:f}'.format(sol.t)

fig.animate(update, interval=1)
fig.show()
