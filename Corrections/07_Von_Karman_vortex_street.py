import numpy as np
import sympy as sp
import pyLBM

X, Y, LA = sp.symbols('X, Y, LA')
rho, qx, qy = sp.symbols('rho, qx, qy')

def bc_in(f, m, x, y, scheme):
    m[:, 1] = m[:, 0] * vin
    m[:, 2] = 0.
    scheme.equilibrium(m)
    scheme.m2f(m, f)

def vorticity(sol):
    ux = sol.m[0][1] / sol.m[0][0]
    uy = sol.m[0][2] / sol.m[0][0]
    V = np.abs(uy[2:,1:-1] - uy[0:-2,1:-1] - ux[1:-1,2:] + ux[1:-1,0:-2])/(2*sol.domain.dx)
    return -V

# parameters
rayon = 0.05
Re = 500
dx = 1./128  # spatial step
la = 1.      # velocity of the scheme
Tf = 500     # final time of the simulation
vin = la/20  # maximal velocity obtained in the middle of the channel
rhoo = 1.    # mean value of the density
mu = 1.e-3   # bulk viscosity
eta = rhoo*vin*2*rayon/Re  # shear viscosity
# initialization
xmin, xmax, ymin, ymax = 0., 3., 0., 1.
dummy = 3.0/(la*rhoo*dx)
s_mu = 1.0/(0.5+mu*dummy)
s_eta = 1.0/(0.5+eta*dummy)
s_q = s_eta
s_es = s_mu
s  = [0.,0.,0.,s_mu,s_es,s_q,s_q,s_eta,s_eta]
dummy = 1./(LA**2*rho)
qx2 = dummy*qx**2
qy2 = dummy*qy**2
q2  = qx2+qy2
qxy = dummy*qx*qy

print "Reynolds number: {0:10.3e}".format(Re)
print "Bulk viscosity : {0:10.3e}".format(mu)
print "Shear viscosity: {0:10.3e}".format(eta)
print "relaxation parameters: {0}".format(s)

dico = {
    'box':{'x':[xmin, xmax], 'y':[ymin, ymax], 'label':[0,2,0,0]},
    'elements':[pyLBM.Circle([.3, 0.5*(ymin+ymax)+2*dx], rayon, label=1)],
    'space_step':dx,
    'scheme_velocity':la,
    'parameters':{LA:la},
    'schemes':[
        {
            'velocities':range(9),
            'conserved_moments':[rho, qx, qy],
            'polynomials':[
                1, LA*X, LA*Y,
                3*(X**2+Y**2)-4,
                (9*(X**2+Y**2)**2-21*(X**2+Y**2)+8)/2,
                3*X*(X**2+Y**2)-5*X, 3*Y*(X**2+Y**2)-5*Y,
                X**2-Y**2, X*Y
            ],
            'relaxation_parameters':s,
            'equilibrium':[
                rho, qx, qy,
                -2*rho + 3*q2,
                rho+3/2*q2,
                -qx/LA, -qy/LA,
                qx2-qy2, qxy
            ],
            'init':{rho:rhoo, qx:0., qy:0.},
        },
    ],
    'boundary_conditions':{
        0:{'method':{0: pyLBM.bc.bouzidi_bounce_back}, 'value':bc_in},
        1:{'method':{0: pyLBM.bc.bouzidi_bounce_back}, 'value':None},
        2:{'method':{0: pyLBM.bc.neumann_vertical}, 'value':None},
    },
    'generator': pyLBM.generator.CythonGenerator,
}

sol = pyLBM.Simulation(dico)
viewer = pyLBM.viewer.matplotlibViewer
fig = viewer.Fig()
ax = fig[0]
im = ax.image(vorticity(sol).transpose(), clim = [-3., 0])
ax.ellipse([.3/dx, (0.5*(ymin+ymax)+2*dx)/dx], [rayon/dx,rayon/dx], 'r')
ax.title = 'Von Karman vortex street at t = {0:f}'.format(sol.t)

def update(iframe):
    for k in range(32):
        sol.one_time_step()
    sol.f2m()
    im.set_data(vorticity(sol).transpose())
    ax.title = 'Von Karman vortex street at t = {0:f}'.format(sol.t)
    print "image {0:d}".format(iframe)

fig.animate(update, interval=1)
fig.show()
