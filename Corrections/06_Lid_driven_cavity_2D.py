import numpy as np
import sympy as sp
import pylab as plt
import pyLBM

X, Y, LA = sp.symbols('X, Y, LA')
rho, qx, qy = sp.symbols('rho, qx, qy')

def bc(f, m, x, y, scheme):
    m[:, 1] = m[:, 0] * vup
    m[:, 2] = 0.
    scheme.equilibrium(m)
    scheme.m2f(m, f)

def plot(sol):
    plt.clf()
    pas = 2
    y, x = np.meshgrid(sol.domain.x[1][1:-1:pas], sol.domain.x[0][1:-1:pas])
    u = sol.m[0][1][1:-1:pas,1:-1:pas] / sol.m[0][0][1:-1:pas,1:-1:pas]
    v = sol.m[0][2][1:-1:pas,1:-1:pas] / sol.m[0][0][1:-1:pas,1:-1:pas]
    nv = np.sqrt(u**2+v**2)
    normu = nv.max()
    u = u / (nv+1e-5)
    v = v / (nv+1e-5)
    plt.quiver(x, y, u, v, nv, pivot='mid')
    plt.title('Solution at t={0:8.2f}'.format(sol.t))
    plt.pause(1.e-3)

# parameters
Re = 5000
dx = 1./128  # spatial step
la = 1.      # velocity of the scheme
Tf = 500     # final time of the simulation
vup = la/10    # maximal velocity obtained in the middle of the channel
rhoo = 1.    # mean value of the density
mu = 1.e-3   # bulk viscosity
eta = rhoo*vup/Re  # shear viscosity
# initialization
xmin, xmax, ymin, ymax = 0., 1., 0., 1.
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
    'box':{'x':[xmin, xmax], 'y':[ymin, ymax], 'label':[0,0,0,1]},
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
        0:{'method':{0: pyLBM.bc.bouzidi_bounce_back}, 'value':None},
        1:{'method':{0: pyLBM.bc.bouzidi_bounce_back}, 'value':bc}
    },
    'generator': pyLBM.generator.CythonGenerator,
}

sol = pyLBM.Simulation(dico)

plot(sol)
compt = 0
while (sol.t<Tf):
    sol.one_time_step()
    compt += 1
    if (compt%64==0):
        compt = 0
        sol.f2m()
        plot(sol)

plt.show()
