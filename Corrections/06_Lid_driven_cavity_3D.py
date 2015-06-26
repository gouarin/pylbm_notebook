import numpy as np
import sympy as sp
import pylab as plt
import pyLBM

X, Y, Z, LA = sp.symbols('X, Y, Z, LA')
rho, qx, qy, qz = sp.symbols('rho, qx, qy, qz')

def bc(f, m, x, y, z, scheme):
    m[:, 3] = m[:, 0] * vup
    m[:, 5] = 0.
    m[:, 7] = 0.
    scheme.equilibrium(m)
    scheme.m2f(m, f)

def plot(sol):
    plt.clf()
    pas = 2
    nz = int(sol.domain.N[1] / 2) + 1
    y, x = np.meshgrid(sol.domain.x[2][1:-1:pas], sol.domain.x[0][1:-1:pas])
    u = sol.m[0][3][1:-1:pas,nz,1:-1:pas] / sol.m[0][0][1:-1:pas,nz,1:-1:pas]
    v = sol.m[0][7][1:-1:pas,nz,1:-1:pas] / sol.m[0][0][1:-1:pas,nz,1:-1:pas]
    nv = np.sqrt(u**2+v**2)
    normu = nv.max()
    u = u / (nv+1e-5)
    v = v / (nv+1e-5)
    plt.quiver(x, y, u, v, nv, pivot='mid')
    plt.title('Solution at t={0:9.3f}'.format(sol.t))
    plt.pause(1.e-3)

# parameters
Re = 2000
dx = 1./128  # spatial step
la = 1.      # velocity of the scheme
Tf = 10      # final time of the simulation
vup = la/10  # maximal velocity obtained in the middle of the channel
rhoo = 1.    # mean value of the density
eta = rhoo*vup/Re  # shear viscosity
# initialization
xmin, xmax, ymin, ymax, zmin, zmax = 0., 1., 0., 1., 0., 1.
dummy = 3.0/(la*rhoo*dx)

s1 = 1.6
s2 = 1.2
s4 = 1.6
s9 = 1./(.5+dummy*eta)
s11 = s9
s14 = 1.2
s  = [0, s1, s2, 0, s4, 0, s4, 0, s4, s9, s9, s11, s11, s11, s14]

r = X**2+Y**2+Z**2

print "Reynolds number: {0:10.3e}".format(Re)
print "Shear viscosity: {0:10.3e}".format(eta)

dico = {
    'box':{
        'x':[xmin, xmax],
        'y':[ymin, ymax],
        'z':[zmin, zmax],
        'label':[0,0,0,0,0,1]
    },
    'space_step':dx,
    'scheme_velocity':la,
    'parameters':{LA:la},
    'schemes':[
        {
            'velocities':range(7) + range(19,27),
            'conserved_moments':[rho, qx, qy, qz],
            'polynomials':[
                1,
                r - 2, .5*(15*r**2-55*r+32),
                X, .5*(5*r-13)*X,
                Y, .5*(5*r-13)*Y,
                Z, .5*(5*r-13)*Z,
                3*X**2-r, Y**2-Z**2,
                X*Y, Y*Z, Z*X,
                X*Y*Z
            ],
            'relaxation_parameters':s,
            'equilibrium':[
                rho,
                -rho + qx**2 + qy**2 + qz**2,
                -rho,
                qx,
                -7./3*qx,
                qy,
                -7./3*qy,
                qz,
                -7./3*qz,
                1./3*(2*qx**2-(qy**2+qz**2)),
                qy**2-qz**2,
                qx*qy,
                qy*qz,
                qz*qx,
                0
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
    if (compt%16==0):
        compt = 0
        sol.f2m()
        plot(sol)

plt.show()
