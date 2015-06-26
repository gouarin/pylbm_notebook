import numpy as np
import pylab as plt

def mesh(N, h):
    xmin, xmax = 0., 2.*np.pi
    dx = (xmax-xmin)/N
    x = np.linspace(xmin-.5*dx+h*dx, xmax+.5*dx+h*dx, N+2)
    return x

def equilibrium(m0, c):
    return .5*c**2*m0

def initialize(mesh, c, la):
    m0 = np.sin(mesh)
    m1 = np.zeros(mesh.shape)
    m2 = equilibrium(m0, c)
    f0 = np.empty(m0.shape)
    f1 = np.empty(m0.shape)
    f2 = np.empty(m0.shape)
    return f0, f1, f2, m0, m1, m2

def f2m(f0, f1, f2, m0, m1, m2, la):
    m0[:] = f0 + f1 + f2
    m1[:] = la * (f1 - f2)
    m2[:] = .5* la**2 * (f1 + f2)

def m2f(f0, f1, f2, m0, m1, m2, la):
    f0[:] = m0 - 2./la**2 * m2
    f1[:] = .5/la * m1 + 1/la**2 * m2
    f2[:] = -.5/la * m1 + 1/la**2 * m2

def relaxation(m0, m1, m2, c, s):
    m2[:] *= (1-s)
    m2[:] += s*equilibrium(m0, c)

def transport(f0, f1, f2, bc):
    if bc == 0:
        # anti bounce back boundary conditions
        f1[0] = -f2[1]
        f2[-1] = -f1[-2]
    else:
        # Bouzidi
        # left
        s = .5+h
        if s>0.5:
            coeff = 0.5/s
            f1[0] = -coeff*f2[1] + (1.-coeff)*f1[1]
        else:
            coeff = 2.*s
            f1[0] = -coeff*f2[1] - (1.-coeff)*f2[2]
        # right
        s = .5-h
        if s>0.5:
            coeff = 0.5/s
            f2[-1] = -coeff*f1[-2] + (1.-coeff)*f2[-2]
        else:
            coeff = 2.*s
            f2[-1] = -coeff*f1[-2] - (1.-coeff)*f1[-3]
    # transport
    f1[1:-1] = f1[:-2]
    f2[1:-1] = f2[2:]

def run(h, bc):
    # initialization
    x = mesh(N, h)     # mesh
    dx = x[1]-x[0]  # space step
    dt = dx/la      # time step
    f0, f1, f2, m0, m1, m2 = initialize(x, c, la)
    # time loops
    nt = int(Tf/dt)
    m2f(f0, f1, f2, m0, m1, m2, la)
    for k in range(nt):
        transport(f0, f1, f2, bc)
        f2m(f0, f1, f2, m0, m1, m2, la)
        relaxation(m0, m1, m2, c, s)
        m2f(f0, f1, f2, m0, m1, m2, la)
    plt.plot(x[1:-1], m0[1:-1], label='h={0:5.2f}'.format(h))

# parameters
c = .5       # velocity for the transport equation
Tf = 2*np.pi # final time
N = 32       # number of points in space
la = 1.      # scheme velocity
s = 1.5      # relaxation parameter

plt.figure(1)
plt.clf()
plt.subplot(2,1,1)
for h in [-0.25, 0., 0.25, 0.5]:
    run(h, 0)
plt.title('Anti bounce back')
plt.legend()
plt.subplot(2,1,2)
for h in [-0.25, 0., 0.25, 0.5]:
    run(h, 1)
plt.title('Bouzidi')
plt.legend()
plt.show()
