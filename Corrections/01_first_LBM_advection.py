import numpy as np
import pylab as plt

def mesh(N):
    xmin, xmax = 0., 1.
    dx = 1./N
    x = np.linspace(xmin-.5*dx, xmax+.5*dx, N+2)
    return x

def equilibrium(m0, c):
    return c*m0

def initialize(mesh, c, la):
    xmin, xmax = .5*np.sum(mesh[:2]), .5*np.sum(mesh[-2:])
    xm = .5*(xmin+xmax)
    L = .125*(xmax-xmin)
    m0 = np.maximum(-(mesh-xm-L)*(mesh-xm+L)/L**2, 0.)
    m1 = equilibrium(m0, c)
    f0 = np.empty(m0.shape)
    f1 = np.empty(m0.shape)
    return f0, f1, m0, m1

def f2m(f0, f1, m0, m1, la):
    m0[:] = f0 + f1
    m1[:] = la * (f1-f0)

def m2f(f0, f1, m0, m1, la):
    f0[:] = .5 * (m0 - m1/la)
    f1[:] = .5 * (m0 + m1/la)

def relaxation(m0, m1, c, s):
    m1[:] *= (1-s)
    m1[:] += s*equilibrium(m0, c)

def transport(f0, f1):
    # periodic boundary conditions
    f0[-1] = f0[1]
    f1[0] = f1[-2]
    # transport
    f0[1:-1] = f0[2:]
    f1[1:-1] = f1[:-2]

# parameters
c = .5  # velocity for the transport equation
Tf = 2. # final time
N = 128 # number of points in space
la = 1. # scheme velocity
s = 1.8 # relaxation parameter
# initialization
x = mesh(N)     # mesh
dx = x[1]-x[0]  # space step
dt = dx/la      # time step
f0, f1, m0, m1 = initialize(x, c, la)
plt.figure(1)
plt.plot(x[1:-1], m0[1:-1], 'r', label='m0')
plt.plot(x[1:-1], m1[1:-1], 'b', label='m1')
plt.title('Initial moments')
plt.legend(loc='best')
# time loops
t = 0.
while (t<Tf):
    t += dt
    relaxation(m0, m1, c, s)
    m2f(f0, f1, m0, m1, la)
    transport(f0, f1)
    f2m(f0, f1, m0, m1, la)
plt.figure(2)
plt.plot(x[1:-1], m0[1:-1], 'r', label='m0')
plt.plot(x[1:-1], m1[1:-1], 'b', label='m1')
plt.title('Final moments')
plt.legend(loc='best')
plt.show()
