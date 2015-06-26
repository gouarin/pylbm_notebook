import numpy as np
import pylab as plt
import sympy as sp
import pyLBM

u, X = sp.symbols('u, X')

def scheme_constructor(c, s):
    dico = {
        'dim':1,
        'scheme_velocity':1.,
        'schemes':[
            {
            'velocities':range(1, 3),
            'conserved_moments':u,
            'polynomials':[1, X],
            'relaxation_parameters':[0., s],
            'equilibrium':[u, c*u],
            'init':{u:0.},
            },
        ],
        'stability':{},
    }
    return pyLBM.Scheme(dico)

c, s = 0.3, 1.5
S = scheme_constructor(c, s)
print S.amplification_matrix_relaxation

M = np.array([[1.,1.],[1.,-1.]])
invM = np.array([[.5,.5],[.5,-.5]])
iD = np.eye(2)
E = np.array([[1.,0.],[c,0.]])
S = np.array([[0.,0.],[0.,s]])
R =  np.dot(np.dot(invM, iD + np.dot(S, E - iD)), M)
print R

plt.figure(1)
plt.clf()
plt.title('Stability of the linear D1Q2')
plt.xlabel('Velocity')
plt.ylabel('Realaxation parameter')
plt.axis('equal')
plt.hold(True)
Nc, Ns = 16, 16
vc = np.linspace(0., 1.2, Nc+1)
vs = np.linspace(-.5, 2.5, Ns+1)
for c in vc:
    for s in vs:
        S = scheme_constructor(c, s)
        if S.is_monotonically_stable():
            plt.scatter([c, -c], [s, s], c = 'k', marker = 'o')
        else:
            plt.scatter([c, -c], [s, s], c = 'w', marker = 'o')
        plt.pause(1.e-5)
plt.hold(False)
plt.show()
