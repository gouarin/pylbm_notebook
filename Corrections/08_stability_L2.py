import numpy as np
import pylab as plt
import sympy as sp
import pyLBM

u, X = sp.symbols('u, X')

def scheme_constructor(c, s, test=False):
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
        'stability':{
            'test_L2_stability':test,
        },
    }
    return pyLBM.Scheme(dico)

def vp_plot(c):
    Nk = 100
    vkx = np.linspace(0., 2*np.pi, Nk+1)
    Ns = 50
    plt.figure(1)
    for s in np.linspace(0., 2., Ns+1):
        S = scheme_constructor(c, s)
        R = 1.
        plt.clf()
        plt.axis([-1.1, 1.1, -1.1, 1.1])
        plt.hold(True)
        plt.plot(np.cos(np.linspace(0., 2.*np.pi, 200)),
            np.sin(np.linspace(0., 2.*np.pi, 200)), 'r')
        for k in range(Nk):
            vp = S.vp_amplification_matrix((vkx[k],))
            rloc = max(np.abs(vp))
            plt.plot(vp.real, vp.imag, 'ko')
            if rloc>R+1.e-14:
                R = rloc
        if R>1+1.e-14:
            print "instable scheme for s={0:5.3f}".format(s)
        plt.hold(False)
        plt.title('eigenvalues for $s = {0:5.3f}$'.format(s))
        plt.pause(1.e-1)

if __name__ == "__main__":
    S = scheme_constructor(0.75, 1.9, test=True)
    vp_plot(.75)
