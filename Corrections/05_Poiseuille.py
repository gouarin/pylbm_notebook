import numpy as np
import sympy as sp
import pylab as plt
import pyLBM

X, Y, LA = sp.symbols('X, Y, LA')
rho, qx, qy = sp.symbols('rho, qx, qy')

def bc(f, m, x, y, scheme):
    # condition for Cython generator
    m[:, 1] = rhoo * vmax * (1.-4.*y**2/W**2)
    m[:, 2] = 0.
    # condition for Numpy generator
    # m[1] = rhoo * vmax * (1.-4.*y**2/W**2)
    # m[2] = 0.
    scheme.equilibrium(m)
    scheme.m2f(m, f)

def plot_coupe(sol):
    ax1.cla()
    ax2.cla()
    mx = int(sol.domain.N[0]/2-1)
    my = int(sol.domain.N[1]/2-1)
    x = sol.domain.x[0][1:-1]
    y = sol.domain.x[1][1:-1]
    u = sol.m[0][1][1:-1,1:-1] / rhoo
    for i in [0,mx,-1]:
        ax1.plot(y+x[i], u[i, :], 'b')
    for j in [0,my,-1]:
        ax1.plot(x+y[j], u[:,j], 'b')
    ax1.set_ylabel('velocity', color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
    ax1.set_ylim(-.5*rhoo*vmax, 1.5*rhoo*vmax)
    p = sol.m[0][0][1:-1,my] * la**2 / 3.0
    p -= np.average(p)
    ax2.plot(x, p, 'r')
    ax2.set_ylabel('pressure', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    ax2.set_ylim(pressure_gradient*L, -pressure_gradient*L)
    plt.title('Poiseuille flow at t = {0:f}'.format(sol.t))
    plt.draw()
    plt.pause(1.e-3)

if __name__ == "__main__":
    # parameters
    dx = 1./128  # spatial step
    la = 1.      # velocity of the scheme
    Tf = 50      # final time of the simulation
    L = 2        # length of the domain
    W = 1        # width of the domain
    vmax = 0.1   # maximal velocity obtained in the middle of the channel
    rhoo = 1.    # mean value of the density
    mu = 1.e-2   # bulk viscosity
    eta = 1.e-2  # shear viscosity
    pressure_gradient = - vmax * 8.0 / W**2 * eta
    # initialization
    xmin, xmax, ymin, ymax = 0.0, L, -0.5*W, 0.5*W
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

    dico = {
        'box':{'x':[xmin, xmax], 'y':[ymin, ymax], 'label':0},
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
            0:{'method':{0: pyLBM.bc.bouzidi_bounce_back}, 'value':bc}
        },
        'generator': pyLBM.generator.CythonGenerator,
    }

    sol = pyLBM.Simulation(dico)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    plot_coupe(sol)
    compt = 0
    while (sol.t<Tf):
        sol.one_time_step()
        compt += 1
        if (compt%32==0):
            compt = 0
            sol.f2m()
            plot_coupe(sol)

    ny = int(sol.domain.N[1]/2)
    num_pressure_gradient = (sol.m[0][0][-2,ny] - sol.m[0][0][1,ny]) / (xmax-xmin) * la**2/ 3.0
    print "Exact pressure gradient    : {0:10.3e}".format(pressure_gradient)
    print "Numerical pressure gradient: {0:10.3e}".format(num_pressure_gradient)
    plt.show()
