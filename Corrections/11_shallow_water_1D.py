import sympy as sp
import pyLBM

def Riemann_pb(x, u_L, u_R):
    xm = 0.5*(xmin+xmax)
    return u_L*(x<xm) + u_R*(x>xm) + 0.5*(u_L+u_R)*(x==xm)

# parameters
h, q, X, LA, g = sp.symbols('h, q, X, LA, g')
xmin, xmax = 0., 1.  # bounds of the domain
dx = 1./512          # spatial step
la = 2.              # velocity of the scheme
s_h, s_q = 1.7, 1.5  # relaxation parameter
Tf = 0.25            # final time

h_L, h_R, q_L, q_R = 1., .25, 0.10, 0.10
yminh, ymaxh, yminq, ymaxq = 0.2, 1.2, 0., .5

dico = {
    'box':{'x':[xmin, xmax], 'label':0},
    'space_step':dx,
    'scheme_velocity':la,
    'schemes':[
        {
            'velocities':[1,2],
            'conserved_moments':h,
            'polynomials':[1, LA*X],
            'relaxation_parameters':[0, s_h],
            'equilibrium':[h, q],
            'init':{h:(Riemann_pb, (h_L, h_R))},
        },
        {
            'velocities':[1,2],
            'conserved_moments':q,
            'polynomials':[1, LA*X],
            'relaxation_parameters':[0, s_q],
            'equilibrium':[q, q**2/h+.5*g*h**2],
            'init':{q:(Riemann_pb, (q_L, q_R))},
        },
    ],
    'boundary_conditions':{
        0:{'method':{0: pyLBM.bc.neumann, 1: pyLBM.bc.neumann}, 'value':None},
    },
    'parameters':{LA:la, g:1.},
}

sol = pyLBM.Simulation(dico)
viewer = pyLBM.viewer.matplotlibViewer
fig = viewer.Fig(1,2)
ax_h = fig[0,0]
ax_q = fig[0,1]

x = sol.domain.x[0][1:-1]
lh = ax_h.plot(x, sol.m[0][0][1:-1], width=2, color='r', label=r'$h$')[0]
lq = ax_q.plot(x, sol.m[1][0][1:-1], width=2, color='b', label=r'$q$')[0]
ax_h.axis(xmin, xmax, yminh, ymaxh)
ax_q.axis(xmin, xmax, yminq, ymaxq)
ax_h.title = r'$h$ at t = {0:f}'.format(sol.t)
ax_q.title = r'$q$ at t = {0:f}'.format(sol.t)

def update(iframe):
    if (sol.t<Tf):
        sol.one_time_step()
        lh.set_data(x, sol.m[0][0][1:-1])
        lq.set_data(x, sol.m[1][0][1:-1])
        ax_h.title = r'$h$ at t = {0:f}'.format(sol.t)
        ax_q.title = r'$q$ at t = {0:f}'.format(sol.t)

fig.animate(update, interval = 1)
fig.show()
