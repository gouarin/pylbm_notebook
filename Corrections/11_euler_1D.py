import sympy as sp
import pyLBM

def Riemann_pb(x, u_L, u_R):
    xm = 0.5*(xmin+xmax)
    return u_L*(x<xm) + u_R*(x>xm) + 0.5*(u_L+u_R)*(x==xm)

# parameters
rho, q, E, X, LA = sp.symbols('rho,q,E,X,LA')
# parameters
gamma = 1.4
xmin, xmax = 0., 1.
dx = 1.e-3 # spatial step
la = 3. # velocity of the scheme
rho_L, rho_R, p_L, p_R, u_L, u_R = 1., 1./8., 1., 0.1, 0., 0.
q_L = rho_L*u_L
q_R = rho_R*u_R
E_L = rho_L*u_L**2 + p_L/(gamma-1.)
E_R = rho_R*u_R**2 + p_R/(gamma-1.)
Tf = 0.14 # final time
s_rho, s_q, s_E = 1.9, 1.5, 1.4

dico = {
    'box':{'x':[xmin, xmax], 'label':0},
    'space_step':dx,
    'scheme_velocity':la,
    'schemes':[
        {
            'velocities':[1,2],
            'conserved_moments':rho,
            'polynomials':[1, LA*X],
            'relaxation_parameters':[0, s_rho],
            'equilibrium':[rho, q],
            'init':{rho:(Riemann_pb, (rho_L, rho_R))},
        },
        {
            'velocities':[1,2],
            'conserved_moments':q,
            'polynomials':[1, LA*X],
            'relaxation_parameters':[0, s_q],
            'equilibrium':[q, (gamma-1.)*E+0.5*(3.-gamma)*q**2/rho],
            'init':{q:(Riemann_pb, (q_L, q_R))},
        },
        {
            'velocities':[1,2],
            'conserved_moments':E,
            'polynomials':[1, LA*X],
            'relaxation_parameters':[0, s_E],
            'equilibrium':[E, gamma*E*q/rho-0.5*(gamma-1.)*q**3/rho**2],
            'init':{E:(Riemann_pb, (E_L, E_R))},
        },
    ],
    'boundary_conditions':{
        0:{
            'method':{
                0: pyLBM.bc.neumann,
                1: pyLBM.bc.neumann,
                2: pyLBM.bc.neumann
            },
            'value':None
        },
    },
    'parameters':{LA:la},
    'generator': pyLBM.generator.CythonGenerator,
}

sol = pyLBM.Simulation(dico)
x = sol.domain.x[0][1:-1]
rho = sol.m[0][0][1:-1]
q = sol.m[1][0][1:-1]
E = sol.m[2][0][1:-1]
u = q/rho
p = (gamma-1.)*(E - .5*rho*u**2)
e = E/rho - .5*u**2
viewer = pyLBM.viewer.matplotlibViewer
fig = viewer.Fig(2,3)
lrho = fig[0,0].plot(x, rho, width=2, color='k')[0]
fig[0,0].axis(xmin, xmax, 0., 1.2)
fig[0,0].title = 'mass'
lu = fig[0,1].plot(x, u, width=2, color='k')[0]
fig[0,1].axis(xmin, xmax, 0., 1.)
fig[0,1].title = 'velocity'
lp = fig[0,2].plot(x, p, width=2, color='k')[0]
fig[0,2].axis(xmin, xmax, 0., 1.2)
fig[0,2].title = 'pressure'
lE = fig[1,0].plot(x, E, width=2, color='k')[0]
fig[1,0].axis(xmin, xmax, 0., 3.)
fig[1,0].title = 'energy'
lq = fig[1,1].plot(x, q, width=2, color='k')[0]
fig[1,1].axis(xmin, xmax, 0., .5)
fig[1,1].title = 'momentum'
le = fig[1,2].plot(x, e, width=2, color='k')[0]
fig[1,2].axis(xmin, xmax, 1., 3.)
fig[1,2].title = 'internal energy'

def update(iframe):
    if (sol.t<Tf):
        sol.one_time_step()
        sol.f2m()
        u = q/rho
        p = (gamma-1.)*(E - .5*rho*u**2)
        e = E/rho - .5*u**2
        lrho.set_data(x, rho)
        lu.set_data(x, u)
        lp.set_data(x, p)
        lE.set_data(x, E)
        lq.set_data(x, q)
        le.set_data(x, e)

fig.animate(update, interval = 1)
fig.show()
