import numpy as np
import matplotlib.pyplot as plt

# particles per cell
ppc      = 2                 

# mass tolerance
tol = 1e-12

# Domain
L = 25

# Material properties
E = 100
rho = 1

# Computational grid

nelements = 13 # number of elements
dx = L / nelements # element length

# Create equally spaced nodes
x_n = np.linspace(0, L, nelements+1)
nnodes = len(x_n)

# Set-up a 2D array of elements with node ids
elements = np.zeros((nelements, 2), dtype = int)
for nid in range(nelements):
    elements[nid, :] = np.array([nid, nid+1])

# Loading conditions
v0 = 0.1             # initial velocity
c  = np.sqrt(E/rho)  # speed of sound
b1 = np.pi / (2 * L) # beta1
w1 = b1 * c          # omega1

# Create material points 
nparticles = nelements * ppc  # number of particles
# Id of the particle in the central element
pmid = int(np.floor((nparticles/2)))

# Material point properties
x_p      = np.zeros(nparticles)       # positions
vol_p    = np.ones(nparticles)*dx/ppc # volume
mass_p   = vol_p * rho                # mass
stress_p = np.zeros(nparticles)       # stress
vel_p    = np.zeros(nparticles)       # velocity
vol0_p   = vol_p                      # initial volume
defg_p   = np.ones(nparticles)

# Generate particles
x_p, mpoints = generate_particles(ppc, elements, x_n)

for i in range(nparticles):
    # set initial velocities
    vel_p[i] = v0 * np.sin(b1 * x_p[i])

# Time steps and duration
duration = 100
dt_crit = dx / c
dt = 0.1 * dt_crit
t = 0
nsteps = int(duration / dt)

tt, vt, xt = [], [], []

for step in range(nsteps):
    # reset nodal values
    mass_n  = np.zeros(nnodes)  # mass
    mom_n   = np.zeros(nnodes)  # momentum
    fint_n = np.zeros(nnodes)  # internal force

    # iterate through each element
    for eid in range(nelements):
        # get nodal ids
        nid1, nid2 = elements[eid]
        # get particle ids associated with the element 
        mpts = mpoints[eid]
        # iterate through all particles in the element
        for pid in mpts:
            # compute shape functions and derivatives
            N1 = 1 - abs(x_p[pid] - x_n[nid1]) / dx
            N2 = 1 - abs(x_p[pid] - x_n[nid2]) / dx
            dN1 = -1/dx
            dN2 = 1/dx

            # map particle mass and momentum to nodes
            mass_n[nid1] += N1 * mass_p[pid]
            mass_n[nid2] += N2 * mass_p[pid]
            mom_n[nid1]  += N1 * mass_p[pid] * vel_p[pid]
            mom_n[nid2]  += N2 * mass_p[pid] * vel_p[pid]

            # compute nodal internal force
            fint_n[nid1] -= vol_p[pid] * stress_p[pid] * dN1
            fint_n[nid2] -= vol_p[pid] * stress_p[pid] * dN2

    # apply boundary conditions
    mom_n[0]  = 0  # Nodal velocity v = 0 in m * v at node 0.
    fint_n[0] = 0  # Nodal force f = m * a, where a = 0 at node 0.

    # update nodal momentum
    for nid in range(nnodes):
        mom_n[nid] += fint_n[nid] * dt

    # update particle velocity position and stress
    # iterate through each element
    for eid in range(nelements):
        # get nodal ids
        nid1, nid2 = elements[eid]
        # get particle ids associated with the element 
        mpts = mpoints[eid]
        # iterate through all particles in the element
        for pid in mpts:
            # compute shape functions and derivatives
            N1 = 1 - abs(x_p[pid] - x_n[nid1]) / dx
            N2 = 1 - abs(x_p[pid] - x_n[nid2]) / dx
            dN1 = -1/dx
            dN2 = 1/dx

            # compute particle velocity
            if (mass_n[nid1]) > tol:
                vel_p[pid] += dt * N1 * fint_n[nid1] / mass_n[nid1]
            if (mass_n[nid2]) > tol:
                vel_p[pid] += dt * N2 * fint_n[nid2] / mass_n[nid2]
            
            # update particle position based on nodal momentum
            x_p[pid] += dt * (N1 * mom_n[nid1]/mass_n[nid1] + N2 * mom_n[nid2]/mass_n[nid2])

            # nodal velocity
            nv1 = mom_n[nid1]/mass_n[nid1]
            nv2 = mom_n[nid2]/mass_n[nid2]

            # Apply boundary condition
            # Rendundant, since momentum and forces are already set to zero
            # if (nid1 == 0): nv1 = 0

            # rate of strain increment
            grad_v = dN1 * nv1 + dN2 * nv2
            # particle dstrain
            dstrain = grad_v * dt
            # particle volume
            vol_p[pid] *= (1 + dstrain)        
            # update stress using linear elastic model
            stress_p[pid] += E * dstrain

    # update plot params
    tt.append(t)
    vt.append(vel_p[pmid])
    xt.append(x_p[pmid])

    t = t + dt