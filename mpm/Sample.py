import numpy as np
import matplotlib.pyplot as plt

import numpy as np

# Analytical solution
def analytical_vibration(E, rho, v0, x_loc, duration, dt, L):
    nsteps = int(duration/dt)
    tt, vt, xt = [], [], []
    t = 0
    for _ in range(nsteps):
        omega = 1. / L * np.sqrt(E / rho)
        v = v0 * np.cos(omega * t)
        x = x_loc * np.exp(v0 / (L * omega) * np.sin(omega * t))
        vt.append(v)
        xt.append(x)
        tt.append(t)
        t += dt
    return tt, vt, xt

# Computational grid
L          = 1                       # domain size
nodes      = np.array([0, L])        # nodal coordinates
nnodes     = len(nodes)              # number of nodes
nelements  = 1                       # number of elements
nparticles = 1                       # number of particles
el_length  = L / nelements           # element length

# Initial conditions 
v0         = 0.1                     # initial velocity
x_loc      = 0.5                     # location to get analytical solution

# Material property
E          = 4 * (np.pi)**2          # Young's modulus
rho        = 1.                      # Density

# Material points
x_p        = 0.5 * el_length         # position
mass_p     = 1.                      # mass
vol_p      = el_length / nparticles  # volume
vel_p      = v0                      # velocity
stress_p   = 0.                      # stress
strain_p   = 0.                      # strain
mv_p       = mass_p * vel_p          # momentum = m * v

# Time
duration   = 10
dt         = 0.01
time       = 0
nsteps     = int(duration/dt)

# Store time, velocity and position with time
time_t, vel_t, x_t, se_t, ke_t, te_t = [], [], [], [], [], []

for _ in range(nsteps):
    # shape function and derivative
    N  = np.array([1 - abs(x_p - nodes[0])/L, 1 - abs(x_p - nodes[1])/L])
    dN = np.array([-1/L, 1/L])

    # map particle mass and momentum to nodes
    mass_n = N * mass_p
    mv_n   = N * mv_p
    
    # apply boundary condition: velocity at left node is zero
    mv_n[0] = 0

    # external force at nodes
    f_ext_n = np.array([0, 0])

    # compute internal force at nodes
    f_int_n = - dN * vol_p * stress_p 

    # total forces at nodes
    f_total_n = f_int_n + f_ext_n

    # apply boundary condition: left node has no acceleration (f = m * a, and a = 0)
    f_total_n[0] = 0

    # update nodal momentum
    mv_n += f_total_n * dt

    # update particle position and velocity
    for i in range(nnodes):
        vel_p += dt * N[i] * f_total_n[i] / mass_n[i]
        x_p += dt * N[i] * mv_n[i] / mass_n[i]

    # update particle momentum
    mv_p = mass_p * vel_p

    # map nodal velocity
    vel_n = mass_p * vel_p * np.divide(N, mass_n)
    # Apply boundary condition and set left nodal velocity to zero
    vel_n[0] = 0

    # compute strain rate at the particle
    strain_rate_p = np.dot(dN, vel_n) 
    # compute strain increament 
    dstrain_p = strain_rate_p * dt
    # compute strain
    strain_p += dstrain_p
    # compute stress
    stress_p += E * dstrain_p

    # store properties for plotting
    time_t.append(time)
    vel_t.append(vel_p)
    x_t.append(x_p)

    # Energies
    strain_energy = 0.5 * stress_p * strain_p * vol_p
    kinetic_energy = 0.5 * vel_p**2 * mass_p**2
    total_energy = strain_energy + kinetic_energy
    se_t.append(strain_energy)
    ke_t.append(kinetic_energy)
    te_t.append(total_energy)
    
    # update time
    time += dt


ta, va, xa = analytical_vibration(E, rho, v0, x_loc, duration, dt, L)

plt.plot(ta, va, 'r', linewidth=2,label='analytical')
plt.plot(time_t, vel_t, 'ob', markersize = 2, label='mpm')

plt.xlabel('time (s)')
plt.ylabel('velocity (m/s)')
plt.legend()
plt.show()