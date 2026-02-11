#!/usr/bin/env python3
"""
MPM 1D Python Implementation
Implements the same Material Point Method as the C++ version
for direct comparison and validation
"""

import numpy as np
import time

class MPM1D:
    def __init__(self, E, rho, L, v0, dt, duration, xloc, n_nodes=2, n_mp_per_ele=1):
        self.E = E
        self.rho = rho
        self.L = L
        self.v0 = v0
        self.dt = dt
        self.duration = duration
        self.xloc = xloc
        self.n_nodes = n_nodes
        self.n_mp_per_ele = n_mp_per_ele
        
        # Check CFL condition
        c = np.sqrt(E / rho)
        dt_crit = L / c
        assert dt <= dt_crit / 10.0, f"Time step too large! dt={dt}, dt_crit/10={dt_crit/10.0}"
        
        self.n_steps = int(duration / dt)
        self.n_elements = n_nodes - 1
        self.n_mps = self.n_elements * n_mp_per_ele
        
        # Node positions
        self.node_x = np.linspace(0, L, n_nodes)
        
        # Initialize MPs
        self.mass_p = rho * L / self.n_mps
        self.volume_p = L / self.n_mps
        
        # MP arrays
        self.x_p = np.zeros(self.n_mps)
        self.v_p = np.full(self.n_mps, v0)
        self.stress_p = np.zeros(self.n_mps)
        self.strain_p = np.zeros(self.n_mps)
        
        # Nodal arrays
        self.mass_n = np.zeros(n_nodes)
        self.momentum_n = np.zeros(n_nodes)
        self.velocity_n = np.zeros(n_nodes)
        self.force_n = np.zeros(n_nodes)
        
        # Initialize MP positions
        for e in range(self.n_elements):
            x1, x2 = self.node_x[e], self.node_x[e+1]
            for i in range(n_mp_per_ele):
                mp_idx = e * n_mp_per_ele + i
                # Place MP at element center (for 1 MP per element)
                self.x_p[mp_idx] = (x1 + x2) / 2.0
    
    def shape_function(self, x_p, x1, x2):
        """Linear shape functions"""
        L_e = x2 - x1
        N1 = 1.0 - abs(x_p - x1) / L_e
        N2 = 1.0 - abs(x_p - x2) / L_e
        return N1, N2
    
    def shape_function_derivative(self, x1, x2):
        """Shape function derivatives"""
        L_e = x2 - x1
        dN1_dx = -1.0 / L_e
        dN2_dx = 1.0 / L_e
        return dN1_dx, dN2_dx
    
    def find_element(self, x_p):
        """Find which element contains the material point"""
        for e in range(self.n_elements):
            if self.node_x[e] <= x_p <= self.node_x[e+1]:
                return e
        return -1
    
    def p2n(self):
        """Particle to node: map mass and momentum"""
        self.mass_n.fill(0.0)
        self.momentum_n.fill(0.0)
        
        for p in range(self.n_mps):
            x_p = self.x_p[p]
            e = self.find_element(x_p)
            
            if e >= 0:
                x1, x2 = self.node_x[e], self.node_x[e+1]
                N1, N2 = self.shape_function(x_p, x1, x2)
                
                momentum_p = self.mass_p * self.v_p[p]
                
                self.mass_n[e] += N1 * self.mass_p
                self.mass_n[e+1] += N2 * self.mass_p
                self.momentum_n[e] += N1 * momentum_p
                self.momentum_n[e+1] += N2 * momentum_p
        
        # Apply boundary condition: node 0 has zero momentum
        self.momentum_n[0] = 0.0
    
    def nodal_equilibrium(self):
        """Calculate nodal forces"""
        self.force_n.fill(0.0)
        
        for p in range(self.n_mps):
            x_p = self.x_p[p]
            e = self.find_element(x_p)
            
            if e >= 0:
                x1, x2 = self.node_x[e], self.node_x[e+1]
                dN1_dx, dN2_dx = self.shape_function_derivative(x1, x2)
                
                # Internal force
                self.force_n[e] -= self.volume_p * dN1_dx * self.stress_p[p]
                self.force_n[e+1] -= self.volume_p * dN2_dx * self.stress_p[p]
        
        # Apply boundary condition: node 0 has zero force
        self.force_n[0] = 0.0
    
    def n2p(self):
        """Node to particle: update particle state"""
        # Update nodal momentum
        for i in range(self.n_nodes):
            if self.mass_n[i] > 1e-14:
                self.momentum_n[i] += self.force_n[i] * self.dt
        
        # Update particle velocity and position
        for p in range(self.n_mps):
            x_p = self.x_p[p]
            e = self.find_element(x_p)
            
            if e >= 0:
                x1, x2 = self.node_x[e], self.node_x[e+1]
                N1, N2 = self.shape_function(x_p, x1, x2)
                
                # Update velocity (PIC)
                if self.mass_n[e] > 1e-14:
                    acc1 = self.force_n[e] / self.mass_n[e]
                else:
                    acc1 = 0.0
                
                if self.mass_n[e+1] > 1e-14:
                    acc2 = self.force_n[e+1] / self.mass_n[e+1]
                else:
                    acc2 = 0.0
                
                self.v_p[p] += self.dt * (N1 * acc1 + N2 * acc2)
                
                # Update position
                if self.mass_n[e] > 1e-14 and self.mass_n[e+1] > 1e-14:
                    v_grid = N1 * self.momentum_n[e] / self.mass_n[e] + \
                             N2 * self.momentum_n[e+1] / self.mass_n[e+1]
                    self.x_p[p] += v_grid * self.dt
        
        # Update nodal velocity for strain calculation
        for p in range(self.n_mps):
            x_p = self.x_p[p]
            e = self.find_element(x_p)
            
            if e >= 0:
                x1, x2 = self.node_x[e], self.node_x[e+1]
                N1, N2 = self.shape_function(x_p, x1, x2)
                
                if self.mass_n[e] > 1e-14:
                    self.velocity_n[e] += self.mass_p * self.v_p[p] * N1 / self.mass_n[e]
                if self.mass_n[e+1] > 1e-14:
                    self.velocity_n[e+1] += self.mass_p * self.v_p[p] * N2 / self.mass_n[e+1]
        
        # Apply velocity BC
        self.velocity_n[0] = 0.0
        
        # Update strain and stress
        for p in range(self.n_mps):
            x_p = self.x_p[p]
            e = self.find_element(x_p)
            
            if e >= 0:
                x1, x2 = self.node_x[e], self.node_x[e+1]
                dN1_dx, dN2_dx = self.shape_function_derivative(x1, x2)
                
                strain_rate = dN1_dx * self.velocity_n[e] + dN2_dx * self.velocity_n[e+1]
                d_strain = strain_rate * self.dt
                self.strain_p[p] += d_strain
                
                # Linear elastic constitutive law
                self.stress_p[p] += self.E * d_strain
    
    def reset_grid(self):
        """Reset grid quantities for next step"""
        self.mass_n.fill(0.0)
        self.momentum_n.fill(0.0)
        self.velocity_n.fill(0.0)
        self.force_n.fill(0.0)
    
    def run(self):
        """Run the simulation"""
        history = []
        
        # Record initial state
        x_num = self.x_p[0]
        v_num = self.v_p[0]
        history.append([0.0, x_num, v_num])
        
        start_time = time.time()
        
        for step in range(self.n_steps):
            self.p2n()
            self.nodal_equilibrium()
            self.n2p()
            self.reset_grid()
            
            # Record state
            t = (step + 1) * self.dt
            x_num = self.x_p[0]
            v_num = self.v_p[0]
            history.append([t, x_num, v_num])
            
            if step < 5 or step == self.n_steps - 1:
                print(f"Step {step} t={t:.6f} | x={x_num:.6f} v={v_num:.6f}")
        
        elapsed = time.time() - start_time
        print(f"\nComputation time: {elapsed:.6f} seconds")
        
        return np.array(history)


def analytical_solution(t, L, v0, E, rho, xloc):
    """Analytical solution for vibrating bar"""
    omega = (1.0 / L) * np.sqrt(E / rho)
    v_ana = v0 * np.cos(omega * t)
    x_ana = xloc * np.exp((v0 / (L * omega)) * np.sin(omega * t))
    return x_ana, v_ana


def main():
    print("=" * 70)
    print("MPM 1D Python Implementation")
    print("=" * 70)
    
    # Parameters (matching C++ test)
    L = 1.0
    v0 = 0.1
    E = 4.0 * np.pi**2
    rho = 1.0
    dt = 0.01
    duration = 10.0
    xloc = 0.5
    
    # Create and run simulation
    mpm = MPM1D(E, rho, L, v0, dt, duration, xloc, n_nodes=2, n_mp_per_ele=1)
    history = mpm.run()
    
    # Compute analytical solution
    print("\nComputing analytical solution...")
    t_array = history[:, 0]
    x_ana, v_ana = analytical_solution(t_array, L, v0, E, rho, xloc)
    
    # Save results
    output = np.column_stack([t_array, history[:, 1], history[:, 2], x_ana, v_ana])
    np.savetxt('mpm1D_python.txt', output, 
               header='time\tx_num\tv_num\tx_ana\tv_ana',
               fmt='%.10f', delimiter='\t')
    print(f"Results saved to: mpm1D_python.txt")
    
    # Compute errors
    x_error = np.abs(history[:, 1] - x_ana)
    v_error = np.abs(history[:, 2] - v_ana)
    
    print("\n" + "=" * 70)
    print("Python Results Summary:")
    print("=" * 70)
    print(f"Position max error: {np.max(x_error):.6e}")
    print(f"Position mean error: {np.mean(x_error):.6e}")
    print(f"Velocity max error: {np.max(v_error):.6e}")
    print(f"Velocity mean error: {np.mean(v_error):.6e}")
    print(f"\nFinal state (t={t_array[-1]:.1f}):")
    print(f"  x_num={history[-1, 1]:.6f}, x_ana={x_ana[-1]:.6f}")
    print(f"  v_num={history[-1, 2]:.6f}, v_ana={v_ana[-1]:.6f}")
    print("=" * 70)


if __name__ == "__main__":
    main()
