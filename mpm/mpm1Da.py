"""Compare the one-particle 1-D MPM implementations with the exact solution."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


L = 1.0
E = 4.0 * np.pi**2
RHO = 1.0
V0 = 0.1
X0 = 0.5
DT = 0.01
DURATION = 10.0
MASS = RHO * L
VOLUME0 = L


def analytical_vibration(times):
    omega = np.sqrt(E / RHO) / L
    velocity = V0 * np.cos(omega * times)
    position = X0 * np.exp(V0 / (L * omega) * np.sin(omega * times))
    strain = np.log(position / X0)
    kinetic = 0.5 * MASS * velocity**2
    strain_energy = 0.5 * E * strain**2 * VOLUME0
    return position, velocity, kinetic, strain_energy, kinetic + strain_energy


def python_mpm():
    nodes = np.array([0.0, L])
    x_p, vel_p = X0, V0
    mass_p, volume_p = MASS, VOLUME0
    stress_p = strain_p = 0.0
    momentum_p = mass_p * vel_p

    # Include the true initial state, then label every updated state by t + dt.
    times = np.arange(int(DURATION / DT) + 1, dtype=float) * DT
    rows = [[0.0, x_p, vel_p, 0.5 * mass_p * vel_p**2, 0.0]]

    for step in range(len(times) - 1):
        shape = np.array([1.0 - abs(x_p - nodes[0]) / L,
                          1.0 - abs(x_p - nodes[1]) / L])
        grad_shape = np.array([-1.0 / L, 1.0 / L])

        nodal_mass = shape * mass_p
        nodal_momentum = shape * momentum_p
        nodal_momentum[0] = 0.0

        nodal_force = -grad_shape * volume_p * stress_p
        nodal_force[0] = 0.0
        nodal_momentum += nodal_force * DT

        # FLIP velocity and updated-grid particle position.
        vel_p += DT * np.sum(shape * nodal_force / nodal_mass)
        x_p += DT * np.sum(shape * nodal_momentum / nodal_mass)
        momentum_p = mass_p * vel_p

        nodal_velocity = shape * momentum_p / nodal_mass
        nodal_velocity[0] = 0.0
        strain_increment = np.dot(grad_shape, nodal_velocity) * DT
        strain_p += strain_increment
        stress_p += E * strain_increment

        kinetic = 0.5 * mass_p * vel_p**2
        # Reference volume is appropriate for this small-strain formulation.
        strain_energy = 0.5 * stress_p * strain_p * VOLUME0
        rows.append([times[step + 1], x_p, vel_p, kinetic, strain_energy])

    data = np.asarray(rows)
    total = data[:, 3] + data[:, 4]
    return np.column_stack((data, total))


def load_cpp_history(path):
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path.name}. Compile/run mpm1Da.cpp from the mpm directory first."
        )
    return np.loadtxt(path, comments="#")


def error_summary(name, numerical, exact):
    labels = ("position", "velocity", "total energy")
    columns = (1, 2, 5)
    print(f"\n{name} maximum absolute errors")
    for label, column in zip(labels, columns):
        print(f"  {label:12s}: {np.max(np.abs(numerical[:, column] - exact[:, column])):.6e}")


def main():
    here = Path(__file__).resolve().parent
    python = python_mpm()
    cpp = load_cpp_history(here / "mpm1Da_history.txt")

    times = python[:, 0]
    xa, va, kea, sea, tea = analytical_vibration(times)
    exact = np.column_stack((times, xa, va, kea, sea, tea))

    # C++ stores total energy directly; build the same layout as Python.
    cpp_compact = np.column_stack(
        (cpp[:, 0], cpp[:, 1], cpp[:, 2],
         np.zeros_like(cpp[:, 0]), np.zeros_like(cpp[:, 0]), cpp[:, 3])
    )
    error_summary("Python", python, exact)
    error_summary("C++", cpp_compact, exact)
    print(f"\nMax |Python - C++|: x={np.max(np.abs(python[:, 1] - cpp[:, 1])):.6e}, "
          f"v={np.max(np.abs(python[:, 2] - cpp[:, 2])):.6e}, "
          f"E={np.max(np.abs(python[:, 5] - cpp[:, 3])):.6e}")

    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    axes[0].plot(times, va, "k-", lw=2, label="Analytical")
    axes[0].plot(cpp[:, 0], cpp[:, 2], "C0--", label="C++")
    axes[0].plot(times, python[:, 2], "C1:", lw=2, label="Python")
    axes[0].set_ylabel("Velocity (m/s)")

    # Plot displacement u=x-x0, while retaining position in the data/error report.
    axes[1].plot(times, xa - X0, "k-", lw=2, label="Analytical")
    axes[1].plot(cpp[:, 0], cpp[:, 1] - X0, "C0--", label="C++")
    axes[1].plot(times, python[:, 1] - X0, "C1:", lw=2, label="Python")
    axes[1].set_ylabel("Displacement (m)")

    axes[2].plot(times, tea, "k-", lw=2, label="Analytical total")
    axes[2].plot(cpp[:, 0], cpp[:, 3], "C0--", label="C++ total")
    axes[2].plot(times, python[:, 5], "C1:", lw=2, label="Python total")
    axes[2].set_ylabel("Energy (J)")
    axes[2].set_xlabel("Time (s)")

    for axis in axes:
        axis.grid(True, alpha=0.3)
        axis.legend(ncol=3, fontsize=8)
    fig.tight_layout()
    output = here / "mpm1Da_comparison.png"
    fig.savefig(output, dpi=180)
    print(f"\nSaved plot to {output}")


if __name__ == "__main__":
    main()
