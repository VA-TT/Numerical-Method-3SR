"""Compare Python/C++ MPM with the first axial vibration mode."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


L = 25.0
E = 100.0
RHO = 1.0
V0 = 0.1
NELEMENTS = 13
DX = L / NELEMENTS
C = np.sqrt(E / RHO)
DT = 0.1 * DX / C
DURATION = 100.0
BETA1 = np.pi / (2.0 * L)
OMEGA1 = BETA1 * C
TOL = 1.0e-12


def axial_vibration_bar1d(times, x):
    velocity = V0 * np.cos(OMEGA1 * times) * np.sin(BETA1 * x)
    displacement = V0 / OMEGA1 * np.sin(OMEGA1 * times) * np.sin(BETA1 * x)
    return velocity, displacement


def python_mpm():
    """One material point per element, matching mpm1Db.cpp."""
    nodes = np.linspace(0.0, L, NELEMENTS + 1)
    # One Gauss point in a linear element is its midpoint.
    x_initial = 0.5 * (nodes[:-1] + nodes[1:])
    x_particle = x_initial.copy()
    volume = np.full(NELEMENTS, DX)
    mass = RHO * volume
    stress = np.zeros(NELEMENTS)
    velocity = V0 * np.sin(BETA1 * x_particle)

    tracked = NELEMENTS // 2
    nsteps = int(DURATION / DT)
    times = np.arange(nsteps + 1, dtype=float) * DT
    history = np.empty((nsteps + 1, 3))
    history[0] = (0.0, x_particle[tracked], velocity[tracked])

    for step in range(nsteps):
        nodal_mass = np.zeros(NELEMENTS + 1)
        nodal_momentum = np.zeros(NELEMENTS + 1)
        nodal_force = np.zeros(NELEMENTS + 1)

        # Particle-to-grid mapping.
        element_ids = np.clip((x_particle / DX).astype(int), 0, NELEMENTS - 1)
        for particle, element in enumerate(element_ids):
            left, right = element, element + 1
            n_left = (nodes[right] - x_particle[particle]) / DX
            n_right = (x_particle[particle] - nodes[left]) / DX
            nodal_mass[left] += n_left * mass[particle]
            nodal_mass[right] += n_right * mass[particle]
            nodal_momentum[left] += n_left * mass[particle] * velocity[particle]
            nodal_momentum[right] += n_right * mass[particle] * velocity[particle]
            nodal_force[left] += volume[particle] * stress[particle] / DX
            nodal_force[right] -= volume[particle] * stress[particle] / DX

        nodal_momentum[0] = 0.0
        nodal_force[0] = 0.0
        nodal_momentum += nodal_force * DT
        active = nodal_mass > TOL
        nodal_acceleration = np.zeros_like(nodal_mass)
        nodal_velocity = np.zeros_like(nodal_mass)
        nodal_acceleration[active] = nodal_force[active] / nodal_mass[active]
        nodal_velocity[active] = nodal_momentum[active] / nodal_mass[active]
        nodal_acceleration[0] = 0.0
        nodal_velocity[0] = 0.0

        # FLIP particle velocity and grid-velocity position update.
        for particle, element in enumerate(element_ids):
            left, right = element, element + 1
            n_left = (nodes[right] - x_particle[particle]) / DX
            n_right = (x_particle[particle] - nodes[left]) / DX
            velocity[particle] += DT * (
                n_left * nodal_acceleration[left]
                + n_right * nodal_acceleration[right]
            )
            x_particle[particle] += DT * (
                n_left * nodal_velocity[left] + n_right * nodal_velocity[right]
            )

        # MUSL remapping used by the C++ implementation for strain/stress.
        remapped_mass = np.zeros_like(nodal_mass)
        remapped_momentum = np.zeros_like(nodal_momentum)
        element_ids = np.clip((x_particle / DX).astype(int), 0, NELEMENTS - 1)
        for particle, element in enumerate(element_ids):
            left, right = element, element + 1
            n_left = (nodes[right] - x_particle[particle]) / DX
            n_right = (x_particle[particle] - nodes[left]) / DX
            remapped_mass[left] += n_left * mass[particle]
            remapped_mass[right] += n_right * mass[particle]
            remapped_momentum[left] += n_left * mass[particle] * velocity[particle]
            remapped_momentum[right] += n_right * mass[particle] * velocity[particle]
        remapped_velocity = np.zeros_like(remapped_mass)
        active = remapped_mass > TOL
        remapped_velocity[active] = remapped_momentum[active] / remapped_mass[active]
        remapped_velocity[0] = 0.0

        for particle, element in enumerate(element_ids):
            left, right = element, element + 1
            strain_increment = (
                remapped_velocity[right] - remapped_velocity[left]
            ) / DX * DT
            volume[particle] *= 1.0 + strain_increment
            stress[particle] += E * strain_increment

        history[step + 1] = (
            times[step + 1],
            x_particle[tracked],
            velocity[tracked],
        )

    return history, x_initial[tracked]


def main():
    here = Path(__file__).resolve().parent
    cpp_path = here / "mpm1Db_history.txt"
    if not cpp_path.exists():
        raise FileNotFoundError("Run mpm1Db.cpp before plotting.")

    cpp = np.loadtxt(cpp_path, comments="#")
    python, tracked_x0 = python_mpm()
    times = python[:, 0]
    analytical_velocity, analytical_displacement = axial_vibration_bar1d(
        times, tracked_x0
    )

    cpp_displacement = cpp[:, 1] - cpp[0, 1]
    python_displacement = python[:, 1] - python[0, 1]

    print("Maximum absolute error against analytical solution")
    print(
        f"  Python velocity:     "
        f"{np.max(np.abs(python[:, 2] - analytical_velocity)):.6e}"
    )
    print(
        f"  C++ velocity:        "
        f"{np.max(np.abs(cpp[:, 2] - analytical_velocity)):.6e}"
    )
    print(
        f"  Python displacement: "
        f"{np.max(np.abs(python_displacement - analytical_displacement)):.6e}"
    )
    print(
        f"  C++ displacement:    "
        f"{np.max(np.abs(cpp_displacement - analytical_displacement)):.6e}"
    )

    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    axes[0].plot(times, analytical_velocity, "k-", lw=2, label="Analytical")
    axes[0].plot(cpp[:, 0], cpp[:, 2], "C0--", label="C++")
    axes[0].plot(times, python[:, 2], "C1:", lw=2, label="Python")
    axes[0].set_ylabel("Velocity (m/s)")

    axes[1].plot(times, analytical_displacement, "k-", lw=2, label="Analytical")
    axes[1].plot(cpp[:, 0], cpp_displacement, "C0--", label="C++")
    axes[1].plot(times, python_displacement, "C1:", lw=2, label="Python")
    axes[1].set_ylabel("Displacement (m)")
    axes[1].set_xlabel("Time (s)")

    for axis in axes:
        axis.grid(True, alpha=0.3)
        axis.legend(ncol=3)
    fig.tight_layout()
    output = here / "mpm1Db_comparison.png"
    fig.savefig(output, dpi=180)
    print(f"Saved plot to {output}")


if __name__ == "__main__":
    main()
