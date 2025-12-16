#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt

# 0. Constitutive law
def Constitutivelaw(law, l, l0, alpha):
    if (law == 1):
        A = alpha * (l - l0) / (l0 * l)
        dA = alpha / l**2
    elif (law == 2):
        A = alpha * (l**2 - l0**2) / (2 * l0**2 * l)
        dA = alpha / (2 * l0**2) + alpha / (2 * l**2)
    else:
        A = alpha * np.log(l / l0) / l
        dA = alpha * (1 - np.log(l / l0)) / l**2
    return A, dA

# 1. Geometry and material properties
# 1.1. Given parameters
nnodes = 6  # Number of nodes
a = 1
x0 = np.array([[0, 0,                 # node 1
               a, 0,                 # node 2
               0, a,                 # node 3
               a, a,                 # node 4
               0, 2 * a,             # node 5
               a, 2 * a]]).transpose()  # node 6

IPN = [1, 2]  # Imposed position nodes (Constraint)

nbars = 8  # Number of bars
no = [1, 3, 2, 3, 3, 4, 4, 5]  # Origin nodes
ne = [3, 2, 4, 4, 5, 5, 6, 6]  # End nodes

Forces = np.array([[0, 0,               # node 1
                   0, 0,               # node 2
                   0, 0,               # node 3
                   0, 0,               # node 4
                   0, 0,               # node 5
                   0, -5000000]]).transpose()  # node 6

E = 200 * 10**9  # Young Modulus
S = 0.0004  # Section
alpha = E * S  # Stiffness

tol = 10**-7  # Tolerance

epsilon = 10**-25  # Used in penalization method

# 1.2. Calculated parameters
L0 = np.zeros(nbars)  # Neutral length L0 of each bar in row
for i in range(nbars):
    a = x0[2 * ne[i] - 2, 0] - x0[2 * no[i] - 2, 0]
    b = x0[2 * ne[i] - 1, 0] - x0[2 * no[i] - 1, 0]
    L0[i] = np.sqrt(a**2 + b**2)

# 2. Solution
# Initialize some variables
count = 0
normU = 2 * tol
x = x0.copy().astype(float)
ID = np.identity(2 * nnodes)
rb = np.zeros(2)
L = np.zeros(nbars)
Total_U = 0
displacement_norms = []  # List to store the norm of displacement at each iteration

while (normU > tol) and (count < 10000):
    count += 1
    K = np.zeros((2 * nnodes, 2 * nnodes))  # Global stiffness matrix
    Nr = Forces.copy()  # Right-hand side vector

    '''2.1. Global stiffness matrix'''
    for i in range(nbars):
        j = 2 * ne[i] - 2  # Index for the end node
        k = 2 * no[i] - 2  # Index for the origin node

        rb[0] = x[j, 0] - x[k, 0]  # Relative x-displacement
        rb[1] = x[j + 1, 0] - x[k + 1, 0]  # Relative y-displacement
        L[i] = np.linalg.norm(rb)  # Current length of the bar

        # Calculate A and dA using the chosen constitutive law
        [A, dA] = Constitutivelaw(1, L[i], L0[i], alpha)

        CE = np.array([ID[j], ID[j + 1]])  # C matrix for the end node
        CO = np.array([ID[k], ID[k + 1]])  # C matrix for the origin node

        # Local stiffness matrix for the bar
        Kb = np.multiply(A, np.identity(2)) + np.multiply(dA / L[i], np.outer(rb, rb))

        # Add contribution to the global stiffness matrix
        K += np.dot(np.subtract(CE, CO).transpose(), np.dot(Kb, np.subtract(CE, CO)))

        '''2.2. Right term'''
        Nr[j, 0]  += -np.multiply(A, rb[0])
        Nr[j + 1, 0] += -np.multiply(A, rb[1])
        Nr[k, 0]  += np.multiply(A, rb[0])
        Nr[k + 1, 0] += np.multiply(A, rb[1])

    '''2.3. Penalization method'''
    for i in range(np.size(IPN)):
        j = 2 * IPN[i] - 2
        K[j, j] += 1 / epsilon
        K[j + 1, j + 1] += 1 / epsilon

    # Solve for displacements
    U = np.linalg.solve(K, Nr)
    Total_U += U
    normU = np.linalg.norm(U)  # Compute norm of displacement vector
    displacement_norms.append(normU)  # Track norm at each iteration
    x += U  # Update node positions

# Print results
print('The displacement of nodes after deformation are:', Total_U)
print('The position of nodes after deformation:', x)
print('The number of iterations:', count)

'''3. Plotting Convergence and Results'''
# Plot convergence: Displacement norm vs Iteration
iterations = list(range(count))
plt.figure(figsize=(8, 6))
plt.plot(iterations, displacement_norms, marker='o', label="Displacement Norm")
plt.axhline(tol, color='r', linestyle='--', label="Tolerance")
plt.xlabel("Iteration")
plt.ylabel("Displacement Norm")
plt.title("Convergence of Displacement Norm")
plt.legend()
plt.grid(True)
plt.show()

# 4. Function to plot truss structure
def plot_truss(x0, x, no, ne, title="Truss Structure"):
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plotting initial bars and nodes
    for i in range(len(no)):
        # Coordinates of the start and end nodes for the initial bar
        x0_start = x0[2 * no[i] - 2:2 * no[i]].flatten()
        x0_end = x0[2 * ne[i] - 2:2 * ne[i]].flatten()

        # Draw the initial bar with a dashed line
        ax.plot([x0_start[0], x0_end[0]], [x0_start[1], x0_end[1]], 'b--', label="Initial" if i == 0 else "")

    # Plotting the bars and nodes after deformation
    for i in range(len(no)):
        # Coordinates of the start and end nodes for the bar after deformation
        x_start = x[2 * no[i] - 2:2 * no[i]].flatten()
        x_end = x[2 * ne[i] - 2:2 * ne[i]].flatten()

        # Draw the deformed bar with a solid line
        ax.plot([x_start[0], x_end[0]], [x_start[1], x_end[1]], 'r-', label="Deformed" if i == 0 else "")

    # Plot the initial and deformed nodes
    ax.plot(x0[0::2], x0[1::2], 'bo', label="Nodes Initial")
    ax.plot(x[0::2], x[1::2], 'ro', label="Nodes Deformed")

    # Set up the plot
    ax.set_title(title)
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.legend()
    ax.grid(True)
    plt.axis('equal')
    plt.show()

# Call the function to plot the truss structure
plot_truss(x0, x, no, ne, "Truss Structure Before and After Deformation")


# In[ ]:




