from pathlib import Path

import numpy as np 
import matplotlib.pyplot as plt
from tqdm import tqdm

def elasticity_matrix(youngs_modulus: float, poissons_ratio: float, stress_state: str) -> np.ndarray:
    """Compute the elasticity (stiffness) matrix for an isotropic elastic material.
    
    Args:
        youngs_modulus: Young's modulus of the material.
        poissons_ratio: Poisson's ratio of the material.
        stress_state: Type of stress state. Must be one of "PLANE_STRESS", 
                      "PLANE_STRAIN", or "3D".
    
    Returns:
        Elasticity matrix as a numpy array.
        
    Raises:
        ValueError: If an invalid stress_state is provided.
    """
    if stress_state == 'PLANE_STRESS':
        C = youngs_modulus / (1 - poissons_ratio**2) * np.array([
            [1, poissons_ratio, 0],
            [poissons_ratio, 1, 0],
            [0, 0, (1 - poissons_ratio) / 2]
        ])
    elif stress_state == 'PLANE_STRAIN':
        C = youngs_modulus / ((1 + poissons_ratio) * (1 - 2 * poissons_ratio)) * np.array([
            [1 - poissons_ratio, poissons_ratio, 0],
            [poissons_ratio, 1 - poissons_ratio, 0],
            [0, 0, 0.5 - poissons_ratio]
        ])
    else:  # 3D stress state
        C = np.zeros((6, 6))
        C[:3, :3] = youngs_modulus / ((1 + poissons_ratio) * (1 - 2 * poissons_ratio)) * np.array([
            [1 - poissons_ratio, poissons_ratio, poissons_ratio],
            [poissons_ratio, 1 - poissons_ratio, poissons_ratio],
            [poissons_ratio, poissons_ratio, 1 - poissons_ratio]
        ])
        C[3:, 3:] = youngs_modulus / (2 * (1 + poissons_ratio)) * np.eye(3)
    return C


def generate_rectangular_particles(origin, length, height, spacing):
    """
    Generate a numpy array that discretizes a rectangular shaped mass into equally spaced material points.
    
    Parameters:
    origin (tuple): (x, y) coordinates of the bottom left corner of the rectangle
    length (float): Length of the rectangle
    width (float): Width of the rectangle
    spacing (float): Distance between adjacent particles
    
    Returns:
    numpy.ndarray: Array of shape (n, 2) containing the coordinates of material points
    """
    
    # Calculate the number of points along each dimension
    nx = int(np.ceil(length / spacing)) + 1
    ny = int(np.ceil(height / spacing)) + 1
    
    # Generate equally spaced points along x and y axes
    x = np.linspace(origin[0], origin[0] + length, nx)
    y = np.linspace(origin[1], origin[1] + height, ny)
    
    # Create a mesh grid
    xx, yy = np.meshgrid(x, y)
    
    # Reshape the meshgrid into a 2D array of points
    points = np.column_stack((xx.ravel(), yy.ravel()))
    
    return points


def particle_element_mapping(particle_positions: np.ndarray, cell_width: float, cell_height: float,
                             n_ele_x: int, n_elements: int):
    """Map each particle to an element based on its position.
    
    Args:
        particle_positions: Particle positions as an array of shape (n_particles, 2).
        cell_width: Element (cell) width in the x-direction.
        cell_height: Element (cell) height in the y-direction.
        n_ele_x: Number of elements in the x-direction.
        n_elements: Total number of elements.
    
    Returns:
        ele_ids: Array of element ids for each particle.
        particle_ids_in_elements: List of particle indices for each element.
    """
    x_indices = np.clip(
        np.floor(particle_positions[:, 0] / cell_width).astype(int), 0, n_ele_x - 1
    )
    n_ele_y = n_elements // n_ele_x
    y_indices = np.clip(
        np.floor(particle_positions[:, 1] / cell_height).astype(int), 0, n_ele_y - 1
    )
    ele_ids = x_indices + n_ele_x * y_indices

    particle_ids_in_elements = [[] for _ in range(n_elements)]
    for p, e in enumerate(ele_ids):
        particle_ids_in_elements[e].append(p)

    return ele_ids, particle_ids_in_elements


def lagrange_basis(element_type: str, local_coord: np.ndarray):
    """Compute the shape functions and their derivatives for a finite element.
    
    Args:
        element_type: Type of element. Currently only "Q4" is supported.
        local_coord: Array of local coordinates [xi, eta].
    
    Returns:
        N: Shape function values.
        dN_dxi: Derivatives of the shape functions with respect to xi and eta.
        
    Raises:
        ValueError: If an unsupported element type is provided.
    """
    if element_type == 'Q4':
        xi, eta = local_coord
        N = np.array([
            (1 - xi) * (1 - eta),
            (1 + xi) * (1 - eta),
            (1 + xi) * (1 + eta),
            (1 - xi) * (1 + eta)
        ]) / 4.0

        dN_dxi = np.array([
            [-(1 - eta), -(1 - xi)],
            [ (1 - eta), -(1 + xi)],
            [ (1 + eta),  (1 + xi)],
            [-(1 + eta),  (1 - xi)]
        ]) / 4.0
    else:
        raise ValueError(f"Unsupported element type: {element_type}")

    return N, dN_dxi


def compute_local_coord(particle_pos: np.ndarray, element_coords: np.ndarray,
                        cell_width: float, cell_height: float) -> np.ndarray:
    """Compute the local (ξ, η) coordinates for a particle within an element.
    
    Args:
        particle_pos: [x, y] coordinates of the particle.
        element_coords: Coordinates of the element nodes (4 x 2 array).
        cell_width: Element width in the x-direction.
        cell_height: Element height in the y-direction.
    
    Returns:
        Local coordinates [ξ, η].
    """
    xi = (2 * particle_pos[0] - (element_coords[0, 0] + element_coords[1, 0])) / cell_width
    eta = (2 * particle_pos[1] - (element_coords[1, 1] + element_coords[2, 1])) / cell_height
    return np.array([xi, eta])

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Temporary hard-coded variable
nu = 0.3

def principal_stresses(stress_tensor):
    """
    Calculate the principal stresses sigma_1 and sigma_2 for a 2D stress tensor using eigen decomposition.

    Parameters:
    stress_tensor (list): A list of the form [sigma_xx, sigma_yy, sigma_xy].

    Returns:
    tuple: Principal stresses (sigma_1, sigma_2) with sigma_1 >= sigma_2.
    """
    # Unpack the stress components
    sigma_xx, sigma_yy, sigma_xy = stress_tensor

    # Form the 2D stress tensor matrix
    stress_matrix = np.array([[sigma_xx, sigma_xy],
                              [sigma_xy, sigma_yy]])

    # Perform eigen decomposition to get the eigenvalues (principal stresses)
    eigenvalues, _ = np.linalg.eig(stress_matrix)

    # Sort eigenvalues in descending order: sigma_1 >= sigma_2
    sigma_1, sigma_2 = np.sort(eigenvalues)[::-1]

    return sigma_1, sigma_2


# Function to compute alpha and k from friction angle and cohesion
def drucker_prager_parameters(friction_angle_deg, cohesion):
    # Convert friction angle to radians
    phi = np.radians(friction_angle_deg)
    
    # Compute alpha and k
    alpha = (2 * np.sin(phi)) / (np.sqrt(3) * (3 - np.sin(phi)))
    # alpha = 0
    k = (6 * cohesion * np.cos(phi)) / (np.sqrt(3) * (3 - np.sin(phi)))
    
    return alpha, k

# Elasticity matrix (2D plane strain assumption)
def elasticity_matrix(E, nu):
    D = (E / ((1 + nu) * (1 - 2 * nu))) * np.array([
        [1 - nu,     nu,          0],
        [nu,         1 - nu,      0],
        [0,          0,           (1 - 2 * nu) / 2]
    ])
    return D

# Invariants of the stress tensor
def compute_invariants(stress, stress_condition="plane_strain"):
    sigma_xx, sigma_yy, tau_xy = stress
    
    if stress_condition != "plane_strain":
        # First invariant
        I1 = sigma_xx + sigma_yy  
        
        # J2 invariant
        # Deviatoric stresses
        s_xx = sigma_xx - I1 / 2
        s_yy = sigma_yy - I1 / 2
        s_xy = tau_xy
        
        # J2 = 0.5 * s:s
        J2 = 0.5 * (s_xx**2 + s_yy**2 + 2 * s_xy**2)
        # The above is equivalent to:
        # J2 = (1/4) * (sigma_x - sigma_y)**2 + tau_xy**2
        
    else:  # Plain strain constraint
        sigma_zz = nu * (sigma_xx + sigma_yy)
        I1 = sigma_xx + sigma_yy + sigma_zz
        
        # Deviatoric stresses
        s_xx = sigma_xx - I1 / 3
        s_yy = sigma_yy - I1 / 3
        s_zz = sigma_zz - I1 / 3
        s_xy = tau_xy
        
        # J2 = 0.5 s:s
        J2 = 0.5 * (s_xx**2 + s_yy**2 + s_zz**2 + 2 * s_xy**2)
    
    return I1, J2

# Drucker-Prager yield function
def drucker_prager_yield(stress, alpha, k):
    I1, J2 = compute_invariants(stress)
    return alpha * I1 + np.sqrt(J2) - k

# Gradient of yield function df/dsigma
def df_dsigma(stress, alpha):
    sigma_x, sigma_y, tau_xy = stress
    I1, J2 = compute_invariants(stress)
    root_j2 = np.sqrt(J2)
    if root_j2 <= 1.0e-14:
        return np.array([alpha * (1 + nu), alpha * (1 + nu), 0.0])

    sigma_z = nu * (sigma_x + sigma_y)
    mean = (sigma_x + sigma_y + sigma_z) / 3.0
    s_x, s_y, s_z = sigma_x - mean, sigma_y - mean, sigma_z - mean
    one_plus_nu = 1.0 + nu
    dJ2dx = (
        s_x * (1.0 - one_plus_nu / 3.0)
        - s_y * one_plus_nu / 3.0
        + s_z * (nu - one_plus_nu / 3.0)
    )
    dJ2dy = (
        -s_x * one_plus_nu / 3.0
        + s_y * (1.0 - one_plus_nu / 3.0)
        + s_z * (nu - one_plus_nu / 3.0)
    )
    dfdx = alpha * one_plus_nu + dJ2dx / (2.0 * root_j2)
    dfdy = alpha * one_plus_nu + dJ2dy / (2.0 * root_j2)
    dftau = tau_xy / root_j2
    return np.array([dfdx, dfdy, dftau])

# Update stress and strain process
def update_stress_strain(stress_n, strain_n, strain_increment, D, alpha, k):
    # Step 1: Compute trial stress
    stress_trial = stress_n + np.dot(D, strain_increment)
    
    # Step 2: Check yield condition
    f_trial = drucker_prager_yield(stress_trial, alpha, k)
    # print(f_trial)
    
    if f_trial <= 0:
        # Elastic step: no plastic deformation
        return stress_trial, strain_n + strain_increment, 0.0  # No plastic multiplier
    
    # Iterative return mapping because the Drucker-Prager normal changes along
    # the correction direction.
    stress_updated = stress_trial.copy()
    delta_lambda = 0.0
    tolerance = 1.0e-10 * (1.0 + abs(k) + np.linalg.norm(stress_trial))
    for _ in range(25):
        f_value = drucker_prager_yield(stress_updated, alpha, k)
        if f_value <= tolerance:
            break
        df_sigma = df_dsigma(stress_updated, alpha)
        D_df = D @ df_sigma
        H = df_sigma @ D_df
        if abs(H) <= 1.0e-14:
            break
        increment = f_value / H
        stress_updated -= increment * D_df
        delta_lambda += increment
    
    # Store total strain. Plastic strain requires a separate state variable.
    strain_updated = strain_n + strain_increment
    
    return stress_updated, strain_updated, delta_lambda

def elasticity_matrix(youngs_modulus: float, poissons_ratio: float, stress_state: str) -> np.ndarray:
    """Compute the elasticity (stiffness) matrix for an isotropic elastic material.
    
    Args:
        youngs_modulus: Young's modulus of the material.
        poissons_ratio: Poisson's ratio of the material.
        stress_state: Type of stress state. Must be one of "PLANE_STRESS", 
                      "PLANE_STRAIN", or "3D".
    
    Returns:
        Elasticity matrix as a numpy array.
        
    Raises:
        ValueError: If an invalid stress_state is provided.
    """
    if stress_state == 'PLANE_STRESS':
        C = youngs_modulus / (1 - poissons_ratio**2) * np.array([
            [1, poissons_ratio, 0],
            [poissons_ratio, 1, 0],
            [0, 0, (1 - poissons_ratio) / 2]
        ])
    elif stress_state == 'PLANE_STRAIN':
        C = youngs_modulus / ((1 + poissons_ratio) * (1 - 2 * poissons_ratio)) * np.array([
            [1 - poissons_ratio, poissons_ratio, 0],
            [poissons_ratio, 1 - poissons_ratio, 0],
            [0, 0, 0.5 - poissons_ratio]
        ])
    else:  # 3D stress state
        C = np.zeros((6, 6))
        C[:3, :3] = youngs_modulus / ((1 + poissons_ratio) * (1 - 2 * poissons_ratio)) * np.array([
            [1 - poissons_ratio, poissons_ratio, poissons_ratio],
            [poissons_ratio, 1 - poissons_ratio, poissons_ratio],
            [poissons_ratio, poissons_ratio, 1 - poissons_ratio]
        ])
        C[3:, 3:] = youngs_modulus / (2 * (1 + poissons_ratio)) * np.eye(3)
    return C


def generate_rectangular_particles(origin, length, height, spacing):
    """
    Generate a numpy array that discretizes a rectangular shaped mass into equally spaced material points.
    
    Parameters:
    origin (tuple): (x, y) coordinates of the bottom left corner of the rectangle
    length (float): Length of the rectangle
    width (float): Width of the rectangle
    spacing (float): Distance between adjacent particles
    
    Returns:
    numpy.ndarray: Array of shape (n, 2) containing the coordinates of material points
    """
    
    # Calculate the number of points along each dimension
    nx = int(np.ceil(length / spacing)) + 1
    ny = int(np.ceil(height / spacing)) + 1
    
    # Generate equally spaced points along x and y axes
    x = np.linspace(origin[0], origin[0] + length, nx)
    y = np.linspace(origin[1], origin[1] + height, ny)
    
    # Create a mesh grid
    xx, yy = np.meshgrid(x, y)
    
    # Reshape the meshgrid into a 2D array of points
    points = np.column_stack((xx.ravel(), yy.ravel()))
    
    return points


def particle_element_mapping(particle_positions: np.ndarray, cell_width: float, cell_height: float,
                             n_ele_x: int, n_elements: int):
    """Map each particle to an element based on its position.
    
    Args:
        particle_positions: Particle positions as an array of shape (n_particles, 2).
        cell_width: Element (cell) width in the x-direction.
        cell_height: Element (cell) height in the y-direction.
        n_ele_x: Number of elements in the x-direction.
        n_elements: Total number of elements.
    
    Returns:
        ele_ids: Array of element ids for each particle.
        particle_ids_in_elements: List of particle indices for each element.
    """
    x_indices = np.clip(
        np.floor(particle_positions[:, 0] / cell_width).astype(int), 0, n_ele_x - 1
    )
    n_ele_y = n_elements // n_ele_x
    y_indices = np.clip(
        np.floor(particle_positions[:, 1] / cell_height).astype(int), 0, n_ele_y - 1
    )
    ele_ids = x_indices + n_ele_x * y_indices

    particle_ids_in_elements = [[] for _ in range(n_elements)]
    for p, e in enumerate(ele_ids):
        particle_ids_in_elements[e].append(p)

    return ele_ids, particle_ids_in_elements


def lagrange_basis(element_type: str, local_coord: np.ndarray):
    """Compute the shape functions and their derivatives for a finite element.
    
    Args:
        element_type: Type of element. Currently only "Q4" is supported.
        local_coord: Array of local coordinates [xi, eta].
    
    Returns:
        N: Shape function values.
        dN_dxi: Derivatives of the shape functions with respect to xi and eta.
        
    Raises:
        ValueError: If an unsupported element type is provided.
    """
    if element_type == 'Q4':
        xi, eta = local_coord
        N = np.array([
            (1 - xi) * (1 - eta),
            (1 + xi) * (1 - eta),
            (1 + xi) * (1 + eta),
            (1 - xi) * (1 + eta)
        ]) / 4.0

        dN_dxi = np.array([
            [-(1 - eta), -(1 - xi)],
            [ (1 - eta), -(1 + xi)],
            [ (1 + eta),  (1 + xi)],
            [-(1 + eta),  (1 - xi)]
        ]) / 4.0
    else:
        raise ValueError(f"Unsupported element type: {element_type}")

    return N, dN_dxi


def compute_local_coord(particle_pos: np.ndarray, element_coords: np.ndarray,
                        cell_width: float, cell_height: float) -> np.ndarray:
    """Compute the local (ξ, η) coordinates for a particle within an element.
    
    Args:
        particle_pos: [x, y] coordinates of the particle.
        element_coords: Coordinates of the element nodes (4 x 2 array).
        cell_width: Element width in the x-direction.
        cell_height: Element height in the y-direction.
    
    Returns:
        Local coordinates [ξ, η].
    """
    xi = (2 * particle_pos[0] - (element_coords[0, 0] + element_coords[1, 0])) / cell_width
    eta = (2 * particle_pos[1] - (element_coords[1, 1] + element_coords[2, 1])) / cell_height
    return np.array([xi, eta])

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Temporary hard-coded variable
nu = 0.3

def principal_stresses(stress_tensor):
    """
    Calculate the principal stresses sigma_1 and sigma_2 for a 2D stress tensor using eigen decomposition.

    Parameters:
    stress_tensor (list): A list of the form [sigma_xx, sigma_yy, sigma_xy].

    Returns:
    tuple: Principal stresses (sigma_1, sigma_2) with sigma_1 >= sigma_2.
    """
    # Unpack the stress components
    sigma_xx, sigma_yy, sigma_xy = stress_tensor

    # Form the 2D stress tensor matrix
    stress_matrix = np.array([[sigma_xx, sigma_xy],
                              [sigma_xy, sigma_yy]])

    # Perform eigen decomposition to get the eigenvalues (principal stresses)
    eigenvalues, _ = np.linalg.eig(stress_matrix)

    # Sort eigenvalues in descending order: sigma_1 >= sigma_2
    sigma_1, sigma_2 = np.sort(eigenvalues)[::-1]

    return sigma_1, sigma_2


# Function to compute alpha and k from friction angle and cohesion
def drucker_prager_parameters(friction_angle_deg, cohesion):
    # Convert friction angle to radians
    phi = np.radians(friction_angle_deg)
    
    # Compute alpha and k
    alpha = (2 * np.sin(phi)) / (np.sqrt(3) * (3 - np.sin(phi)))
    # alpha = 0
    k = (6 * cohesion * np.cos(phi)) / (np.sqrt(3) * (3 - np.sin(phi)))
    
    return alpha, k

# Elasticity matrix (2D plane strain assumption)
def elasticity_matrix(E, nu):
    D = (E / ((1 + nu) * (1 - 2 * nu))) * np.array([
        [1 - nu,     nu,          0],
        [nu,         1 - nu,      0],
        [0,          0,           (1 - 2 * nu) / 2]
    ])
    return D

# Invariants of the stress tensor
def compute_invariants(stress, stress_condition="plane_strain"):
    sigma_xx, sigma_yy, tau_xy = stress
    
    if stress_condition != "plane_strain":
        # First invariant
        I1 = sigma_xx + sigma_yy  
        
        # J2 invariant
        # Deviatoric stresses
        s_xx = sigma_xx - I1 / 2
        s_yy = sigma_yy - I1 / 2
        s_xy = tau_xy
        
        # J2 = 0.5 * s:s
        J2 = 0.5 * (s_xx**2 + s_yy**2 + 2 * s_xy**2)
        # The above is equivalent to:
        # J2 = (1/4) * (sigma_x - sigma_y)**2 + tau_xy**2
        
    else:  # Plain strain constraint
        sigma_zz = nu * (sigma_xx + sigma_yy)
        I1 = sigma_xx + sigma_yy + sigma_zz
        
        # Deviatoric stresses
        s_xx = sigma_xx - I1 / 3
        s_yy = sigma_yy - I1 / 3
        s_zz = sigma_zz - I1 / 3
        s_xy = tau_xy
        
        # J2 = 0.5 s:s
        J2 = 0.5 * (s_xx**2 + s_yy**2 + s_zz**2 + 2 * s_xy**2)
    
    return I1, J2

# Drucker-Prager yield function
def drucker_prager_yield(stress, alpha, k):
    I1, J2 = compute_invariants(stress)
    return alpha * I1 + np.sqrt(J2) - k

# Gradient of yield function df/dsigma
def df_dsigma(stress, alpha):
    sigma_x, sigma_y, tau_xy = stress
    I1, J2 = compute_invariants(stress)
    root_j2 = np.sqrt(J2)
    if root_j2 <= 1.0e-14:
        return np.array([alpha * (1 + nu), alpha * (1 + nu), 0.0])

    sigma_z = nu * (sigma_x + sigma_y)
    mean = (sigma_x + sigma_y + sigma_z) / 3.0
    s_x, s_y, s_z = sigma_x - mean, sigma_y - mean, sigma_z - mean
    one_plus_nu = 1.0 + nu
    dJ2dx = (
        s_x * (1.0 - one_plus_nu / 3.0)
        - s_y * one_plus_nu / 3.0
        + s_z * (nu - one_plus_nu / 3.0)
    )
    dJ2dy = (
        -s_x * one_plus_nu / 3.0
        + s_y * (1.0 - one_plus_nu / 3.0)
        + s_z * (nu - one_plus_nu / 3.0)
    )
    dfdx = alpha * one_plus_nu + dJ2dx / (2.0 * root_j2)
    dfdy = alpha * one_plus_nu + dJ2dy / (2.0 * root_j2)
    dftau = tau_xy / root_j2
    return np.array([dfdx, dfdy, dftau])

# Update stress and strain process
def update_stress_strain(stress_n, strain_n, strain_increment, D, alpha, k):
    # Step 1: Compute trial stress
    stress_trial = stress_n + np.dot(D, strain_increment)
    
    # Step 2: Check yield condition
    f_trial = drucker_prager_yield(stress_trial, alpha, k)
    # print(f_trial)
    
    if f_trial <= 0:
        # Elastic step: no plastic deformation
        return stress_trial, strain_n + strain_increment, 0.0  # No plastic multiplier
    
    # Iterative return mapping because the Drucker-Prager normal changes along
    # the correction direction.
    stress_updated = stress_trial.copy()
    delta_lambda = 0.0
    tolerance = 1.0e-10 * (1.0 + abs(k) + np.linalg.norm(stress_trial))
    for _ in range(25):
        f_value = drucker_prager_yield(stress_updated, alpha, k)
        if f_value <= tolerance:
            break
        df_sigma = df_dsigma(stress_updated, alpha)
        D_df = D @ df_sigma
        H = df_sigma @ D_df
        if abs(H) <= 1.0e-14:
            break
        increment = f_value / H
        stress_updated -= increment * D_df
        delta_lambda += increment
    
    # Store total strain. Plastic strain requires a separate state variable.
    strain_updated = strain_n + strain_increment
    
    return stress_updated, strain_updated, delta_lambda

    # Analysis setting
time = 0.1
t = 0
interval = 200
g = 9.81

# Material parameters matching mpm2Da.cpp.
E = 100000.0
nu = 0.3
rho = 3600.0
mu = 0.385
phi = 30.0
c = 1.0
K0 = 0.5
dtime = 0.0001
pic_ratio = 1.0

# length of domain
l = 1

# number of elements for each dim
n_ele_x = 20
n_ele_y = 20

# cell size
delta_x = l / n_ele_x
delta_y = l / n_ele_y

# number of grid nodes
n_node_x = n_ele_x + 1
n_node_y = n_ele_y + 1
n_nodes = n_node_x * n_node_y

# Generate mesh nodes
x_coords = np.linspace(0, l, n_node_x)
y_coords = np.linspace(0, l, n_node_y)
X, Y = np.meshgrid(x_coords, y_coords)
nodes = np.vstack([X.flatten(), Y.flatten()]).T

# Generate elements (connectivity)
elements = []
for j in range(n_node_y - 1):
    for i in range(n_node_x - 1):
        n1 = j * n_node_x + i
        n2 = n1 + 1
        n3 = n2 + n_node_x
        n4 = n1 + n_node_x
        elements.append([n1, n2, n3, n4])
elements = np.array(elements)
elements_coordinate_map = nodes[elements]
n_elements = elements.shape[0]

# Boundary nodes
bottomNodes = np.where(nodes[:, 1] < 1e-8)[0]
upperNodes = np.where(nodes[:, 1] > l - 1e-8)[0]
leftNodes = np.where(nodes[:, 0] < 1e-8)[0]
rightNodes = np.where(nodes[:, 0] > l - 1e-8)[0]

# Initialize node quantities
n_masses = np.zeros(n_nodes)  # nodal masses
n_momentums = np.zeros((n_nodes, 2))  # nodal momentums
n_iforces = np.zeros((n_nodes, 2))  # Internal forces
n_eforces = np.zeros((n_nodes, 2))  # External forces. Not used in this practice
# Initialize grid state variables
gStress = np.zeros((3, n_nodes))  # [sig_xx, sig_yy, sig_xy]
gDisp = np.zeros((2, n_nodes))  # [disp_x, disp_y]

# Inputs for particles representing granular column
origin = [0.0125, 0.0125]
lengths = [0.300, 0.300]
n_particle_per_cell_per_dim = 2

# Generate 12 x 12 cell-centred particles, matching the C++ particle domain.
spacing = delta_x / n_particle_per_cell_per_dim
particle_x = np.arange(origin[0], lengths[0], spacing)
particle_y = np.arange(origin[1], lengths[1], spacing)
particle_xx, particle_yy = np.meshgrid(particle_x, particle_y)
xp = np.column_stack((particle_xx.ravel(), particle_yy.ravel()))

# particle states
vp = np.zeros(xp.shape)
s = np.zeros((len(xp), 3))  # [sigma_xx, sigma_yy, sigma_xy]
eps = np.zeros((len(xp), 3))  # [epsilon_xx, epsilon_yy, gamma_xy]
Fp = np.tile([1.0, 0.0, 0.0, 1.0], (len(xp), 1))

# Find elements to which particles belong
ele_ids_of_particles, p_ids_in_eles = particle_element_mapping(
    xp, delta_x, delta_y, n_ele_x, n_elements)
active_elements = np.unique(ele_ids_of_particles)
active_nodes = np.unique(elements[active_elements, :])

# Compute initial particle volume
Vp = np.zeros(len(xp))
# Volume (area) of each background cell
for p_ids_in_ele in p_ids_in_eles:
    n_mp_in_element = len(p_ids_in_ele)
    if n_mp_in_element > 0:
        volume_per_mp = (delta_x * delta_y) / n_mp_in_element
        Vp[p_ids_in_ele] = volume_per_mp
Vp0 = Vp.copy()  # save initial volume
Mp = Vp * rho  # mass of particles

# Initialize stress
height = lengths[1] + spacing / 2
# Tension-positive convention: geostatic compression is negative.
s[:, 1] = -g * rho * (height - xp[:, 1])
s[:, 0] = K0 * s[:, 1]

# Plot mesh, particles
meshes = [(nodes, elements)]
def plot_mesh_and_particles(nodes, elements, particles):
    fig, ax = plt.subplots()

    # Plot the mesh
    for element in elements:
        # Get the coordinates of the element's nodes
        quad_coords = nodes[element]
        # Repeat the first point to close the quad
        quad_coords = np.vstack([quad_coords, quad_coords[0]])
        ax.plot(quad_coords[:, 0], quad_coords[:, 1], 'b-', linewidth=1.5)  # Blue lines for the mesh

    # Plot the particles
    ax.scatter(particles[:, 0], particles[:, 1], color='red', marker='o', label='Particles')

    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Mesh and Particles Visualization')
    ax.legend()

    # Equal aspect ratio
    ax.set_aspect('equal', 'box')

    plt.show()
plot_mesh_and_particles(nodes, elements, xp)

# Track the particle whose initial position matches the C++ history.
cpp_history_path = Path(__file__).resolve().parent / "mpm2Da_history.txt"
cpp_initial = np.loadtxt(cpp_history_path, comments="#", max_rows=1)
tracked_initial = cpp_initial[1:3]
tracked_particle = np.argmin(np.linalg.norm(xp - tracked_initial, axis=1))
tracked_history = [[0.0, *xp[tracked_particle], *vp[tracked_particle]]]

# Results to save
ta = []  # time
ka = []  # kinetic energy
pa = []  # potential energy
da = []  # dissipation energy
# sa = []  # strain energy
pos = []  # positions
vel = []  # velocities
stresses = []
strains = []

nsteps = int(time / dtime)

istep = 0
potential_e_t0 = np.sum(Mp * g *  xp[:, 1])


for istep in tqdm(range(nsteps), desc="Simulating"):
    # print(f"Time: {t:.5f}/{time}")
    
    # Compute potential energy
    potential_e = np.sum(Mp * g *  xp[:, 1])
    kinetic_e = np.sum(0.5 * Mp * np.linalg.norm(vp, axis=1)**2)
    dissipation_e = potential_e_t0 - kinetic_e - potential_e

    # Refresh the nodal values
    n_masses.fill(0)
    n_momentums.fill(0)
    n_iforces.fill(0)
    n_eforces.fill(0)
    gStress.fill(0)
    
    # Iterate over the computational cells (i.e., elements)
    for ele_id in active_elements:
        # node ids of the current element
        node_ids = elements[ele_id]
        # coords of the current element
        node_coords = nodes[node_ids, :]  
        # particle ids inside the current element
        material_points = p_ids_in_eles[ele_id]

        # Iterate over the material points in the current cell
        for p in material_points:
            
            # Convert global coordinate "p=(x, y)" to local coordiate "pt=(xi, eta)".
            pt = np.array(
                [(2 * xp[p, 0] - (node_coords[0, 0] + node_coords[1, 0])) / delta_x,
                 (2 * xp[p, 1] - (node_coords[1, 1] + node_coords[2, 1])) / delta_y])

            # Evaluate shape function and its derivatives with respect to local coords (xi, eta) at (x, y).
            N, dNdxi = lagrange_basis("Q4", pt)
            # Evaluate the Jacobian at current the current local coords (xi, eta).
            J0 = node_coords.T @ dNdxi
            # Get the inverse of Jacobian
            invJ0 = np.linalg.inv(J0)
            # Get the derivative of shape function with respect to global coords, i.e., (x, y)
            dNdx = dNdxi @ invJ0

            # Current stress of material points
            stress = s[p, :]
            
            # Iterate over the nodes of the current element & update nodal values by interpolating material point values
            for i, node_id in enumerate(node_ids):
                dNIdx = dNdx[i, 0]
                dNIdy = dNdx[i, 1]
                n_masses[node_id] += N[i] * Mp[p]
                n_momentums[node_id, :] += N[i] * Mp[p] * vp[p, :]
                n_iforces[node_id, 0] -= Vp[p] * (stress[0] * dNIdx + stress[2] * dNIdy)
                n_iforces[node_id, 1] -= Vp[p] * (stress[2] * dNIdx + stress[1] * dNIdy)
                n_eforces[node_id, 1] -= N[i] * Mp[p] * g

    # Total nodal force
    nforce = n_iforces + n_eforces
    
    # Compute nodal accelerations
    n_accelerations = np.zeros((n_nodes, 2))
    # valid_mass = n_masses > 1e-12
    # n_accelerations[valid_mass, :] = nforce[valid_mass, :] / n_masses[valid_mass, None]
    n_accelerations[active_nodes, :] = nforce[active_nodes, :] / n_masses[active_nodes, None]
    
    
    # Boundary conditions
    # 0. Update nodal kinematics (nodal mass, momentum -> particles)
    # 1. Velocity & acceleration to 0
    # 2. Compute nodal force.
    #   * First, map particle body force to nodal body force
    #   * Next, map particle internal force to nodel internal force
    # 3. compute_particle_kinematics.
    #   * compute_acceleration_velocity
    #       * accel = (f_e + f_i) / mass
    #       * apply_friction_constraints
    #       * vel += accel * dt
    #       * apply_velocity_constraints
    #   * compute_updated_position
    
    
    # Function to apply frictional boundary condition on given nodes
    def apply_frictional_boundary_condition(node_ids, dir_n, sign_dir_n):
        dir_t = 1 - dir_n  # Tangential direction

        # Normal and tangential accelerations
        acc_n = n_accelerations[node_ids, dir_n]
        acc_t = n_accelerations[node_ids, dir_t]

        # Tangential velocities
        vel_t = n_momentums[node_ids, dir_t] / n_masses[node_ids]

        # Determine nodes where particles are acting towards the boundary
        positive_movement_towards_boundary = (acc_n * sign_dir_n) > 0.0

        # Apply frictional boundary condition
        for idx in np.where(positive_movement_towards_boundary)[0]:
            node_id = node_ids[idx]
            acc_n_i = acc_n[idx]
            acc_t_i = acc_t[idx]
            vel_t_i = vel_t[idx]
            # print(vel_t_i)
            
            # Determine static or kinetic friction
            if vel_t_i != 0.0:  # Kinetic friction
                # compute tangential velocity at next timestep
                vel_net = dtime * acc_t_i + vel_t_i
                vel_frictional = dtime * mu * abs(acc_n_i)
                if abs(vel_net) <= vel_frictional:
                    # friction stops the particle
                    acc_t_i = -vel_t_i / dtime
                else:
                    # friction reduces the tangential acceleration
                    acc_t_i -= np.sign(vel_net) * mu * abs(acc_n_i)
                # acc_t_i -= np.sign(vel_net) * mu * abs(acc_n_i)
            else:  # Static friction
                if abs(acc_t_i) <= mu * abs(acc_n_i):
                    acc_t_i = 0.0
                else:
                    acc_t_i -= np.sign(acc_t_i) * mu * abs(acc_n_i)

            # Update tangential acceleration
            n_accelerations[node_id, dir_t] = acc_t_i

            # Update nodal force at this node
            nforce[node_id, :] = n_accelerations[node_id, :] * n_masses[node_id]
    
    apply_frictional_boundary_condition(
        np.intersect1d(bottomNodes, active_nodes), dir_n=1, sign_dir_n=-1.0)
    
    # Clamp the bottom in both directions, matching mpm2Da.cpp.
    n_momentums[bottomNodes, 0] = 0
    n_momentums[bottomNodes, 1] = 0
    nforce[bottomNodes, 0] = 0
    nforce[bottomNodes, 1] = 0
    
    # Update nomal momentum
    n_momentums += nforce * dtime    
    
    k = 0
    u = 0
    
    # Iterate over the computational cells (i.e., elements)
    for ele_id in active_elements:
        # node ids of the current element
        node_ids = elements[ele_id]
        # coords of the current element
        node_coords = nodes[node_ids, :]  
        # particle ids inside the current element
        material_points = p_ids_in_eles[ele_id]

        # Iterate over the material points in the current cell
        for p in material_points:
            # Convert global coordinate "p=(x, y)" to local coordiate "pt=(xi, eta)".
            pt = np.array(
                [(2 * xp[p, 0] - (node_coords[0, 0] + node_coords[1, 0])) / delta_x,
                 (2 * xp[p, 1] - (node_coords[1, 1] + node_coords[2, 1])) / delta_y])

            N, dNdxi = lagrange_basis("Q4", pt)
            J0 = node_coords.T @ dNdxi
            invJ0 = np.linalg.inv(J0)
            dNdx = dNdxi @ invJ0
            Lp = np.zeros((2, 2))
            acceleration_p = np.zeros(2)
            velocity_pic = np.zeros(2)
            
            for i, node_id in enumerate(node_ids):  
                vI = np.zeros(2)
                acceleration_p += N[i] * nforce[node_id, :] / n_masses[node_id]
                xp[p, :] += dtime * N[i] * n_momentums[node_id, :] / n_masses[node_id]
                vI = n_momentums[node_id, :] / n_masses[node_id]  # nodal velocity
                velocity_pic += N[i] * vI
                Lp += vI.reshape(2, 1) @ dNdx[i, :].reshape(1, 2)  # particle velocity gradient

            velocity_flip = vp[p, :] + dtime * acceleration_p
            vp[p, :] = pic_ratio * velocity_pic + (1.0 - pic_ratio) * velocity_flip

            F = Fp[p, :].reshape(2, 2) @ (np.eye(2) + Lp * dtime)
            Fp[p, :] = F.flatten()
            Vp[p] = np.linalg.det(F) * Vp0[p]
            dEps = dtime * 0.5 * (Lp + Lp.T)
            
            D = elasticity_matrix(E, nu)
            alpha, k = drucker_prager_parameters(phi, c)
            new_sig, new_epsE, _ = update_stress_strain(
                stress_n=s[p, :],
                strain_n=eps[p, :],
                strain_increment=np.array(
                    [dEps[0, 0], dEps[1, 1], 2.0 * dEps[0, 1]]
                ),
                D=D,
                alpha=alpha,
                k=k)

            s[p, :] = new_sig
            eps[p, :] = new_epsE

            for i, node_id in enumerate(node_ids):
                gStress[:, node_id] += N[i] * s[p, :]

    stresses.append(s.copy())
    strains.append(eps.copy())
    pos.append(xp.copy())
    vel.append(vp.copy())
    tracked_history.append([t + dtime, *xp[tracked_particle], *vp[tracked_particle]])
    
    # Update particle-element mapping (elements to which particles belong)
    ele_ids_of_particles, p_ids_in_eles = particle_element_mapping(
        xp, delta_x, delta_y, n_ele_x, n_elements)
    active_elements = np.unique(ele_ids_of_particles)
    active_nodes = np.unique(elements[active_elements, :])

    if istep % interval == 0:
        # result = {"positions": xp, "stress": s}
        # save_dir = f"./results_vtk_column/"
        # visualizer.save_vtk(istep, result, save_dir)
        
        ta.append(t)
        ka.append(kinetic_e)
        pa.append(potential_e)
        da.append(dissipation_e)

        
    t += dtime

# Compare the tracked particle. There is no closed-form analytical solution for
# granular-column collapse, so C++ and Python are compared directly.
tracked_history = np.asarray(tracked_history)
cpp_history = np.genfromtxt(
    cpp_history_path, comments="#", invalid_raise=False, filling_values=np.nan
)
cpp_history = cpp_history[np.all(np.isfinite(cpp_history), axis=1)]
cpp_history = cpp_history[cpp_history[:, 0] <= time]
python_at_cpp_time = np.column_stack(
    [np.interp(cpp_history[:, 0], tracked_history[:, 0], tracked_history[:, i])
     for i in range(1, 5)]
)

fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
labels = (("Velocity x (m/s)", 3, 2), ("Velocity y (m/s)", 4, 3),
          ("Displacement x (m)", 1, 0), ("Displacement y (m)", 2, 1))
for axis, (ylabel, cpp_col, python_col) in zip(axes.ravel(), labels):
    cpp_values = cpp_history[:, cpp_col]
    python_values = python_at_cpp_time[:, python_col]
    if "Displacement" in ylabel:
        cpp_values = cpp_values - cpp_values[0]
        python_values = python_values - python_values[0]
    axis.plot(cpp_history[:, 0], cpp_values, "C0-", label="C++ cubic B-spline")
    axis.plot(cpp_history[:, 0], python_values, "C1--", label="Python Q4")
    axis.set_ylabel(ylabel)
    axis.grid(True, alpha=0.3)
    axis.legend()
axes[1, 0].set_xlabel("Time (s)")
axes[1, 1].set_xlabel("Time (s)")
fig.tight_layout()
comparison_path = Path(__file__).resolve().parent / "mpm2Da_comparison.png"
fig.savefig(comparison_path, dpi=180)
print(f"Saved comparison to {comparison_path}")
raise SystemExit

# Plot kinetic and strain energy
plt.figure(figsize=(10, 6))
plt.plot(ta, ka, label='Kinetic Energy (ka)', color='b', linewidth=2)
plt.plot(ta, pa, label='Potential Energy (sa)', color='r', linewidth=2)
plt.plot(ta, da, label='Dissipation Energy (sa)', color='k', linewidth=2)
plt.xlabel('Time (ta)', fontsize=14)
plt.ylabel('Energy', fontsize=14)
plt.title('Kinetic and Strain Energy Over Time', fontsize=16)
plt.legend()
plt.grid(True)
plt.show()

from matplotlib.animation import FuncAnimation
from IPython.display import HTML

# Create figure and axis
fig, ax = plt.subplots(figsize=(5, 3.5))
ax.set_xlim(0, l)
ax.set_ylim(0, l)
ax.set_xlabel('X Position')
ax.set_ylabel('Y Position')
ax.set_title('Granular Column Collapse')

# Initialize scatter plot
scatter = ax.scatter([], [])

# Animation update function
def update(frame):
    scatter.set_offsets(pos[frame * interval])
    ax.set_title(f'Time: {ta[frame]:.3f}s')
    ax.set_aspect('equal')
    return scatter,

# Create animation
# Using the same interval as your save interval
anim = FuncAnimation(
    fig, update, 
    frames=len(ta),  # number of frames matches saved timesteps
    interval=50,     # 50ms between frames
    blit=True
)

# Display animation in notebook
HTML(anim.to_jshtml())

# Create animation of particles colored by yy stress
fig, ax = plt.subplots(figsize=(10, 7))

# Set axis properties
ax.set_aspect('equal')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('Particle Distribution with YY Stress')
ax.grid(True)

# Set fixed axis limits based on domain size
ax.set_xlim(0, l)
ax.set_ylim(0, l)

# Create a normalized colormap - use consistent limits for better visualization
vmin = min([np.min(stresses[0][:, 1]) for i in range(len(stresses))])
vmax = max([np.max(stresses[0][:, 1]) for i in range(len(stresses))])

# Create a ScalarMappable for the colorbar outside the animation function
sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm.set_array([])  # You need to set an array (empty in this case)
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('YY Stress Component (σyy)')

def update_stress_plot(frame):
    # Clear previous frame
    ax.clear()
    
    # Get yy stress component (s[:, 1]) for current frame
    yy_stress = stresses[frame * interval][:, 1]
    
    # Plot particles colored by stress
    scatter = ax.scatter(
        pos[frame * interval][:, 0], 
        pos[frame * interval][:, 1],
        c=yy_stress, 
        cmap='coolwarm', 
        vmin=vmin,
        vmax=vmax,
        s=20,
        alpha=0.8
    )
    
    # Set axis properties
    ax.set_aspect('equal')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(f'Particle Distribution with YY Stress (Step {frame * interval})')
    ax.grid(True)
    ax.grid(True, linestyle='-', linewidth=0.5, alpha=0.7)
    ax.set_xticks(np.arange(0, l, l/n_ele_x)) 
    ax.set_yticks(np.arange(0, l, l/n_ele_y)) 

    # Set fixed axis limits based on domain size
    ax.set_xlim(0, l)
    ax.set_ylim(0, l)
    plt.tight_layout()
    
    return scatter,

# Create animation
from matplotlib.animation import FuncAnimation

anim = FuncAnimation(
    fig, 
    update_stress_plot, 
    frames=len(ta),
    interval=100,  # milliseconds between frames
    blit=False
)

# Display the animation in the notebook
from IPython.display import HTML
HTML(anim.to_jshtml())

# Create animation of particles colored by yy strain
fig, ax = plt.subplots(figsize=(10, 7))

# Set axis properties
ax.set_aspect('equal')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('Particle Distribution with YY Strain')
ax.grid(True)

# Set fixed axis limits based on domain size
ax.set_xlim(0, l)
ax.set_ylim(0, l)

# Create a normalized colormap - use consistent limits for better visualization
vmin = 0
vmax = 0.08

# Create a ScalarMappable for the colorbar outside the animation function
sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm.set_array([])  # You need to set an array (empty in this case)
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('YY Strain Component (εyy)')

def update_strain_plot(frame):
    # Clear previous frame
    ax.clear()
    
    # Get yy strain component (eps[:, 1]) for current frame
    xy_strain = strains[frame * interval][:, 2]
    
    # Plot particles colored by strain
    scatter = ax.scatter(
        pos[frame * interval][:, 0], 
        pos[frame * interval][:, 1],
        c=xy_strain, 
        cmap='coolwarm', 
        vmin=vmin,
        vmax=vmax,
        s=20,
        alpha=0.8
    )
    
    # Set axis properties
    ax.set_aspect('equal')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(f'Particle Distribution with YY Strain (Step {frame * interval})')
    ax.grid(True)
    ax.grid(True, linestyle='-', linewidth=0.5, alpha=0.7)
    ax.set_xticks(np.arange(0, l, l/n_ele_x)) 
    ax.set_yticks(np.arange(0, l, l/n_ele_y)) 

    # Set fixed axis limits based on domain size
    ax.set_xlim(0, l)
    ax.set_ylim(0, l)
    plt.tight_layout()
    
    return scatter,

# Create animation
from matplotlib.animation import FuncAnimation

strain_anim = FuncAnimation(
    fig, 
    update_strain_plot, 
    frames=len(ta),
    interval=100,  # milliseconds between frames
    blit=False
)

# Display the animation in the notebook
from IPython.display import HTML
HTML(strain_anim.to_jshtml())
