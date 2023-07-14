import math
import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import Delaunay, distance
from sklearn.manifold import MDS

"""
Script that handles the flattening of the surface
(Projection from 3D space onto 2D plane)
"""


def plot_mesh_2D(uv, centroid=False):
    """
    Plot the 2D mesh of the given uv coordinates
    :param uv: 2D coordinates of the mesh
    :param centroid: whether to plot the centroid of the mesh
    """

    triangles = Delaunay(uv).simplices

    plt.figure()

    # add edges
    for tri in triangles:
        i, j, k = tri
        # Use magenta thin lines for the edges
        plt.plot([uv[i, 0], uv[j, 0]], [uv[i, 1], uv[j, 1]], color='m', linewidth=0.5)
        plt.plot([uv[j, 0], uv[k, 0]], [uv[j, 1], uv[k, 1]], color='m', linewidth=0.5)
        plt.plot([uv[k, 0], uv[i, 0]], [uv[k, 1], uv[i, 1]], color='m', linewidth=0.5)

    # scatter plot of vertices
    plt.scatter(uv[:, 0], uv[:, 1])

    if centroid is True:
        centroid = np.mean(uv, axis=0)
        plt.scatter(centroid[0], centroid[1], color='g', s=18)

    plt.show()


def plot_meshes_2D(uv1, uv2):
    """
    Plot two meshes against each other in 2D plot
    (Used to create the mapping from one mesh to another)
    :param uv1: 2D coordinates of the first mesh
    :param uv2: 2D coordinates of the second mesh
    """
    triangles1, triangles2 = Delaunay(uv1).simplices, Delaunay(uv2).simplices

    plt.figure()
    # add edges
    for tri in triangles1:
        i, j, k = tri
        # Use magenta thin lines for the edges
        plt.plot([uv1[i, 0], uv1[j, 0]], [uv1[i, 1], uv1[j, 1]], color='m', linewidth=0.5)
        plt.plot([uv1[j, 0], uv1[k, 0]], [uv1[j, 1], uv1[k, 1]], color='m', linewidth=0.5)
        plt.plot([uv1[k, 0], uv1[i, 0]], [uv1[k, 1], uv1[i, 1]], color='m', linewidth=0.5)

        # add edges
    for tri in triangles2:
        i, j, k = tri
        # Use magenta thin lines for the edges
        plt.plot([uv2[i, 0], uv2[j, 0]], [uv2[i, 1], uv2[j, 1]], color='g', linewidth=0.5)
        plt.plot([uv2[j, 0], uv2[k, 0]], [uv2[j, 1], uv2[k, 1]], color='g', linewidth=0.5)
        plt.plot([uv2[k, 0], uv2[i, 0]], [uv2[k, 1], uv2[i, 1]], color='g', linewidth=0.5)

    # scatter plot of vertices
    plt.scatter(uv1[:, 0], uv1[:, 1], color='m')
    plt.scatter(uv2[:, 0], uv2[:, 1], color='g')

    plt.show()


def orthographic_projection(sensor_positions):
    """
    Project the 3D sensor positions onto a 2D plane
    Simply remove the z coordinate flattening the surface
    :param sensor_positions: 3D coordinates of the sensor positions
    """
    return sensor_positions[:, :2]


def surface_parametrization(sensor_positions):
    """
    Find a function that maps the 3D surface onto a 2D plane (the tactile map)
    while minimizing distortions.
    :param sensor_positions: 3D coordinates of the sensor positions to be mapped onto a 2D plane
    """

    # Compute the pairwise distances between all points in the mesh
    distances = sensor_positions[:, np.newaxis, :] - sensor_positions[np.newaxis, :, :]
    distances = np.sqrt((distances ** 2).sum(axis=-1))

    # Use multidimensional scaling to reduce the dimensionality to 2
    mds = MDS(n_components=2, dissimilarity='precomputed')
    coords_2d = mds.fit_transform(distances)
    return coords_2d


def compute_centroid(points):
    """
    Compute the centroid of a set of points
    :param points: A 2D array of points
    :return: The centroid point
    """
    return np.mean(points, axis=0)


def projection_onto_grid_positions(sensor_positions):
    """
    Project the sensor positions onto the grid positions from a real robotic arm
    :param sensor_positions: 3D coordinates of the sensor positions to be mapped onto hardcoded real arm positions
    :return: 3D coordinates of the sensor positions mapped onto the grid positions
    """

    nds_frontX = [0.000, 25.981, 51.962, 77.943, 12.990, 38.971, 64.952, 90.933, 0.000, 25.981, 51.962, 77.943, 12.990,
                  38.971, 64.952, 90.933, 0.000, 25.981, 51.962, 77.943, 12.990, 38.971, 64.952, 90.933, 0.000, 25.981,
                  51.962, 77.943, 12.990, 38.971, 64.952, 90.933, 0.000, 25.981, 51.962, 77.943, 12.990, 38.971, 64.952,
                  90.933, 0.000, 25.981, 51.962, 77.943, 12.990, 38.971, 64.952, 90.933]
    nds_frontY = [0.000, 0.000, 0.000, 0.000, 7.500, 7.500, 7.500, 7.500, 15.000, 15.000, 15.000, 15.000, 22.500,
                  22.500,
                  22.500, 22.500, 30.000, 30.000, 30.000, 30.000, 37.500, 37.500, 37.500, 37.500, 45.000, 45.000,
                  45.000,
                  45.000, 52.500, 52.500, 52.500, 52.500, 60.000, 60.000, 60.000, 60.000, 67.500, 67.500, 67.500,
                  67.500,
                  75.000, 75.000, 75.000, 75.000, 82.500, 82.500, 82.500, 82.500]

    sensor_positions_flattened = orthographic_projection(sensor_positions)

    # Find the closest grid position for each sensor position
    grid_positions = np.array([nds_frontX, nds_frontY]).T

    # Display the initial meshes
    plot_meshes_2D(sensor_positions_flattened, grid_positions)

    # Map the sensor positions to the grid positions

    # SCALING
    # Compute the scale factors
    scale_x = (np.max(grid_positions[:, 0]) - np.min(grid_positions[:, 0])) / (
                np.max(sensor_positions_flattened[:, 0]) - np.min(sensor_positions_flattened[:, 0]))
    scale_y = (np.max(grid_positions[:, 1]) - np.min(grid_positions[:, 1])) / (
                np.max(sensor_positions_flattened[:, 1]) - np.min(sensor_positions_flattened[:, 1]))

    # Scale the grid
    sensor_positions_flattened[:, 0] = sensor_positions_flattened[:, 0] * scale_x
    sensor_positions_flattened[:, 1] = sensor_positions_flattened[:, 1] * scale_y

    plot_meshes_2D(sensor_positions_flattened, grid_positions)

    # TRANSLATION
    # Compute the centroids
    centroid_sensor = compute_centroid(sensor_positions_flattened)
    centroid_grid = compute_centroid(grid_positions)

    # Compute the translation factors
    translation_x = centroid_grid[0] - centroid_sensor[0]
    translation_y = centroid_grid[1] - centroid_sensor[1]

    # Translate the sensor positions
    sensor_positions_flattened[:, 0] = sensor_positions_flattened[:, 0] + translation_x
    sensor_positions_flattened[:, 1] = sensor_positions_flattened[:, 1] + translation_y

    plot_meshes_2D(sensor_positions_flattened, grid_positions)

    plot_mesh_2D(sensor_positions_flattened, centroid=True)

    # Find the closest grid position for each sensor position
    # Start from the middle of the grid because that is where there is less distortion

    # Calculate the centroid of the grid
    centroid_grid = compute_centroid(grid_positions)

    # Sort grid positions based on their distance to the centroid
    distances_to_centroid = np.linalg.norm(grid_positions - centroid_grid, axis=1)
    sorted_grid_indices = np.argsort(distances_to_centroid)

    mapped_positions = np.zeros_like(sensor_positions_flattened)

    for i in sorted_grid_indices:
        # Compute the distance matrix between the current sensor position and all the grid positions
        dist_matrix = np.linalg.norm(sensor_positions_flattened[i] - grid_positions, axis=1)

        # Get the index of the closest grid point
        closest_grid_idx = np.argmin(dist_matrix)

        mapped_positions[i] = grid_positions[closest_grid_idx]

        # Remove the assigned grid position from the grid positions
        grid_positions = np.delete(grid_positions, closest_grid_idx, axis=0)

    # Convert the mapped positions to a numpy array
    mapped_positions = np.array(mapped_positions)

    plot_mesh_2D(mapped_positions)

    return mapped_positions
