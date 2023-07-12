import math

import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import Delaunay


def encloses(triangle, point):
    """
    Check if a point lays within a 2D triangle
    :param triangle: Consists of 3 points in 2D space
    :param point: Point in 2D space
    :return: True if point lays within triangle, False otherwise
    """

    # Use baycentric coordinates to check if point lays within triangle

    # Any point P within the triangle can be expressed as:
    # P = A + u * (C - A) + v * (B - A)

    # Where A, B, C are the vertices of the triangle
    # u, v are the barycentric coordinates of P
    # Solve for u and v

    # Get the vectors of the triangle
    v0 = triangle[2] - triangle[0]
    v1 = triangle[1] - triangle[0]
    v2 = point - triangle[0]

    u = (np.dot(v1, v1) * np.dot(v2, v0) - np.dot(v1, v0) * np.dot(v2, v1)) \
        / (np.dot(v0, v0) * np.dot(v1, v1) - np.dot(v0, v1) * np.dot(v1, v0))

    v = (np.dot(v0, v0) * np.dot(v2, v1) - np.dot(v0, v1) * np.dot(v2, v0)) \
        / (np.dot(v0, v0) * np.dot(v1, v1) - np.dot(v0, v1) * np.dot(v1, v0))

    # Check if point lays within triangle
    if v < 0 or u < 0 or u > 1 or v > 1:
        return False
    if u + v > 1:
        return False
    return True


def get_enclosing_triangle(triangles, uv, x_rc):
    """
    Find the triangle in which a given point lays
    :param triangles: list of delaunay simplices
    :param uv: 2D mapping of the mesh
    :param x_rc: point in 2D space that we need to find the enclosing triangle for
    :return: position of the triangle, pressure of the triangle
    """

    # Find triangle in which given point lays
    for triangle in triangles:

        # Get real position of the triangle with respect to 2D mapping and resolution
        position = []
        for index in triangle:
            position.append(uv[index])

        # Find out if x_rc lays within this triangle
        if encloses(position, x_rc):
            return triangle


def get_pixel_color(sensor_readings, p_rc):
    # Get the maximum pressure registered in the sensors
    P = max(sensor_readings)
    if P == 0:
        return 0
    else:
        return math.floor((255 / P) * p_rc)


def get_area(p1, p2, p3):
    # Compute the area of a triangle in 2D
    return abs(0.5 * np.linalg.norm(np.cross(p2 - p1, p3 - p1)))


def get_pressure(triangles, uv, x_rc, sensor_readings):
    """
    Compute the pressure in p[x(r,c)] using the barycentric interpolation
    :param triangles: list of delaunay simplices
    :param uv: 2D mapping of the mesh
    :param x_rc: point in 2D space that we need to find the pressure for
    :param sensor_readings: list of sensor outputs
    :return: pressure in p[x(r,c)]
    """

    # Find the enclosing triangle of x(r,c)
    # Consists of 3 mapped points onto 2D surface
    triangle = get_enclosing_triangle(triangles, uv, x_rc)
    if triangle is None:
        return 0.0

    # Get the positions of the triangle points
    m_j, m_k, m_h = uv[triangle]

    # Get the pressure registered in the corresponding sensors
    p_j, p_k, p_h = sensor_readings[triangle]

    # Get the areas of the triangles
    A_total = get_area(m_j, m_k, m_h)
    A1 = get_area(x_rc, m_k, m_h)
    A2 = get_area(m_j, x_rc, m_h)
    A3 = get_area(m_j, m_k, x_rc)

    # Compute the pressure in p[x(r,c)] using the barycentric interpolation
    return (A1 * p_j + A2 * p_k + A3 * p_h) / A_total


def show_image(image, uv=None, triangles=None):

    plt.imshow(image, cmap='plasma')

    if uv is not None and triangles is not None:
        # Normalize the uv coordinates to the size of the image
        x_min, x_max = np.min(uv[:, 0]), np.max(uv[:, 0])
        y_min, y_max = np.min(uv[:, 1]), np.max(uv[:, 1])

        uv[:, 0] = (uv[:, 0] - x_min) / (x_max - x_min) * image.shape[1]
        uv[:, 1] = (uv[:, 1] - y_min) / (y_max - y_min) * image.shape[0]

        alpha, line_width = 0.1, 0.5

        # add edges
        for tri in triangles:
            i, j, k = tri[0], tri[1], tri[2]
            plt.plot(
                [uv[i, 0], uv[j, 0]], [uv[i, 1], uv[j, 1]], color='black', linewidth=line_width, alpha=alpha
            )
            plt.plot(
                [uv[j, 0], uv[k, 0]], [uv[j, 1], uv[k, 1]], color='black', linewidth=line_width, alpha=alpha
            )
            plt.plot(
                [uv[k, 0], uv[i, 0]], [uv[k, 1], uv[i, 1]], color='black', linewidth=line_width, alpha=alpha
            )

        # scatter plot of vertices
        plt.scatter(uv[:, 0], uv[:, 1], s=2, color='black', alpha=alpha)

    plt.show()


def plot_mesh_2D(uv, triangles):
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

    plt.show()


def orthographic_projection(sensor_positions):
    return sensor_positions[:, :2]


def form_contact_image(sensor_positions, sensor_readings, resolution):
    """
    Form the contact image from the sensor positions and the sensor readings

    Made based on the paper:
    Human Hand Recognition From Robotic Skin Measurements in Human-Robot Physical Interactions

    :param sensor_positions: list of sensor positions (x, y, z)
    :param sensor_readings: list of sensor stress readings
    :param resolution: resolution of the contact image
    :return: the contact image for this reading
    """

    # Get the projection of the sensor positions onto 2D space
    uv = orthographic_projection(sensor_positions)

    # Perform a 3D Delaunay triangulation on the sensor positions
    triangles = Delaunay(uv).simplices

    # plot_mesh_2D(uv, triangles)

    # Create a grid for the image
    x_min, y_min = min(uv[:, 0]), min(uv[:, 1])
    x_max, y_max = max(uv[:, 0]), max(uv[:, 1])

    resolution_x = int((x_max - x_min) * resolution)
    resolution_y = int((y_max - y_min) * resolution)
    print("image dimensions: ", resolution_x, "x", resolution_y, "pixels")

    # Create a meshgrid for the image
    X = np.linspace(start=x_min, stop=x_max, num=resolution_x)
    Y = np.linspace(start=y_min, stop=y_max, num=resolution_y)

    # Create an empty image
    image = np.zeros((resolution_x, resolution_y))

    # Iterate through the grid and compute the pressure at each point x_rc
    for r in range(len(X)):
        for c in range(len(Y)):
            # Get the coordinates of the point x_rc
            x_rc = np.array([X[r], Y[c]])
            # Find the pressure at x_rc by means of barycentric interpolation
            p_rc = get_pressure(triangles, uv, x_rc, sensor_readings)
            image[r, c] = get_pixel_color(sensor_readings, p_rc)

    # Map the grid onto the image
    show_image(image)

    return image

