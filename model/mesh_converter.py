import pyvista
from model.fenics import SfepyMesh
from abc import abstractmethod
import numpy as np
import meshio
import pyvista as pv

from model.mesh_helper import *


def convert_to_vtk(path):
    """
    Convert a mesh file to a vtk object (UnstructuredGrid)
    :param path: Path to the mesh file
    """
    return pv.UnstructuredGrid(path)


class MeshBoost:

    # This is a parent class for all the meshes in the project (the VTK version of the mesh)
    # All meshes consist of hexahedrons (cells) and vertices

    def __init__(self):
        # Path to the .mesh file
        self.path = '../meshes/mesh.mesh'

        # The initial mesh object in VTK format
        self.initial_vtk, self.top_region_ids, self.bottom_region_ids = self.create_mesh()

        # Represents the current displacements of the mesh
        self.current_vtk = self.initial_vtk.copy()

        # The mesh object in Sfepy format for solvers
        self.sfepy_mesh = self.create_sfepy_mesh()

    @abstractmethod
    def create_mesh(self):
        """
        TODO override in subclasses
        :return: The mesh object in VTK format
        """

    def override_mesh(self, u):
        """
        Override the vtk version of the mesh with a new one.
        :param u: The displacements of the mesh to be added to the initial mesh

        If the displacement distorts the mesh too much, the override is not performed.
        """

        # Create the copy of the old mesh
        mesh_copy = self.current_vtk.copy()
        # Override the vtk version of the mesh
        self.current_vtk.points = self.initial_vtk.points.copy() + u

        # Create the sfepy version of the mesh for the solver
        sfepy_mesh = self.create_sfepy_mesh()

        if not sfepy_mesh:
            # If the mesh was invalid, restore the old mesh
            self.current_vtk = mesh_copy
        else:
            # Otherwise, override the sfepy version of the mesh
            self.sfepy_mesh = sfepy_mesh

    def update_mesh(self, u):
        """
        Update the vtk version of the mesh with a new one.
        :param u: The displacements of the mesh to be added to the current mesh

        If the displacement distorts the mesh too much, the update is not performed.
        """

        # Create the copy of the old mesh
        mesh_copy = self.current_vtk.copy()
        # Update the vtk version of the mesh
        self.current_vtk.points = self.current_vtk.points.copy() + u

        # Create the sfepy version of the mesh for the solver
        sfepy_mesh = self.create_sfepy_mesh()

        if not sfepy_mesh:
            # If the mesh was invalid, restore the old mesh
            self.current_vtk = mesh_copy
        else:
            # Otherwise, override the sfepy version of the mesh
            self.sfepy_mesh = sfepy_mesh

    def create_sfepy_mesh(self):
        """
        Create a mesh in Sfepy format. The mesh is created from the current vtk mesh.
        """

        sfepy_mesh = SfepyMesh(self)
        # sfepy_mesh.create()
        # return sfepy_mesh

        # Check the validity of the mesh
        if not sfepy_mesh.create():
            return None

        if sfepy_mesh.validate():
            return sfepy_mesh
        else:
            print('Mesh is invalid')
            return None

    def get_vertex_ids_from_coords(self, cell_coords):
        """
        Given a mesh and cell coordinates, find the matching vertex IDs in the mesh.
        """
        mesh_points = self.current_vtk.points
        vertex_ids = []

        for cell_point in cell_coords:
            for i, mesh_point in enumerate(mesh_points):
                if np.allclose(cell_point, mesh_point):
                    vertex_ids.append(i)
                    break

        return vertex_ids


class GridMesh(MeshBoost):

    # Thickness of the mesh
    THICKNESS = 2.73

    def __init__(self, width, height, z_function=flat, layers=2):
        """
        Defines the mesh as a grid of vertices

        :param width: dimension of the grid
        :param height: dimension of the grid
        :param z_function: function that defines the height of the grid
        :param layers: number of layers in the grid (in z direction)
        """
        self.width, self.height = width, height
        self.z_function, self.layers = z_function, layers

        super().__init__()

    def create_mesh(self):

        # Create the vertices of the tetrahedrons
        # Define two regions of the mesh
        all_layers = []
        for i in range(self.layers):
            all_layers.append([])

        # Add the vertices to the regions
        for i in range(self.height):
            for j in range(self.width):
                all_layers[0].append([
                    i, j, self.z_function(i, j, self.width, self.height)
                ])
                for k in range(1, self.layers):
                    all_layers[k].append([
                        i, j, (self.layers - k - 1) * self.THICKNESS
                    ])

        # ADJUST THE TOP LAYER (where the deformations will be applied)

        # Convert all_layers[0] to numpy array for vectorized operation
        top_layer = np.array(all_layers[0])

        # Adjust the z-coordinates and ensure all are non-negative
        top_layer[:, 2] = np.maximum(top_layer[:, 2] - np.amin(top_layer[:, 2]), 0)

        # Add the thickness to the top vertices
        top_layer[:, 2] += (self.layers - 1) * self.THICKNESS

        # Convert the adjusted numpy array back to list and assign it back to all_layers[0]
        all_layers[0] = top_layer.tolist()

        # Combine all the vertices
        vertices = [vertex for vertices in all_layers for vertex in vertices]
        n = len(vertices) // self.layers

        # Create the cells (give the indices of vertices of each hexahedron)
        cells = []
        for i in range(self.height - 1):
            for j in range(self.width - 1):
                # Go through all the layers
                for k in range(self.layers - 1):
                    a_top, b_top = i * self.width + j + k * n, i * self.width + j + 1 + k * n
                    c_top, d_top = (i + 1) * self.width + j + 1 + k * n, (i + 1) * self.width + j + k * n
                    cells.append([
                        a_top, b_top, c_top, d_top,
                        a_top + n, b_top + n,  c_top + n, d_top + n
                    ])

        # Create a top region (Where the displacements happen)
        top_region_ids = range(n)
        # Create a bottom region (Where the boundary conditions apply so that the positions are fixed)
        bottom_region_ids = range((self.layers - 1) * n, self.layers * n)

        # Create hexahedron mesh in the meshio format
        mesh = meshio.Mesh(points=vertices, cells={"hexahedron": cells})
        meshio.write(self.path, mesh)
        # Read the mesh as a vtk mesh
        return convert_to_vtk(self.path), top_region_ids, bottom_region_ids


class ArmMesh(MeshBoost):

    # Thickness of the mesh
    # Make sure that thickness is not too big ( otherwise the mesh will crash :( )
    THICKNESS = 1.2

    OBJ_PATH = '../meshes/model_kfadrat.obj'

    def __init__(self):
        super().__init__()

    def create_mesh(self):
        """
        Extruding the mesh along a direction perpendicular to the faces
        :return:
        """

        # Load the input mesh
        mesh = pv.read(self.OBJ_PATH)

        # Compute normals for each face
        normals = mesh.face_normals

        # Create a dictionary to handle shared vertices
        # base vertex coordinates  : [extruded vertex coordinates, extruded vertex id]
        vertices_dict = {}

        # Save the top and bottom region ids
        top_region_ids, bottom_region_ids = [], []

        # For each face in the input mesh find the extruded face coordinates
        for i in range(mesh.n_cells):

            # Get the point coordinates for this face
            face_coords = mesh.get_cell(i).points

            # Loop over the vertices of the base (non-extruded) face
            for k in range(4):

                # Create a unique identifier for each vertex:
                # - a tuple of the vertex's coordinates
                identifier = tuple(np.round(face_coords[k], 4))

                if identifier not in vertices_dict:
                    extruded_coords = tuple(np.round(face_coords[k] - normals[i] * self.THICKNESS, 4))
                    vertices_dict[identifier] = extruded_coords

        # Make another pass to form the connections and create the cells

        # Create the vertices and cells for the output mesh (vtk format)
        vertices, cells = [], []
        # For each face in the input mesh, create a hexahedral cell in the output mesh
        for i in range(mesh.n_cells):
            face_coords = mesh.get_cell(i).points

            # The cell is defined by eight points: four points on the base face
            # and their corresponding extruded points
            cell_base, cell_extruded = [], []
            for k in range(4):
                identifier = tuple(np.round(face_coords[k], 4))

                # Handle the base vertices
                # Check if the vertex is already in the list of vertices
                if identifier not in vertices:
                    vertices.append(identifier)
                    vertex_id = vertices.index(identifier)
                else:
                    vertex_id = vertices.index(identifier)
                cell_base.append(vertex_id)  # base point
                top_region_ids.append(vertex_id)

                # Handle the extruded vertices
                extruded_point = vertices_dict[identifier]
                if extruded_point not in vertices:
                    vertices.append(extruded_point)
                    vertex_id = vertices.index(extruded_point)
                else:
                    vertex_id = vertices.index(extruded_point)
                cell_extruded.append(vertex_id)  # extruded point
                bottom_region_ids.append(vertex_id)

            # Combine the front and back faces to form the cell
            cell = cell_base[::-1] + cell_extruded[::-1]
            cells.append(cell)

        mesh = meshio.Mesh(points=vertices, cells={"hexahedron": cells})
        meshio.write(self.path, mesh)
        # Read the mesh as a vtk mesh
        return convert_to_vtk(self.path), top_region_ids, bottom_region_ids


def display_obj_file(path):
    """
    Display the obj file from a given path on the pyvista plot
    """
    reader = pyvista.get_reader(path)
    mesh = reader.read()
    mesh.plot(cpos='yz', show_scalar_bar=False)
