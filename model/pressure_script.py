from abc import abstractmethod
import numpy as np

from fenics import FENICS
from relaxation_script import StressRelaxation

"""

The ForceHandler class is a parent class for all the forces.

    The VolumePressure class is a force that is applied to the whole mesh.
    
    The StimuliPressure class is a force that is applied with a stimulus.
    (Force is dependant on the stimulus shape and position)
"""


class PressureHandler:

    def __init__(self):
        # Parent class for all the forces
        pass

    @abstractmethod
    def get_pressure(self, vertex_coordinates: np.ndarray) -> float:
        """
        TODO override in subclasses
        :param vertex_coordinates: A 3D coordinates of the vertex in the mesh
        :return: A float value representing the pressure acting on the vertex

        Inverse square law
        """


class VolumePressure(PressureHandler):

    def __init__(self, pressure_strength):
        """
        Stable force that is applied to the whole mesh
        :param pressure_strength: float value representing the strength of the pressure applied
        """
        self.pressure_strength = pressure_strength
        super().__init__()

    def get_pressure(self, vertex_coordinates: np.ndarray) -> float:
        # Pressure is the same for all the vertices
        return self.pressure_strength


def is_inside(vertex, vertex_coordinates):
    """
    :param vertex: A 3D vertex in the mesh
    :param vertex_coordinates: 3D coordinates of the rectangle plane
    :return: True if the vertex is inside the plane, False otherwise
    """
    if len(vertex_coordinates) != 4:
        return False

    x, y, z = vertex
    x1, y1, z1 = vertex_coordinates[0]
    x2, y2, z2 = vertex_coordinates[1]
    x3, y3, z3 = vertex_coordinates[2]
    x4, y4, z4 = vertex_coordinates[3]

    # Check if the vertex is within the rectangle's boundaries
    if min(x1, x2, x3, x4) <= x <= max(x1, x2, x3, x4) and \
            min(y1, y2, y3, y4) <= y <= max(y1, y2, y3, y4) and \
            min(z1, z2, z3, z4) <= z <= max(z1, z2, z3, z4):
        return True
    else:
        return False


class StimuliPressure(PressureHandler):

    def __init__(self, stimuli, pressure_strength):
        """
        Force that is applied to the mesh with a stimulus
        :param stimuli: Stimuli object that has its own pressure function specified
        :param pressure_strength: float value representing the strength of the pressure applied
        """
        super().__init__()
        self.stimuli = stimuli
        self.pressure_strength = pressure_strength

    def get_pressure(self, vertex_coordinates: np.ndarray) -> float:
        # Scale the force with the force strength
        pressure = self.stimuli.calculate_pressure(vertex_coordinates) * self.pressure_strength

        return pressure


"""
Handle applying the pressure to the mesh
"""

# I need to hold the reference to the timer class and destroy it
# when the simulation of the relaxation process is over
stress_relaxation_ref = None


def apply_volume_pressure(gui, relaxation=True):
    force_handler = VolumePressure(gui.PRESSURE)
    apply_pressure(gui, force_handler, relaxation)


def apply_stimuli_pressure(gui, stimuli, picker, cell_id, relaxation=True):
    # Recompute the position of the stimuli based on the cell that was picked
    if stimuli.recompute_position(picker, cell_id):
        # This object will handle assigning the pressure to the vertices
        force_handler = StimuliPressure(stimuli, gui.PRESSURE)
        # Apply the pressure to the mesh, start the FENICS computation
        apply_pressure(gui, force_handler, relaxation)


def apply_pressure(gui, force_handler, relaxation):
    """
    :param gui: GUI class
    :param force_handler: ForceHandler class instance
    :param relaxation: boolean that indicates if the stress relaxation process should be started
    Function that applies a vertex specific force or volume (stable) force across the whole mesh
    """
    global stress_relaxation_ref

    fenics = FENICS(gui.mesh_boost.sfepy_mesh, gui.mesh_material, gui.sensors)

    # Calculate the displacement
    u = fenics.apply_pressure(force_handler)
    if u is None:
        return

    # UPDATE plot and meshes
    gui.mesh_boost.update_mesh(u)
    gui.sensors.update_visualization()

    gui.draw_mesh()
    gui.draw_sensors()
    gui.plotter.update()

    if relaxation:
        # Start the stress relaxation process
        # Stop existing relaxation process
        if stress_relaxation_ref is not None:
            stress_relaxation_ref.stop()

        stress_relaxation_ref = StressRelaxation(gui)
        stress_relaxation_ref.initiate(wait=True)
