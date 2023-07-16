import time
import numpy as np


class StressRelaxation:
    """
    Standard Linear Solid (SLS) model for stress relaxation
    Viscosity -  the material's resistance to flow or change shape

    Exponential decay function for relaxation
        R(t) = e^(-t/τ)
    """

    # The force limit to stop the relaxation process
    FORCE_LIMIT = 1e-2

    def __init__(self, gui):
        self.gui = gui

    def exponential_decay(self, t):
        """
        Exponential decay function for relaxation
        R(t) = e^(-t/τ)
        :param t: current time of the simulation
        :return: the relaxation factor between 0 and 1,
            representing the proportion of stress remaining in the material
        """
        return np.exp(-t / self.gui.mesh_material.time_constant)

    def get_relaxation_displacement(self, u):
        """
        Calculate the relaxation displacement of the mesh
        :return: the current displacement dependant of the time t of the simulation
        """
        relaxation_rate = 0.1
        # Iterate through all the points in the mesh
        relaxation = np.zeros_like(u)
        for i in range(len(u)):

            # calculate the current displacement from the initial position
            current_displacement = self.gui.mesh_boost.current_vtk.points[i] - self.gui.mesh_boost.initial_vtk.points[i]
            # calculate the relaxation displacement
            relaxation_displacement = relaxation_rate * current_displacement
            # print(relaxation_displacement)
            relaxation[i] = -relaxation_displacement

        return relaxation

    def relax_iteration(self):
        """

        This function is called every dt milliseconds
        Updates the displacement of the mesh and the GUI
        """

        # Get the overall displacement of the mesh on this iteration
        u = self.gui.mesh_boost.current_vtk.points - self.gui.mesh_boost.initial_vtk.points

        # calculate the displacement from initial mesh to the current mesh
        relaxation = self.get_relaxation_displacement(u)

        return relaxation
