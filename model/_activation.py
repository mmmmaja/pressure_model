import datetime

import vtk
from model.pressure_script import *


class ActivationClass(vtk.vtkInteractorStyleTrackballCamera):
    """
    This class was added to override the default mouse events of the vtkInteractorStyleTrackballCamera class.
    It disables the rotation of the mesh when the left mouse button is pressed.
    Triggered in Activation mode
    """

    def __init__(self, parent=None, gui=None, *args, **kwargs):

        self.gui = gui

        # Mouse pressed flag for mesh activation
        self.mouse_pressed = False

        # Create a cell picker for the mesh (to map the mouse coordinates to the mesh)
        self.picker = vtk.vtkCellPicker()
        self.picker.AddPickList(self.gui.mesh_actor)

        self.mouse_pressed = False

        self.stress_relaxation = StressRelaxation(self.gui)
        # Timer interval in ms between each callback
        self.dt = 50
        self.timer_id = None

        # Override the default mouse events (disable rotation)
        self.AddObserver("LeftButtonPressEvent", self.left_button_press_event)
        self.AddObserver("MiddleButtonPressEvent", self.middle_button_press_event)
        self.AddObserver("RightButtonPressEvent", self.right_button_press_event)
        self.AddObserver("LeftButtonReleaseEvent", self.left_button_release_event)

        super().__init__(*args, **kwargs)

    def start_timer(self):
        # Create a timer id for the timer callback
        self.timer_id = self.GetInteractor().CreateRepeatingTimer(self.dt)
        self.AddObserver("TimerEvent", self.timer_callback)

    def left_button_press_event(self, obj, event):
        self.mouse_pressed = True
        if self.timer_id is None:
            self.start_timer()

    def left_button_release_event(self, obj, event):
        self.mouse_pressed = False

    def middle_button_press_event(self, obj, event):
        # Disable the middle button events
        pass

    def right_button_press_event(self, obj, event):
        # Disable the right button events
        pass

    def timer_callback(self, obj, event):
        """
        Timer callback that is triggered every dt ms.
        Triggers the pick_cell function which applies a force to the cell that was clicked.

        :param obj: object that triggered the event
        :param event: event that was triggered
        """
        relaxation = self.stress_relaxation.relax_iteration()
        new_mesh_coords = self.gui.mesh_boost.current_vtk.points.copy() + relaxation

        if self.mouse_pressed:
            # Get the mouse coordinates and pick the cell
            x, y = self.GetInteractor().GetEventPosition()
            u = self.pick_cell(x, y)
            if u is not None:
                new_mesh_coords += u

        self.update(new_mesh_coords)

    def pick_cell(self, x, y):
        """
        Function that picks the cell that was clicked and applies a force to it.
        :param x: x coordinate of the mouse
        :param y: y coordinate of the mouse
        """

        # Pick the cell that was clicked
        self.picker.Pick(x, y, 0, self.gui.plotter.renderer)
        cell_id = self.picker.GetCellId()

        # If the cell exists
        if cell_id != -1:
            print("Cell ID: ", cell_id)
            u = apply_stimuli_pressure(self.gui, self.picker, cell_id)
            if u is not None:
                self.update_activation_matrix(u)
                return u
        return None

    def update(self, new_mesh_points):

        # UPDATE plot and meshes
        self.gui.mesh_boost.update_mesh(new_mesh_points)
        self.gui.sensors.update_visualization()

        self.gui.draw_mesh()
        self.gui.draw_sensors()
        self.gui.plotter.update()

    def update_activation_matrix(self, u):
        current_time = datetime.datetime.now()
        for i in range(self.gui.mesh_boost.current_vtk.points.shape[0]):
            # Get the pressure that should be applied to the vertex
            displacement = u[i]
            displacement_norm = np.linalg.norm(displacement)

            if displacement_norm > 1e-10:
                self.gui.mesh_boost.last_activation_time[i] = current_time

                if displacement_norm > np.linalg.norm(self.gui.mesh_boost.max_displacement[i]):
                    self.gui.mesh_boost.max_displacement[i] = displacement

