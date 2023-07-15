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

        # Create a timer id for the timer callback
        self.timer_id = None
        # Timer interval in ms between each callback
        self.dt = 50

        # Override the default mouse events (disable rotation)
        self.AddObserver("LeftButtonPressEvent", self.left_button_press_event)
        self.AddObserver("MiddleButtonPressEvent", self.middle_button_press_event)
        self.AddObserver("RightButtonPressEvent", self.right_button_press_event)
        self.AddObserver("LeftButtonReleaseEvent", self.left_button_release_event)
        # self.AddObserver("MouseMoveEvent", self.mouse_move_event)

        super().__init__(*args, **kwargs)

    def left_button_press_event(self, obj, event):

        # Start a timer callback that triggers every 50 ms (you can adjust this value)
        self.timer_id = self.GetInteractor().CreateRepeatingTimer(50)

        # Add an observer for the timer event
        self.AddObserver('TimerEvent', self.timer_callback)

    def left_button_release_event(self, obj, event):
        """
        Function that is triggered when the left mouse button is released.
        It destroys the timer that was created in the left_button_press_event function.

        :param obj: object that triggered the event
        :param event: event that was triggered
        """

        # Destroy the timer when the mouse button is released
        if self.timer_id is not None:
            self.GetInteractor().DestroyTimer(self.timer_id)
            self.timer_id = None
        self.gui.sensors.relax()

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
        # Get the mouse coordinates and pick the cell
        x, y = self.GetInteractor().GetEventPosition()
        self.pick_cell(x, y)

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
            apply_stimuli_pressure(self.gui, self.gui.stimuli, self.picker, cell_id)
