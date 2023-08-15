import vtk
from PyQt5.QtWidgets import QApplication
from GUI_handler import GUI
from pressure_script import *
import io
import sys
from _activation import ActivationClass
from setup import _stimuli, _sensors, _mesh_boost, _material

"""
From this file application is started.
Setup is loaded from setup.py script.
"""

# Set to True to enable the terminal output,
# otherwise the output will be redirected to the log file (maybe it is faster this way)
TERMINAL_OUTPUT = False


class Main:

    def __init__(self, mesh_boost, stimuli, sensors, rank_material):
        """
        Initialize the main class
        :param mesh_boost: Mesh class ( activated object geometry )
        :param stimuli: Stimuli class ( object activating the mesh )
        :param sensors: Sensors class ( sensors attached to the mesh )
        :param rank_material: Material class ( material of the mesh )
        """

        self.stimuli = stimuli
        self.sensors = sensors
        self.gui = GUI(mesh_boost, rank_material, stimuli, sensors)

        self.add_interactive_events()

    def add_interactive_events(self):
        """
        Add the interactive events to the plotter
        When user presses space bar, mode of the simulation is changed
        """

        # Add the event on the press of the space bar, enable user to apply pressure
        self.gui.plotter.add_key_event(
            'space', lambda: apply_volume_pressure(self.gui)
        )

        # If the control button is pressed, the interactive mode is toggled
        self.gui.plotter.add_key_event('Control_L', self.toggle_interactive)
        self.gui.plotter.show()

    def toggle_interactive(self):
        """
        Toggle the interactive mode
        Changes the interactor style of the plotter
        """

        # If the interactor style is the activation class, change it to the default one
        if self.gui.plotter.interactor.GetInteractorStyle().__class__ == ActivationClass:
            self.gui.plotter.interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
            self.gui.add_mode_text("Interactive")
        # If the interactor style is the default one, change it to the activation class
        else:
            self.gui.add_mode_text("Activation")
            self.gui.plotter.interactor.SetInteractorStyle(ActivationClass(
                gui=self.gui
            ))


# If the terminal output is disabled, redirect the output to the log file
if not TERMINAL_OUTPUT:
    # create a text trap and redirect stdout
    text_trap = io.StringIO()
    sys.stdout = text_trap

# Start the application
app = QApplication(sys.argv)

# Create the main class and run simulation
Main(_mesh_boost, _stimuli, _sensors, _material)
app.exec_()