import vtk
from PyQt5.QtWidgets import QApplication
from model.GUI_handler import GUI
from model._sensors import *
from model.material_handler import *
from model.mesh_converter import *
from model.mesh_helper import *
from model.pressure_script import *
from model.recording_manager import *
import io
import sys
from model._activation import ActivationClass
from model.stimulis import *

# Set to True to enable the terminal output,
# otherwise the output will be redirected to the log file (maybe it is faster this way)
TERMINAL_OUTPUT = True


class Main:

    def __init__(self, mesh_boost, stimuli, sensors, rank_material):
        self.stimuli = stimuli
        self.sensors = sensors

        self.gui = GUI(mesh_boost, rank_material, stimuli, sensors)
        self.add_interactive_events()

    def add_interactive_events(self):

        # Add the event on the press of the space bar, apply the force
        self.gui.plotter.add_key_event(
            'space', lambda: apply_volume_pressure(self.gui)
        )

        # If the enter button is pressed, the interactive mode is toggled
        self.gui.plotter.add_key_event('Control_L', self.toggle_interactive)
        self.gui.plotter.show()

    def toggle_interactive(self):

        if self.gui.plotter.interactor.GetInteractorStyle().__class__ == ActivationClass:
            self.gui.plotter.interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
            self.gui.add_mode_text("Interactive")
        else:
            self.gui.add_mode_text("Activation")
            self.gui.plotter.interactor.SetInteractorStyle(ActivationClass(
                gui=self.gui,
            ))


if not TERMINAL_OUTPUT:
    # create a text trap and redirect stdout
    text_trap = io.StringIO()
    sys.stdout = text_trap

app = QApplication(sys.argv)

# _mesh_boost = GridMesh(30, 30, z_function=concave, layers=3)
# _sensors = SensorGrid(10, 10, _mesh_boost)
# _sensors = RandomSensors(20, _mesh_boost)

_mesh_boost = ArmMesh()
_sensors = SensorArm(_mesh_boost)
# _sensors = SensorPatchesFromFile("../patches/circle.csv", _mesh_boost, n_patches=4)

# _stimuli = Cylinder(radius=5.0, height=1.0)
# _stimuli = Cuboid(2.0, 4.0, 2.0)
_stimuli = Sphere(radius=1.6)


# force_handler = pressure_script.StimuliPressure(_stimuli, 10, rubber)
# FENICS(_mesh_boost, rubber, _sensors).apply_pressure(force_handler)

# ArtificialRecording(_sensors, _mesh_boost)

Main(_mesh_boost, _stimuli, _sensors, silicon)
app.exec_()


"""
TODO:
Stress relaxation process
Apply pressure to sensors
Add README
Add random mesh grid
Add maximum displacement (look at the fucking foam at 8N)

Very important!!!
Check why displacement is 0 when there is contact!!!

https://docs.google.com/document/d/1y1DOSzD8cJKVa9my3Vqec7bWuaSYQyegPPBmfde-MQU/edit

"""
