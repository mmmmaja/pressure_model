from PyQt5.QtWidgets import QApplication
from model.GUI_handler import GUI
from model.pressure_script import *
import io
import sys
from model._activation import ActivationClass
from model.setup import _stimuli, _sensors, _mesh_boost, _material

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


Main(_mesh_boost, _stimuli, _sensors, _material)
app.exec_()


"""
TODO:

Add README

Add maximum displacement (look at the fucking foam at 8N)

Very important!!!
Check why displacement is 0 when there is contact!!!

Stress needs to be present during stress relaxation

Add max pressure per material (in the material handler)

Add stress sensors vs strain sensors | DONE
+ investigate displacement

Fill the robotic arm


https://docs.google.com/document/d/1y1DOSzD8cJKVa9my3Vqec7bWuaSYQyegPPBmfde-MQU/edit


"""
