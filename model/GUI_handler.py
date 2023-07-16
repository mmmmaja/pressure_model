import numpy as np
import pyvistaqt as pvqt  # For updating plots real time
from PyQt5.QtWidgets import QAction  # For the custom button

from model.pressure_script import apply_volume_pressure
from model.recording_manager import Recording


"""
Script that handles the GUI of the application
Pyvistaqt plotting is used
"""


class GUI:

    def __init__(self, mesh_boost, mesh_material, stimuli, sensors):

        # Mesh object that has all the information about the mesh
        self.mesh_boost = mesh_boost
        # Rank material for rendering the mesh
        self.mesh_material = mesh_material
        # Stimuli object just for reference for the user
        self.stimuli = stimuli
        # List of sensors that record pressure
        self.sensors = sensors

        # Define the sensor output type
        self.sensor_output_type = 'stress'

        # Define the pressure
        self.PRESSURE = 0.001
        # Change in the pressure when on the event
        self.pressure_dt = 0.01

        # Define all the actors present in the scene
        self.mesh_actor = None
        self.stimuli_actor = None
        self.sensor_actor = None
        self.material_text_actor = None
        self.mode_text_actor = None
        self.force_indicator_actor = None
        self.sensor_output_text_actor = None

        # Define the actions
        self.record_action = None
        self.stop_record_action = None
        # Recording object from recording_manager.py
        self.recording = None

        # Define the plotter (pyvistaqt)
        self.plotter = pvqt.BackgroundPlotter()
        self.text_color = 'eaf5ff'
        self.background_color = '282a36'
        self.plotter.set_background(self.background_color)

        # Add all the actors to the scene
        self.draw_mesh()
        self.draw_stimuli()
        self.draw_sensors()
        self.add_material_text()
        self.add_mode_text('Interactive')
        self.add_recording_actions()
        # Add sensor output switch action
        self.add_sensor_output_switch_action()
        self.add_sensor_output_text(self.sensor_output_type)

        # Add the interactive events to the plotter
        self.add_force_events()
        self.add_pressure_indicator()

        # Update the plotter and show it
        self.plotter.update()
        self.plotter.show()

    def add_material_text(self):
        # Add the material name to the plotter
        text = self.mesh_material.name
        text += '\nE: ' + str(self.mesh_material.young_modulus)
        text += '\nv: ' + str(self.mesh_material.poisson_ratio)

        self.material_text_actor = self.plotter.add_text(
            text, position='lower_left', font_size=8, color=self.text_color
        )

    def draw_mesh(self):
        """
        Adds the mesh in the .vtk format to the plotter
        """
        visual_properties = self.mesh_material.visual_properties
        self.mesh_actor = self.plotter.add_mesh(
            self.mesh_boost.current_vtk,
            show_edges=True,
            smooth_shading=True,
            show_scalar_bar=False,
            edge_color='white',
            color=visual_properties['color'],
            specular=visual_properties['specular'],
            metallic=visual_properties['metallic'],
            roughness=visual_properties['roughness'],
            name='initial_mesh',
            # opacity=0.99,
        )

    def draw_stimuli(self):
        return
        # Draw the stimuli object in the plotter
        stimuli = self.stimuli.create_visualization()

        # Get the minimal x and y coordinates of the mesh
        min_x = np.min(self.mesh_boost.current_vtk.points[:, 0])
        min_y = np.min(self.mesh_boost.current_vtk.points[:, 1])

        # Get the length of the stimuli
        x_range = np.max(stimuli.points[:, 0]) - np.min(stimuli.points[:, 0])
        y_range = np.max(stimuli.points[:, 1]) - np.min(stimuli.points[:, 1])

        # Translate the stimuli to the bottom left corner of the mesh
        translation = (min_x - x_range, min_y - y_range, 0)
        stimuli = stimuli.translate(translation, inplace=False)

        self.stimuli_actor = self.plotter.add_mesh(
            stimuli,
            color=self.stimuli.color,
            name='stimuli',
            show_edges=False,
            smooth_shading=True,
            specular=0.8,
            metallic=0.95,
            roughness=0.0,
        )

    def draw_sensors(self):
        return
        # Draw the sensor points in the plotter
        self.sensor_actor = self.plotter.add_mesh(
            self.sensors.visualization,
            color='white',
            name='sensor_points',
            show_edges=False,
            smooth_shading=True,
            specular=0.8,
            metallic=0.99,
            roughness=0.0,
            render_points_as_spheres=True,
            opacity=0.8,
            point_size=8,
        )

    def add_mode_text(self, text):
        # Indicates which mode the user is in (Activation or Interactive)
        # Remove the text
        if self.mode_text_actor is not None:
            self.plotter.remove_actor(self.mode_text_actor)
        # Add the new text
        self.mode_text_actor = self.plotter.add_text(
            text, position='upper_right', font_size=10, color=self.text_color
        )

    def add_sensor_output_text(self, text):
        # Indicates which mode the user is in (Activation or Interactive)
        # Remove the text
        if self.sensor_output_text_actor is not None:
            self.plotter.remove_actor(self.sensor_output_text_actor)
        # Add the new text
        text = f'Sensor output: {text}'
        self.sensor_output_text_actor = self.plotter.add_text(
            text, position='upper_left', font_size=9, color=self.text_color
        )

    def add_pressure_indicator(self):
        # Add the information about the pressure that is being applied
        if self.force_indicator_actor is not None:
            self.plotter.remove_actor(self.force_indicator_actor)

        text = f'Pressure: {self.PRESSURE} N'
        self.force_indicator_actor = self.plotter.add_text(
            text, position='lower_right', font_size=8, color=self.text_color,
        )

    def increase_force(self):
        self.PRESSURE = round(min(self.PRESSURE + self.pressure_dt, 5.0), 2)
        self.add_pressure_indicator()

    def decrease_force(self):
        self.PRESSURE = round(max(self.PRESSURE - self.pressure_dt, 0.0), 2)
        self.add_pressure_indicator()

    def add_force_events(self):

        # Add key event on the right arrow press
        self.plotter.add_key_event('Right', self.increase_force)

        # Add key event on the left arrow press
        self.plotter.add_key_event('Left', self.decrease_force)

    def add_recording_actions(self):
        self.record_action = QAction('Record', self.plotter.main_menu)
        self.record_action.triggered.connect(self.start_recording)
        self.plotter.main_menu.addAction(self.record_action)

        self.stop_record_action = QAction('Stop Recording', self.plotter.main_menu)
        self.stop_record_action.triggered.connect(self.stop_recording)
        self.stop_record_action.setVisible(False)  # initially hidden
        self.plotter.main_menu.addAction(self.stop_record_action)

    def start_recording(self):
        self.recording = Recording(self.sensors, file_name=None)
        self.recording.start()
        self.update_recording_actions()

    def stop_recording(self):
        self.recording.stop()
        self.recording = None
        self.update_recording_actions()

    def update_recording_actions(self):
        if self.recording:
            self.record_action.setVisible(False)
            self.stop_record_action.setVisible(True)
        else:
            self.record_action.setVisible(True)
            self.stop_record_action.setVisible(False)

    def add_sensor_output_switch_action(self):

        # Create the QAction object
        self.sensor_output_action = QAction('Switch sensor output', self.plotter.main_menu)

        # Connect the action to the switch sensor output method
        self.sensor_output_action.triggered.connect(self.switch_sensor_output)

        # Add the action to the plotter main menu
        self.plotter.main_menu.addAction(self.sensor_output_action)

    def switch_sensor_output(self):
        # Switch the sensor output type when this method is called
        if self.sensor_output_type == 'stress':
            self.sensor_output_type = 'strain'
        else:
            self.sensor_output_type = 'stress'

        self.add_sensor_output_text(self.sensor_output_type)

        # Here you can add the code that will trigger the event when the sensor output type changes
        # For example:
        self.trigger_event(self.sensor_output_type)

        # Update the QAction text to indicate current sensor output type
        self.sensor_output_action.setText(f'Sensor output: {self.sensor_output_type}')

    def trigger_event(self, sensor_output_type):
        self.sensors.assign_sensor_outputs(sensor_output_type)

