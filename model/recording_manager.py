import csv
import random
from datetime import datetime
from PyQt5.QtCore import QTimer

from model.fenics import FENICS
from model.pressure_script import StimuliPressure


class Recording:

    FOLDER_PATH = '../recordings'

    def __init__(self, sensors, dt=100, file_name=None):
        """
        :param sensors: Sensors object with list of sensors
        :param dt: time step in milliseconds
        :param file_name: Name of the file to save the data to (if None, then the current date and time is used)
        """

        self.sensors = sensors
        self.dt = dt
        self.file_name = file_name

        self.sensor_data = []
        self.timer = None

    def start(self):
        print("Recording Started...")
        self.timer = QTimer()
        self.timer.timeout.connect(self.record)
        self.timer.start(self.dt)  # period of dt milliseconds

    def record(self):
        self.sensor_data.append([sensor.stress for sensor in self.sensors.sensor_list])

    def stop(self):
        print("Recording Stopped...")
        self.save_data()
        self.timer.stop()
        self.timer.deleteLater()

    def save_data(self):
        """
        Save the recorded data to a .csv file

        Form of the csv file:
            First line: sensor positions  -> x y z | x y z | x y z | ...
            Rest of the file: sensor data -> stress | stress | stress | ...
            First column: time
        """

        if self.file_name is None:
            # Get current date and time
            file_name = datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
        else:
            file_name = self.file_name

        print("Saving data to: " + file_name)

        time = [i * self.dt for i in range(len(self.sensor_data))]

        # Save the data to a .csv file
        # The first line of the file should be the sensor positions
        sensor_positions = []
        for i in range(len(self.sensors.sensor_list)):
            position = self.sensors.sensor_list[i].initial_position
            sensor_positions.append(
                str(position[0]) + ' ' + str(position[1]) + ' ' + str(position[2])
            )

        # The rest of the file should be the pressure data
        path = self.FOLDER_PATH + '/' + file_name
        with open(path, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(sensor_positions)
            writer.writerows(self.sensor_data)


class ArtificialRecording(Recording):

    PRESSURE = 10

    def __init__(self, sensors, mesh_boost, rank_material, stimuli, dt=100, file_name=None):

        super().__init__(sensors, dt, file_name)

        self.mesh_boost = mesh_boost
        self.rank_material = rank_material
        self.stimuli = stimuli

        self.vertices = self.get_available_vertices()

        self.simulate_random_press()

    def get_available_vertices(self):
        """
        Find the facets of hexahedrons that are on the top face of the mesh
        :return: List of available vertices that are consist of top facets of the mesh
        """

        vertex_id_range = self.mesh_boost.top_region_ids
        for i in vertex_id_range:
            self.mesh_boost.current_vtk

        # available_vertices = []
        # # Get the cell data of the mesh
        # for cell_id in range(mesh_boost.current_vtk.n_cells):
        #     cell = mesh_boost.current_vtk.get_cell(cell_id)
        #     if self.valid_cell(cell, vertex_id_range):
        #         #
        #         available_cells.append(cell)
        # return available_cells

    def simulate_random_press(self):
        # Make a random choice of coordinates for the stimuli
        # Pick a random cell
        vertex = random.choice(self.vertices)

        coords = None

        self.stimuli.position = coords
        self.apply_pressure()
    #
    # def simulate_smooth_press(self):
    #     global ACTIVATION_THRESHOLD
    #
    #     # Number of centimeters per second
    #     # set it to random value
    #     SPEED = random.uniform(10, 20)
    #
    #     sensor_activations = np.array(self.sensor_mesh.get_data())
    #
    #     sensor_num = sensor_activations.shape[1]
    #     time_frame_num = sensor_activations.shape[0]
    #     activations = np.zeros(sensor_num)
    #     for s in range(sensor_num):
    #         for t in range(time_frame_num):
    #             if float(sensor_activations[t, s]) < ACTIVATION_THRESHOLD:
    #                 activations[s] += 1
    #
    #     current_position = np.array(self.stimuli.get_frame_position())
    #
    #     # If all activations are bigger than 3, then all sensors have been activated at least three times
    #     # Then return
    #     if np.all(activations >= 3):
    #         # If all sensors have been activated at least three times, continue in the same direction
    #         target_position = current_position + self.direction * SPEED
    #     else:
    #         # Find the index of the sensor that was activated the least
    #         index = np.argmin(activations)
    #         target_position = np.array(self.sensor_mesh.SENSOR_ARRAY[index].frame_position)
    #         target_position = np.insert(target_position, 2, 0)
    #
    #     # Calculate the direction vector from the stimulus position to the target
    #     direction_vector = target_position - current_position
    #     # Normalize the direction vector
    #     direction_vector /= np.linalg.norm(direction_vector)
    #
    #     # Add a random component to the direction vector
    #     offset = 0.1
    #     direction_vector += np.array([
    #         random.uniform(-offset, offset), random.uniform(-offset, offset), 0])
    #
    #     # # Displacement in this iteration is dependent on the update interval
    #     local_displacement = SPEED
    #
    #     # If the stimulus hits the frame boundary, change direction
    #     if not 0 <= current_position[0] + local_displacement * direction_vector[0] <= self.frame_dim[0] or \
    #             not 0 <= current_position[1] + local_displacement * direction_vector[1] <= self.frame_dim[1]:
    #         direction_vector *= -1
    #
    #     self.direction = direction_vector
    #
    #     # Modify the position to move towards the target sensor
    #     position = [
    #         current_position[0] + local_displacement * direction_vector[0],
    #         current_position[1] + local_displacement * direction_vector[1],
    #         0]
    #
    #     self.stimuli.set_frame_position(position)


    def apply_pressure(self):
        force_handler = StimuliPressure(self.stimuli, self.PRESSURE)
        fenics = FENICS(self.mesh_boost.sfepy_mesh, self.rank_material, self.sensors)

        # Calculate the displacement
        fenics.apply_pressure(force_handler)