import csv
from datetime import datetime
from PyQt5.QtCore import QTimer

"""
This script contains the Recording class, which is used to record the data from the sensors.
"""


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
        self.sensor_data.append([
            str(sensor.output[0]) + ' ' + str(sensor.output[1]) + ' ' + str(sensor.output[2])
            for sensor in self.sensors.sensor_list
        ])

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
            Rest of the file: sensor data -> output (x y z) | output (x y z) | output (x y z) | ...
            LAST column: time
        """

        if self.file_name is None:
            # Get current date and time
            file_name = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        else:
            file_name = self.file_name

        # Add the indication of the sensor output type to the file name
        file_name += '_' + self.sensors.output_type + '.csv'

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
            # In the first line, write the sensor positions
            writer.writerow(sensor_positions)

            for i in range(len(self.sensor_data)):
                # Get all the data from the sensors in the i-th frame and write it to the file as strings
                writer.writerow(self.sensor_data[i] + [str(time[i] / 1000)])

            print("Saving data to: " + file_name)
