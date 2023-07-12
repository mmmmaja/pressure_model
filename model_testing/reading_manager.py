import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import griddata

from model_testing.lukas_handler import read_lukas_recording
from model_testing.pressure_image import form_contact_image


def read_csv(path):
    """
    Read the data from the csv file and save it to the class variables
    :param path: path to the csv file with sensor recordings
    :return: sensor_positions, sensor_readings lists
    """
    index = 0
    sensor_positions, sensor_readings = [], []
    with open(path, 'r', newline='') as file:
        for row in file.readlines():
            # First row contains position of the sensors
            if index == 0:
                sensor_positions = row.split(',')
                # Convert from string to float
                for i in range(len(sensor_positions)):
                    sensor_positions[i] = sensor_positions[i].split(' ')
                    sensor_positions[i] = [float(x) for x in sensor_positions[i]]
            # The rest of the rows contain the sensor readings
            else:
                str_row = row.split(',')
                # Convert from string to float
                float_row = [float(x) for x in str_row]
                sensor_readings.append(float_row)
            index += 1
    return np.array(sensor_positions), np.array(sensor_readings) * (-1)


class ReadingManager:

    def __init__(self, sensor_positions, sensor_reading):
        """
        Visualizes the recording from a csv file
        """
        self.sensor_positions = sensor_positions
        self.sensor_reading = sensor_reading
        # (num_measurements, num_sensors)
        print('shape of readings: ', self.sensor_reading.shape)

    def visualize(self, sleep_time=0.4):
        """
        Visualize the sensor readings with the 3D surface plot
        :param sleep_time: time to sleep between each frame
        """

        # Create a 3D plot
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        # If x is clicked, the plot and app will close
        fig.canvas.mpl_connect('close_event', exit)

        # Get min and max values of Z to set consistent color scale and z-axis limits
        min_z = min([min(x) for x in self.sensor_reading])
        max_z = max([max(x) for x in self.sensor_reading])

        # Set color normalization based on overall min and max z-values
        norm = colors.Normalize(vmin=min_z, vmax=max_z)

        # Extract sensor positions
        x_values = np.array([float(x[0]) for x in self.sensor_positions])
        y_values = np.array([float(y[1]) for y in self.sensor_positions])

        # Create a meshgrid for X and Y
        xi = np.linspace(min(x_values), max(x_values), 100)
        yi = np.linspace(min(y_values), max(y_values), 100)
        X, Y = np.meshgrid(xi, yi)

        for i in range(len(self.sensor_reading)):
            # Clear previous plot
            ax.clear()

            # Set consistent z-axis limits
            ax.set_zlim([min_z, max_z])

            # Extract the sensor readings
            z_values = np.array([float(z) for z in self.sensor_reading[i]])
            # Interpolate Z values over grid

            Z = griddata((x_values, y_values), z_values, (X, Y), method='cubic')

            # Plot the surface
            ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='plasma', edgecolor='none', norm=norm)

            # Add text to the plot with the current iteration
            ax.text2D(0.05, 0.95, f'Frame {i}', transform=ax.transAxes)

            plt.draw()
            plt.pause(sleep_time)

            # Close the plot if all the frames have been shown
            if i == len(self.sensor_reading) - 1:
                plt.close()

        plt.show()

    def create_image(self, resolution=3, frame_index=None):
        if frame_index is None:
            frame_index = len(self.sensor_reading) // 2 + 2
        form_contact_image(self.sensor_positions, self.sensor_reading[frame_index], resolution)

    def identify_faulty_sensors(self, threshold_std=1e-5, threshold_mean_factor=5):

        # shape of readings: (num_measurements, num_sensors)

        # Get the mean and std for all sensors
        means = np.mean(self.sensor_reading, axis=0)
        stds = np.std(self.sensor_reading, axis=0)

        # Get the overall mean and std
        overall_mean = np.mean(means)
        overall_std = np.std(means)

        # Iterate through readings of each sensor
        for i in range(self.sensor_reading.shape[1]):
            # Get the mean and the standard deviation of the readings for each sensor
            mean = means[i]
            std = stds[i]

            # Identify faulty sensors based on std, mean, and outliers
            # Flag the sensor as faulty
            if std < threshold_std and abs(mean - overall_mean) > threshold_mean_factor * overall_std:
                # Change the output of the sensor at each time step to the overall mean
                self.sensor_reading[:, i] = overall_mean
                print(f"Sensor {i}, mean: {mean}, std: {std}")
                print(f"Sensor {i} is faulty and has been replaced with the overall mean")

    def get_descriptive_frame(self):
        # Find a frame with the highest output in sensor readings
        max_frame = np.argmax(np.max(self.sensor_reading, axis=1))
        return max_frame

    def get_statistics(self):

        # Get the mean and std for all sensors
        means = np.mean(self.sensor_reading, axis=0)
        stds = np.std(self.sensor_reading, axis=0)

        # Get the overall mean and std
        overall_mean = np.mean(means)
        overall_std = np.std(means)

        # Get min and max pressure
        min_pressure = np.min(self.sensor_reading)
        max_pressure = np.max(self.sensor_reading)

        print(f"Overall mean: {overall_mean}")
        print(f"Overall std: {overall_std}")
        print(f"Min pressure: {min_pressure}")
        print(f"Max pressure: {max_pressure}")


def distribution_analysis(sensor_readings):

    # Get the sensor with maximum output
    max_sensor = np.argmax(np.max(sensor_readings, axis=0))
    print(f"Sensor with maximum output: {max_sensor}")

    # Plot the histogram of this sensor
    plt.hist(sensor_readings[:, max_sensor])
    plt.title(f"Histogram of sensor {max_sensor}")
    plt.show()



pos, rec = read_csv('../recordings/recording.csv')
reader_maja = ReadingManager(pos, rec)
# reader_maja.get_statistics()
# reader_maja.visualize()


pos, rec = read_lukas_recording(
    "C:/Users/majag/Desktop/arm_data/sphere/1_LIN.xlsx"
)
reader_lukas = ReadingManager(pos, rec)
reader_lukas.identify_faulty_sensors()
# reader_lukas.get_statistics()
reader_lukas.create_image()

distribution_analysis(reader_maja.sensor_reading)
distribution_analysis(reader_lukas.sensor_reading)