import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import griddata
from model_testing.lukas_handler import read_lukas_recording
from model_testing.model_validation import TimeSeriesAnalysis
from model_testing.pressure_image import form_contact_image
from model_testing.surface_flattening import orthographic_projection
from surface_flattening import *
import pandas as pd
import seaborn as sns
import matplotlib.ticker as ticker


def read_csv(path):
    """
    Read the data from the csv file and save it to the class variables
    :param path: path to the csv file with sensor recordings
    :return: sensor_positions, sensor_readings lists
    """
    index = 0
    sensor_positions, sensor_readings, time = [], [], []
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
                # Last column is the time
                str_row = row.split(',')
                # Convert from string to float
                float_row = [float(x) for x in str_row]
                # Convert time to seconds
                time.append(float_row[-1] / 1000)
                sensor_readings.append(float_row[:-1])
            index += 1
    return np.array(sensor_positions), np.array(sensor_readings) * (-1), np.array(time)


class ReadingManager:

    def __init__(self, sensor_positions, sensor_reading, time, mapping=orthographic_projection):
        """
        Visualizes the recording from a csv file
        """
        self.sensor_positions = sensor_positions
        self.sensor_reading = sensor_reading
        self.time = time

        self.sensor_positions_mapping = mapping(self.sensor_positions)

        # (num_measurements, num_sensors)
        print('shape of readings: ', self.sensor_reading.shape)

    def set_mapping(self, mapping):
        """
        Set the mapping from the sensor readings to the pressure image
        :param mapping: mapping from the sensor readings to the pressure image
        """
        self.sensor_positions_mapping = mapping

    def visualize(self, sleep_time=0.4, index=None):
        """
        Visualize the sensor readings with the 3D surface plot
        :param sleep_time: time to sleep between each frame
        :param index: index of the frame to visualize
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
        x_values = np.array([float(x[0]) for x in self.sensor_positions_mapping])
        y_values = np.array([float(y[1]) for y in self.sensor_positions_mapping])

        # Create a meshgrid for X and Y
        xi = np.linspace(min(x_values), max(x_values), 100)
        yi = np.linspace(min(y_values), max(y_values), 100)
        X, Y = np.meshgrid(xi, yi)

        for i in range(len(self.sensor_reading)):

            # If index is specified, only show that frame
            if index is not None:
                i = index

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

            # Add text to the plot with the current time of the recording
            ax.text2D(0.05, 0.95, f"Time: {self.time[i]} s", transform=ax.transAxes)

            plt.draw()
            plt.pause(sleep_time)

            # Close the plot if all the frames have been shown
            if i == len(self.sensor_reading) - 1 and index is None:
                plt.close()

        plt.show()

    def create_image(self, resolution=3, frame_index=None):
        if frame_index is None:
            frame_index = len(self.sensor_reading) // 2 + 2
        form_contact_image(self.sensor_positions_mapping, self.sensor_reading[frame_index], resolution)

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

    def plot_time_series(self, title=None, y_label='Sensor output'):
        sns.set_style('darkgrid')
        sns.set(rc={'figure.figsize': (10, 6)})

        # Create a dataframe with the sensor readings
        sensor_reading_df = pd.DataFrame()
        sensor_reading_df['Time'] = self.time

        # Add the column for average sensor reading for the frame
        sensor_reading_df['Average'] = np.mean(self.sensor_reading, axis=1)

        # Get the average sensor reading for each sensor
        sensor_means = np.mean(self.sensor_reading, axis=0)
        print(sensor_means)
        # Get the index of the sensor with the highest average output
        max_sensor = np.argmax(sensor_means)

        # Add the column for the sensor with the highest average output
        sensor_reading_df['Max'] = self.sensor_reading[:, max_sensor]

        # Add the column for standard deviation of the sensor readings for the frame
        sensor_reading_df['Std'] = np.std(self.sensor_reading, axis=1)

        print(sensor_reading_df)

        palette = sns.color_palette("mako_r", 3)

        ax = sns.lineplot(
            data=sensor_reading_df,
            x='Time', y='Average',
            markers=True, color=palette[0],
            dashes=False, label='Average sensor output', linewidth=2.0
        )
        # Add the standard deviation
        ax.fill_between(
            sensor_reading_df['Time'],
            sensor_reading_df['Average'] - sensor_reading_df['Std'],
            sensor_reading_df['Average'] + sensor_reading_df['Std'],
            alpha=0.2, color=palette[1], label='Standard deviation'
        )

        ax = sns.lineplot(
            data=sensor_reading_df,
            x='Time', y='Max',
            markers=True, color=palette[2],
            dashes=False, label='Sensor with maximum output', linewidth=1.5
        )

        if title is not None:
            plt.title(title)

        plt.ylabel(y_label)
        plt.xlabel('Time [sec]')
        plt.legend()
        plt.show()


def distribution_analysis(sensor_readings):

    # Get the sensor with maximum output
    max_sensor = np.argmax(np.max(sensor_readings, axis=0))
    print(f"Sensor with maximum output: {max_sensor}")

    # Plot the histogram of this sensor
    plt.hist(sensor_readings[:, max_sensor])
    plt.title(f"Histogram of sensor {max_sensor}")
    plt.show()


def get_chosen_time_frames(sensor_readings, readings_time_frame, desired_time_frames):
    """
    Get the sensor readings for the desired time frames
    Find time frames as close as possible to the desired time frames
    :param sensor_readings: (num_measurements, num_sensors)
    :param readings_time_frame:
    :param desired_time_frames:
    :return:
    """

    chosen_time_frames = []
    desired_time_frame_index = 0
    for i in range(sensor_readings.shape[0]):
        if desired_time_frame_index == len(desired_time_frames):
            break
        time_frame = readings_time_frame[i]
        if time_frame >= desired_time_frames[desired_time_frame_index]:
            chosen_time_frames.append(sensor_readings[i])
            desired_time_frame_index += 1

    cropped_recording = np.array(chosen_time_frames)
    return cropped_recording


pos_m, rec_m, time_m = read_csv('../recordings/arm_sphere_contact.csv')
pos_l, rec_l, time_l = read_lukas_recording(
    "C:/Users/majag/Desktop/arm_data/short/sphere_sensor_press_2_LIN.xlsx"
)
# cropped_rec_m = get_chosen_time_frames(rec_m, time_m, time_l)
#
# analysis = TimeSeriesAnalysis(rec_l, cropped_rec_m)
# analysis.multivariate_distance_metrics()
# # analysis.plot_3D_data(pos_l, pos_m)
# analysis.correlation_matrices()


reader_maja = ReadingManager(pos_m, rec_m, time_m)
frame_index_m = reader_maja.get_descriptive_frame()
reader_maja.plot_time_series(title="Model recording")
# reader_maja.visualize()
# reader_maja.create_image(resolution=2, frame_index=frame_index_m)


reader_lukas = ReadingManager(pos_l, rec_l, time_l)
reader_lukas.identify_faulty_sensors()
# frame_index_l = reader_lukas.get_descriptive_frame()
reader_lukas.plot_time_series(title="Robotic arm device recording")
# print(frame_index_l)
# reader_lukas.create_image(resolution=2, frame_index=reader_lukas.get_descriptive_frame())
# reader_lukas.visualize(index=frame_index_l)


# distribution_analysis(reader_maja.sensor_reading)
# distribution_analysis(reader_lukas.sensor_reading)


"""
https://www.researchgate.net/publication/224230704_Tactile-Object_Recognition_From_Appearance_Information

Characterization and simulation of tactile sensors
https://pure.johnshopkins.edu/en/publications/characterization-and-simulation-of-tactile-sensors-4
"""