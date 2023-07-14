import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import griddata


def euclidean_distance(x, y):
    """
    Calculate the euclidean distance between two points
    :param x: first point
    :param y: second point
    :return: euclidean distance between x and y
    """
    return np.sqrt(np.sum((x - y) ** 2))


def mahalanobis_distance(x, y, cov):
    """
    Calculate the mahalanobis distance between two points
    :param x: first point
    :param y: second point
    :param cov: covariance matrix
    :return: mahalanobis distance between x and y
    """
    return np.sqrt(np.dot(np.dot((x - y), np.linalg.inv(cov)), (x - y).T))


def visualize_correlation_matrix(corr_matrix):

    plt.figure(figsize=(10, 10))

    # Plot the correlation matrix for the real data
    plt.subplot(1, 1, 1)
    sns.heatmap(corr_matrix, cmap='plasma', center=0, square=True, linewidths=0.5)
    plt.title('Correlation Matrix')

    plt.tight_layout()
    plt.show()


class TimeSeriesAnalysis:

    def __init__(self, real_time_series, model_time_series):
        """
        Cross correlation analysis
            measure the similarity between two time series as a function
            of the displacement of one relative to the other
        """
        if real_time_series.shape != model_time_series.shape:
            # Shorten the longer time series
            if len(real_time_series) > len(model_time_series):
                real_time_series = real_time_series[:len(model_time_series)]
            else:
                model_time_series = model_time_series[:len(real_time_series)]

        self.real_time_series = real_time_series
        self.model_time_series = model_time_series
        print('shape of real time series: ', self.real_time_series.shape)
        print('shape of model time series: ', self.model_time_series.shape)

    def plot_3D_data(self, real_sensor_positions, model_sensor_positions, sleep_time=0.4):

        n = len(self.real_time_series)

        # Create a 3D plot
        fig = plt.figure()

        # Create first subplot
        ax1 = fig.add_subplot(1, 2, 1, projection='3d')

        # Create second subplot
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')

        # If x is clicked, the plot and app will close
        fig.canvas.mpl_connect('close_event', exit)

        # Get min and max values of Z to set consistent color scale and z-axis limits

        min_z_real = min([min(x) for x in self.real_time_series])
        min_z_model = min([min(x) for x in self.model_time_series])
        print('min_z_real: ', min_z_real)
        print('min_z_model: ', min_z_model)
        min_z = min(min_z_real, min_z_model)

        max_z_real = max([max(x) for x in self.real_time_series])
        max_z_model = max([max(x) for x in self.model_time_series])
        print('max_z_real: ', max_z_real)
        print('max_z_model: ', max_z_model)
        max_z = max(max_z_real, max_z_model)

        # Set color normalization based on overall min and max z-values
        norm = colors.Normalize(vmin=min_z, vmax=max_z)

        # Extract sensor positions
        x_values_real = np.array([float(x[0]) for x in real_sensor_positions])
        y_values_real = np.array([float(y[1]) for y in real_sensor_positions])

        x_values_model = np.array([float(x[0]) for x in model_sensor_positions])
        y_values_model = np.array([float(y[1]) for y in model_sensor_positions])

        # Create a meshgrid for X and Y
        xi_real = np.linspace(min(x_values_real), max(x_values_real), 100)
        yi_real = np.linspace(min(y_values_real), max(y_values_real), 100)
        X_real, Y_real = np.meshgrid(xi_real, yi_real)

        xi_model = np.linspace(min(x_values_model), max(x_values_model), 100)
        yi_model = np.linspace(min(y_values_model), max(y_values_model), 100)
        X_model, Y_model = np.meshgrid(xi_model, yi_model)

        for i in range(n):
            # Clear previous plot
            ax1.clear()
            ax2.clear()

            # Set consistent z-axis limits
            ax1.set_zlim([min_z, max_z])
            ax2.set_zlim([min_z, max_z])

            # Extract the sensor readings
            z_values_real = np.array([float(z) for z in self.real_time_series[i]])
            # Interpolate Z values over grid
            Z_real = griddata((x_values_real, y_values_real), z_values_real, (X_real, Y_model), method='cubic')

            z_values_model = np.array([float(z) for z in self.model_time_series[i]])
            # Interpolate Z values over grid
            Z_model = griddata((x_values_model, y_values_model), z_values_model, (X_model, Y_model), method='cubic')

            # Plot the surface on both subplots
            ax1.plot_surface(
                X_real, Y_real, Z_real, rstride=1, cstride=1, cmap='plasma', edgecolor='none', norm=norm
            )
            ax2.plot_surface(
                X_model, Y_model, Z_model, rstride=1, cstride=1, cmap='plasma', edgecolor='none', norm=norm
            )

            # Add text to the plots with the current time of the recording
            ax1.text2D(0.05, 0.95, f"Frame: {i} real", transform=ax1.transAxes)
            ax2.text2D(0.05, 0.95, f"Frame: {i} model", transform=ax2.transAxes)

            plt.draw()
            plt.pause(sleep_time)

            # Close the plot if all the frames have been shown
            if i == n - 1:
                plt.close()

        plt.show()

    def correlation_matrices(self):
        """
        Show the correlation of pressure readings between each pair of sensors in the real and simulated data
        """
        # Convert the time series to dataframes
        real_df = pd.DataFrame(self.real_time_series)
        model_df = pd.DataFrame(self.model_time_series)

        # Relationship between the readings of different sensors
        real_correlation = real_df.corr()
        model_correlation = model_df.corr()

        corr_diff = real_correlation - model_correlation

        visualize_correlation_matrix(real_correlation)
        visualize_correlation_matrix(model_correlation)
        visualize_correlation_matrix(corr_diff)

    def multivariate_distance_metrics(self):
        overall_distance = 0
        # Iterate through each time step
        for i in range(self.real_time_series.shape[0]):
            # Go through each sensor
            average_distance = 0
            for j in range(self.real_time_series.shape[1]):
                distance = euclidean_distance(
                    self.real_time_series[i][j], self.model_time_series[i][j]
                )
                average_distance += distance
            average_distance /= self.real_time_series.shape[1]
            overall_distance += average_distance
        overall_distance /= self.real_time_series.shape[0]

        print('Overall distance: ', overall_distance)
        return overall_distance
