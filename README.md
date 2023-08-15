# Pressure Model

This project is a pressure simulation model that makes use of Finite Element Analysis to simulate the 
interaction between different stimuli and a different bodies, represented by a 3D mesh. 
Pressure data in the model is recorded by the sensors in the mesh.

Model is implemented in the folder **model**. Additional model testing files are located in the folder **model_testing**.

## Installation and Setup

To run this project, you will need the following libraries:

- PyQt5 [pip](https://pypi.org/project/PyQt5/#:~:text=PyQt5%20is%20a%20comprehensive%20set,including%20iOS%20and%20Android.)
- pyvistaqt [pip](https://pypi.org/project/pyvistaqt/)
- pyvista [pip](https://pypi.org/project/pyvista/)
- sfepy [pip](https://pypi.org/project/sfepy/)
- numpy [pip](https://pypi.org/project/numpy/)
- pandas [pip](https://pypi.org/project/pandas/)
- matplotlib [pip](https://pypi.org/project/matplotlib/)
- scipy [pip](https://pypi.org/project/scipy/)
- vtk [pip](https://pypi.org/project/vtk/)

You can install these libraries using pip:

```bash
pip install PyQt5 pyvistaqt pyvista sfepy numpy pandas matplotlib scipy vtk
```

## How to use

### Simulation Setup
    
You can choose the objects that will be used in the simulation by modifying the **setup.py** file in the model folder.
This includes the following:
- Stimuli
- Mesh
- Sensors
- Material properties

(All of these are explained in more detail in the **setup.py** file)

### Simulation modes

There are two simulation modes available in this project:
- **Interactive Mode** - In this mode user can rotate and inspect the object within the simulation space.
- **Activation Mode** - This mode enables user to apply the force onto the mesh with selected stimuli.

Modes are changes by pressing the control key.

### Recording data

Pressure data from the sensors can be record by clicking the "Record" button in the GUI.
To end the recording, click the **Stop recording** button. 
The csv file containing the data will be saved in the same directory as the project in the folder recordings with the unique name with current date and time.


### Choosing sensor output type

The sensor output type can switch by clicking the **Switch sensor output** button in the GUI.
Current sensor output type is displayed in the GUI.
When recording the data the sensor output type will be added to the CSV file name.

### CSV File Format

This application records sensor data and saves them in a CSV file. The format of the CSV file is as follows:

- The first line contains the sensor positions in the format `x1 y1 z1 | x2 y2 z2 | ...`. Each `x y z` represents the 3D coordinates of a sensor.
- The rest of the file contains the sensor data in the format `sensor_output1 | sensor_output2 | ...`, where each `sensor_output` value corresponds to the sensor at the same position in the first line.
- The last column represents time.

Here is an example of how the data is structured in the CSV file:

| Sensor 1 (x,y,z) | Sensor 2 (x,y,z) | Sensor 3 (x,y,z) | ... | Time   |
|------------------|------------------|------------------|-----|--------|
| sensor_output    | sensor_output    | sensor_output    | ... | time0  |
| sensor_output    | sensor_output    | sensor_output    | ... | time1  |
| sensor_output    | sensor_output    | sensor_output    | ... | time2  |
| ...              | ...              | ...              | ... | ...    |

Note: Each sensor output value in the Sensor columns corresponds to the sensor at the same position in the first line.




