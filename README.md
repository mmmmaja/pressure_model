# Pressure Model

This project is a pressure simulation model that makes use of Finite Element Analysis to simulate the 
interaction between different stimuli and a different bodies, represented by a 3D mesh. 
Pressure data in the model is recorded by the sensors in the mesh.

## Installation and Setup

To run this project, you will need the following libraries:

- PyQt5 [pip](https://pypi.org/project/PyQt5/#:~:text=PyQt5%20is%20a%20comprehensive%20set,including%20iOS%20and%20Android.)
- pyvistaqt [pip](https://pypi.org/project/pyvistaqt/)
- pyvista [pip](https://pypi.org/project/pyvista/)
- sfepy [pip](https://pypi.org/project/sfepy/)

You can install these libraries using pip:

```bash
pip install PyQt5 pyvistaqt pyvista sfepy
```

## How to use

### Recording data

Pressure data from the sensors can be record by clicking the "Record" button in the GUI.
To end the recording, click the "Stop recording" button. 
The csv file containing the data will be saved in the same directory as the project with the unique name with current date and time.

**CSV File Format**

This application records sensor data and saves them in a CSV file. The format of the CSV file is as follows:

- The first line contains the sensor positions in the format `x1 y1 z1 | x2 y2 z2 | ...`. Each `x y z` represents the 3D coordinates of a sensor.
- The rest of the file contains the sensor data in the format `stress1 | stress2 | ...`, where each `stress` value corresponds to the sensor at the same position in the first line.
- The first column represents time.

Here is an example of how the data is structured in the CSV file:

| Time | Sensor 1 (x,y,z) | Sensor 2 (x,y,z) | Sensor 3 (x,y,z) | ... |
|------|------------------|------------------|------------------|-----|
| 0    | stress           | stress           | stress           | ... |
| 1    | stress           | stress           | stress           | ... |
| 2    | stress           | stress           | stress           | ... |
| ...  | ...              | ...              | ...              | ... |

Note: Each stress value in the Sensor columns corresponds to the sensor at the same position in the first line.



### Meshes

### Stimuli

### Sensors