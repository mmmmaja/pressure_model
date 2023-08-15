from stimulis import *
from _sensors import *
from material_handler import *
from mesh_converter import *

"""
In this file you can specify the features of the model
"""

"""
Choose the stimuli from the list of stimuli and specify its size
    Sphere(radius)
    Cylinder(radius, height)
    Cube(width, length, height)
"""
_stimuli = Cylinder(radius=3.6, height=0.5)

"""
Choose the mesh from the list of meshes and specify its features
Types of meshes:
    
    (Mesh defined by the user)
    GridMesh(width, length, z_function, layers)
        width: width of the mesh
        length: length of the mesh
        z_function: function that returns the height of the mesh at a given point (look at mesh_helper.py)
        layers: number of layers of the mesh
    
    (3D model of the robotic arm) 
    ArmMesh()
    
"""
# _mesh_boost = ArmMesh()
_mesh_boost = GridMesh(30, 30, z_function=flat, layers=3)


if _mesh_boost is ArmMesh:
    _sensors = SensorArm(_mesh_boost)
else:
    """
    Choose the sensors from the list of sensors and specify its features
    Types of sensors:
        SensorGrid(width, length, mesh_boost)
            width: width of the sensor grid
            length: length of the sensor grid
            
        RandomSensors(number_of_sensors, mesh_boost)
            number_of_sensors: number of sensors in the grid
            
        SensorPatchesFromFile(file_name, mesh_boost)
            file_name: name of the file with the sensor patches
            
    """
    # Specify params
    _sensors = SensorGrid(10, 10, _mesh_boost)


"""
Choose the material for the mesh
(For more details look at the material_handler.py)

    rubber
    steel
    silicon
    polyurethane foam
"""
_material = rubber
