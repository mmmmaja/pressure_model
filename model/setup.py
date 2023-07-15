from model.stimulis import *
from model._sensors import *
from model.material_handler import *
from model.mesh_converter import *

"""
Choose the stimuli from the list of stimuli and specify its size
    Sphere(radius)
    Cylinder(radius, height)
    Cube(width, length, height)
"""
_stimuli = Sphere(radius=1.6)

"""
Choose the mesh from the list of meshes and specify its features
Types of meshes:

    GridMesh(width, length, z_function, layers)
        width: width of the mesh
        length: length of the mesh
        z_function: function that returns the height of the mesh at a given point (look at mesh_helper.py)
        layers: number of layers of the mesh
        
    ArmMesh()
    
"""
_mesh_boost = GridMesh(30, 30, z_function=wave, layers=3)


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
    # Specify this
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
