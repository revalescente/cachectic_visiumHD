import json
from spatialdata import SpatialData
from spatialdata.transformations import Identity, set_transformation

#
def unify_coordinate_system(sdata: SpatialData, target_cs: str) -> None:
    """
    Updates all elements in the SpatialData object to use the specified 
    target coordinate system with an Identity transformation.
    
    This operation:
      1. Iterates over every element (Images, Shapes, etc.).
      2. Removes all existing transformations.
      3. Sets a single Identity transformation to 'target_cs'.
      4. Writes changes to disk immediately if the SpatialData object is backed.
    
    Parameters
    ----------
    sdata : SpatialData
        The SpatialData object to modify.
    target_cs : str
        The name of the single coordinate system to assign to all elements.
    """
    print(f"--- Unifying coordinate systems to: '{target_cs}' ---")
    
    for _, element_name, element in sdata.gen_elements():
        print(f"Updating {element_name}...")
        
        set_transformation(
            element=element,
            transformation={target_cs: Identity()},
            set_all=True,
            write_to_sdata=sdata
        )
    
    print(f"Successfully updated all elements to '{target_cs}'.\n")

#
def read_from_json(filename):
    """Read from a given json file."""
    with open(filename) as json_file:
        data = json.load(json_file)

    return data
