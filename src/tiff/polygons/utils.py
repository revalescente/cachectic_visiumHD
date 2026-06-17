from spatialdata.models import ShapesModel
from spatialdata.transformations import Identity


def add_gdf_to_sdata(sdata, gdf, key, coordinate_system):
    sdata.shapes[key] = ShapesModel.parse(
        gdf,
        transformations={coordinate_system: Identity()},
    )


# add_gdf_to_sdata(
#      sdata=sdata,
#      gdf=gdf,
#      key="nuclei_boundaries",
#      coordinate_system="blocco1",
#  )
