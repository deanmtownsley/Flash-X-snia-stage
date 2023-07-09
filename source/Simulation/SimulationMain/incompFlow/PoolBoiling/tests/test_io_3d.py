"""Script to create input files"""

import numpy
import h5py
from scipy.stats import qmc
from types import SimpleNamespace

def createHeaterfile(heaterInfo):
    """
    Create heater file
    """
    filename = h5py.File(heaterInfo["name"] + "_hdf5_htr_0001", "w")

    xsite = numpy.ndarray([heaterInfo["numSites"]], dtype=float)
    ysite = numpy.ndarray([heaterInfo["numSites"]], dtype=float)
    zsite = numpy.ndarray([heaterInfo["numSites"]], dtype=float)
    radii = numpy.ndarray([heaterInfo["numSites"]], dtype=float)

    if heaterInfo["numSites"] == 1:
        xsite[:] = 0.0
        ysite[:] = 1e-13
        zsite[:] = 0.0
        radii[:] = 0.2

    else:

        halton = qmc.Halton(d=2, seed=1)
        haltonSample = halton.random(heaterInfo["numSites"])

        xsite[:] = heaterInfo["xmin"] + haltonSample[:, 0] * (
            heaterInfo["xmax"] - heaterInfo["xmin"]
        )
        ysite[:] = 1e-13
        radii[:] = 0.2

        if heaterInfo["dim"] == 1:
            zsite[:] = 0.0
        elif heaterInfo["dim"] == 2:
            zsite[:] = heaterInfo["zmin"] + haltonSample[:, 1] * (
                heaterInfo["zmax"] - heaterInfo["zmin"]
            )
        else:
            raise ValueError("Error in Heater.dim")

    filename.create_dataset(
        "heater/xMin", data=heaterInfo["xmin"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/xMax", data=heaterInfo["xmax"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/zMin", data=heaterInfo["zmin"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/zMax", data=heaterInfo["zmax"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/yMin", data=heaterInfo["ymin"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/yMax", data=heaterInfo["ymax"], shape=(1), dtype="float32"
    )

    filename.create_dataset(
        "heater/wallTemp", data=heaterInfo["wallTemp"], shape=(1), dtype="float32"
    )

    filename.create_dataset(
        "heater/advAngle", data=heaterInfo["advAngle"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/rcdAngle", data=heaterInfo["rcdAngle"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/velContact", data=heaterInfo["velContact"], shape=(1), dtype="float32"
    )
    filename.create_dataset(
        "heater/nucWaitTime", data=heaterInfo["nucWaitTime"], shape=(1), dtype="float32"
    )

    filename.create_dataset(
        "site/num", data=heaterInfo["numSites"], shape=(1), dtype="int32"
    )
    filename.create_dataset(
        "site/x", data=xsite, shape=(heaterInfo["numSites"]), dtype="float32"
    )
    filename.create_dataset(
        "site/y", data=ysite, shape=(heaterInfo["numSites"]), dtype="float32"
    )
    filename.create_dataset(
        "site/z", data=zsite, shape=(heaterInfo["numSites"]), dtype="float32"
    )

    filename.create_dataset(
        "init/radii", data=radii, shape=(heaterInfo["numSites"]), dtype="float32"
    )

    filename.close()

    print(f"Wrote heater information to file {filename}")


if __name__ == "__main__":

    Heater = SimpleNamespace()
  
    # Heater name and number of heaters
    Heater.name = "pool_boiling"
    Heater.dim = 2 
    Heater.xmin = -5.0
    Heater.xmax = 5.0
    Heater.zmin = -5.0
    Heater.zmax = 5.0
    Heater.ymin = 0.0
    Heater.ymax = 1e-13

    # Wall temperature
    Heater.wallTemp = 1.0

    # Contact angle information
    Heater.advAngle = 90.0
    Heater.rcdAngle = 45.0
    Heater.velContact = 0.2
    Heater.nucWaitTime = 0.2
    Heater.nucSeedRadius = 0.1

    # Number of nucleation sites
    Heater.numSites = 600

    createHeaterfile(vars(Heater))
