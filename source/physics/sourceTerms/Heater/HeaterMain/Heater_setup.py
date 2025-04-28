"""
Module for executing pre-processing routines. All Unit_setup.py
modules are called from bin/unitMods.py during setup phase. The main 
routine is used to generate nucleation site and assign properties to a heater.

See Simulation/SimulationMain/incompFlow/FlowBoiling/extraParfiles/example.toml
for an example of a toml inputs file. Note that python modules interfacing
with bin/tomlParser should only have one argument (GVars)

GVars: Setup global variables
"""

import os
import numpy
import h5py
from scipy.stats import qmc  # scipy==1.8.0
from globals import GVarsClass

def main(GVars: GVarsClass):
    """Public interface for Heater setup module"""

    # At present this serves only one functionality, which
    # is to create a heater file using the dictionary encoded
    # in toml file. If this dictionary is present call the
    # private __createHeater. Otherwise continue with other tasks
    if GVars.tomlDict:
        __createHeater(GVars)
    else:
        GVars.out.put("[Heater_setup] Toml dictionary not present will not create a heater file")
     

def __createHeater(GVars):
    # Set a counter to track how many heater files are being written
    # and then loop over items in heater dictionary
    numHeaters = 0
    for key, info in GVars.tomlDict["Heater"].items():
        # if info is of type dictionary we have hit a heater configuration
        # that needs to be written to a file. Start implementing that logic
        if type(info) == dict:
            # Increase heater counter to track number of heaters
            numHeaters = numHeaters + 1

            # Raise error if the heater key does not match
            # expected naming convention
            if str(numHeaters).zfill(4) != key:
                raise ValueError(f'[Heater_setup] Heater "{key}" does not match "{str(numHeaters).zfill(4)}" in tomlfile')

            # Set filename and open the hdf5 file in write mode
            filename = os.path.join(GVars.flashHomeDir,
                                    GVars.objectDir,
                                    GVars.tomlDict["Heater"]["htr_heaterName"] + "_hdf5_htr_" + key)

            hfile = h5py.File(filename, "w")

            xsite = numpy.ndarray([info["numSites"]], dtype=float)
            ysite = numpy.ndarray([info["numSites"]], dtype=float)
            zsite = numpy.ndarray([info["numSites"]], dtype=float)
            radii = numpy.ndarray([info["numSites"]], dtype=float)

            if info["numSites"] == 1:
                xsite[:] = 0.0
                ysite[:] = 0.5*(info["ymax"] + info["ymin"])
                zsite[:] = 0.0
                radii[:] = 0.2

            else:
                halton = qmc.Halton(d=2, seed=1)
                sample = halton.random(info["numSites"])

                xsite[:] = info["xmin"] + sample[:,0]*(info["xmax"]-info["xmin"])
                ysite[:] = 0.5*(info["ymax"] + info["ymin"])
                radii[:] = 0.2

                if info["dim"] == 1:
                    zsite[:] = 0.0
                elif info["dim"] == 2:
                    zsite[:] = info["zmin"] + sample[:,1]*(info["zmax"]-info["zmin"])
                else:
                    raise ValueError(f"[Heater_setup] Error in HEATER.{key}.dim in tomlfile")

            hfile.create_dataset("heater/xMin", data=info["xmin"], shape=(1), dtype="float32")
            hfile.create_dataset("heater/xMax", data=info["xmax"], shape=(1), dtype="float32")
            hfile.create_dataset("heater/zMin", data=info["zmin"], shape=(1), dtype="float32")
            hfile.create_dataset("heater/zMax", data=info["zmax"], shape=(1), dtype="float32")
            hfile.create_dataset("heater/yMin", data=info["ymin"], shape=(1), dtype="float32")
            hfile.create_dataset("heater/yMax", data=info["ymax"], shape=(1), dtype="float32")
            hfile.create_dataset("heater/wallTemp", data=info["wallTemp"], shape=(1), dtype="float32")
            hfile.create_dataset("heater/advAngle", data=info["advAngle"], shape=(1), dtype="float32")
            hfile.create_dataset("heater/rcdAngle", data=info["rcdAngle"], shape=(1), dtype="float32")
            hfile.create_dataset("heater/velContact", data=info["velContact"], shape=(1), dtype="float32")
            hfile.create_dataset("heater/nucWaitTime", data=info["nucWaitTime"], shape=(1), dtype="float32")
            hfile.create_dataset("site/num", data=info["numSites"], shape=(1), dtype="int32")
            hfile.create_dataset("site/x", data=xsite, shape=(info["numSites"]), dtype="float32")
            hfile.create_dataset("site/y", data=ysite, shape=(info["numSites"]), dtype="float32")
            hfile.create_dataset("site/z", data=zsite, shape=(info["numSites"]), dtype="float32")
            hfile.create_dataset("init/radii", data=radii, shape=(info["numSites"]), dtype="float32")

            if GVars.setupVars.get("HeaterFluxBC"):
                hfile.create_dataset("heater/C3",data=info["C3"],shape=(1),dtype="float32")
                hfile.create_dataset("heater/C2",data=info["C2"],shape=(1),dtype="float32")
                hfile.create_dataset("heater/C1",data=info["C1"],shape=(1),dtype="float32")
                hfile.create_dataset("heater/C0",data=info["C0"],shape=(1),dtype="float32")
                hfile.create_dataset("heater/tbl_thickness",data=info["tbl_thickness"],shape=(1),dtype="float32")
                hfile.create_dataset("heater/non_uniform_temp_flag",data=info["non_uniform_temp_flag"],shape=(1),dtype="int32")
                hfile.create_dataset("heater/heat_flux_flag",data=info["heat_flux_flag"],shape=(1),dtype="int32")
                hfile.create_dataset("heater/nd_heat_flux",data=info["nd_heat_flux"],shape=(1),dtype="float32")

            hfile.close()

            GVars.out.put(f'[Heater_setup] Wrote Heater information to file '
                        + f'{filename.replace(GVars.flashHomeDir+os.sep+GVars.objectDir+os.sep,"")}')

    if numHeaters != GVars.tomlDict["Heater"]["htr_numHeaters"]:
        raise ValueError(f"[Heater_setup] Number of heater files not equal to htr_numHeaters in tomlfile")
