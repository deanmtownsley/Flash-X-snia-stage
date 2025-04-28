"""
Module for generating monodispersed bubble distribution in 2D/3D
"""

import os
import numpy
import h5py
from scipy.stats import qmc  # scipy==1.9.0
from globals import GVarsClass


def main(GVars: GVarsClass):
    """Public interface for Simulation setup module"""
    if GVars.tomlDict:
        __createBubbleDist(GVars)
    else:
        GVars.out.put("[Simulation_setup] Toml dictionary not present will not create a bubble distribution file")


def __createBubbleDist(GVars):

    # Set filename and open the hdf5 file in write mode
    filename = os.path.join(GVars.flashHomeDir, GVars.objectDir, "bubble_dist_hdf5")

    hfile = h5py.File(filename, "w")

    xbubble = numpy.ndarray([GVars.tomlDict["Simulation"]["sim_numBubbles"]], dtype=float)
    ybubble = numpy.ndarray([GVars.tomlDict["Simulation"]["sim_numBubbles"]], dtype=float)
    zbubble = numpy.ndarray([GVars.tomlDict["Simulation"]["sim_numBubbles"]], dtype=float)

    if GVars.tomlDict["Simulation"]["sim_numBubbles"] == 1:
        xbubble[:] = 0.0
        ybubble[:] = 0.0
        zbubble[:] = 0.0

    else:
        engine = qmc.PoissonDisk(d=int(GVars.tomlDict["Simulation"]["bubble_dist"]["ndim"]), radius=0.225, seed=1)
        sample = engine.random(GVars.tomlDict["Simulation"]["sim_numBubbles"])

        if numpy.shape(sample)[0] != GVars.tomlDict["Simulation"]["sim_numBubbles"]:
            raise ValueError(f'Cannot add more than {numpy.shape(sample)[0]} bubbles. ' + 
                              'Change PoissonDisk radius to add more.')

        xbubble[:] = GVars.tomlDict["Grid"]["xmin"] + 0.5 + sample[:, 0]*(GVars.tomlDict["Grid"]["xmax"] - 
                                                                          GVars.tomlDict["Grid"]["xmin"] - 1.0)

        ybubble[:] = GVars.tomlDict["Grid"]["ymin"] + 0.5 + sample[:, 1]*(GVars.tomlDict["Grid"]["ymax"] - 
                                                                          GVars.tomlDict["Grid"]["ymin"] - 1.0)

        if GVars.tomlDict["Simulation"]["bubble_dist"]["ndim"] == 3:
            zbubble[:] = GVars.tomlDict["Grid"]["zmin"] + 0.5 + sample[:, 2]*(GVars.tomlDict["Grid"]["zmax"] - 
                                                                              GVars.tomlDict["Grid"]["zmin"] - 1.0)
        else:
            zbubble[:] = 0.0

        hfile.create_dataset("xBubble", data=xbubble, 
                             shape=(GVars.tomlDict["Simulation"]["sim_numBubbles"]), dtype="float32")

        hfile.create_dataset("yBubble", data=ybubble,
                             shape=(GVars.tomlDict["Simulation"]["sim_numBubbles"]), dtype="float32")

        hfile.create_dataset("zBubble", data=zbubble,
                             shape=(GVars.tomlDict["Simulation"]["sim_numBubbles"]), dtype="float32")
        hfile.close()

        GVars.out.put(f"[Simulation_setup] Wrote bubble distribution to file "
                    + f'{filename.replace(GVars.flashHomeDir+os.sep+GVars.objectDir+os.sep,"")}')
