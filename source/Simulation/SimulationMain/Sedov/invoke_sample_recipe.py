#!/usr/bin/env python3
"""
NOTE: This is a sample script that constructs, builds, and compiles
      a recipe for Sedov explosion problem using `Hydro` unit.
      The sample recipes are organized in `load_recipe` functions.

      This script is intended to be executed in the Flash-X object directory.

EXAMPLE: An example setup command would be,
    ```sh
    ./setup Sedov -auto -2d +sparklwf +sqr16 +mh_push +pm4dev -parfile=extraParfiles/test_amr_milhoja_cudart_spark_2d.par
    ```
    or,
    ```sh
    ./setup Sedov -auto -3d +sparklwf +cube16 +mh_push +pm4dev -parfile=extraParfiles/test_amr_milhoja_cudart_spark_3d.par
    ```
"""

import FlashX_RecipeTools as flashx
from loguru import logger
import sys

from matplotlib import pyplot as plt


def load_recipe():
    """
    A sample recipe computes a simple Spark Hydro
    with flux correction
    """

    # initialize recipe
    recipe = flashx.TimeStepRecipe()

    # assuming [ModuleName]_[Routine] format
    # JSON files will be populated by reading [ModuleName]_interface.F90
    # and feeded to milhoja_pypkg
    hydro_begin = recipe.begin_orchestration(after=recipe.root)
    hydro_prepBlock = recipe.add_work("Hydro_prepBlock", after=hydro_begin, map_to="gpu")
    hydro_advance = recipe.add_work("Hydro_advance", after=hydro_prepBlock, map_to="gpu")
    hydro_end = recipe.end_orchestration(begin_node=hydro_begin, after=hydro_advance)

    # Insert a code block after the orchestration (outside Milhoja)
    # for the flux correction. Currently, the ORCHA doesn't support
    # node-to-node communications during the orchestration; thus,
    # the flux correction algorithm is inserted as a CG-Kit template.
    #
    # recipe.add_fluxCorrection is just a thin wrapper for recipe.add_tpl,
    # which is a function that adds a node with a bare CG-Kit template.
    # recipe.add_fluxCorrection provides a dedicated template for flux correction,
    # which is maintained by FlashX_RecipeTools.
    fluxCorrection = recipe.add_fluxCorrection(after=hydro_end)

    return recipe


if __name__ == "__main__":

    # enabling logger to stdout
    logger.remove(0)
    logger.add(sys.stdout, level=0)
    logger.enable(flashx.__name__)

    # recipe
    recipe = load_recipe()

    # compile & generate codes from recipe
    ir = recipe.compile()
    ir.generate_all_codes()

    # print graphs to file
    fig = plt.figure(figsize=(16, 6))
    ax = plt.subplot(211)
    recipe.plot(nodeLabels=True, edgeLabels=True, ax=ax)
    ax = plt.subplot(212)
    ir.flowGraph.plot(nodeLabels=True, edgeLabels=True, ax=ax)
    fig.savefig("__graphs.pdf")


