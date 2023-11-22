# ----- CONTENTS OF RECIPE FILE WRITTEN BY RECIPE BUILDER
# I'm assuming for the moment that each recipe file will be a python file that
# contains the functions
# * get_variable_order and
# * load_recipe
#
# and that recipe builders will provide these with a Flash-X-defined argument
# list and in accord with an interface contract.

from pathlib import Path

from flashx import (
    CPU, GPU
    BlockTypes, BlockLevels,
    TileStepRecipe,
    SetupLogger
)

def get_variable_order(index_space, data_structure_index):
    """
    User-provided in recipe file passed to setup tool.

    .. todo::
        * What should argument list be for this function?
        * Is this an OK means to encode this information?  Will it work with
          setup tool?
        * Do we need to specify variable order for all possible Grid data
          structures?  
        * Can this be the means for expressing to the setup tool what Grid data
          structures are needed?  Or at the very least can this be used to
          generate/override the commands needed to be specified to setup tool to
          determine what Grid data structures are needed?
    """
    variable_order = None
    if index_space.lower() == "center":
        # Assume only one MultiFab for each index space for now
        assert data_structure_index == 1
        # Variable order optimized for smallest data packets.  No clue if this
        # is good for other tasks (e.g., GC filling).
        # - Variables for analytic data never go to GPU, so put at end
        # - GAME is never read from or written to.  Put at end of normal
        #   variables
        # - GAMC is read by hydro code, but is never written to.  Put just
        #   before GAME
        variable_order = [
            "DENS", "PRES", "EINT", "ENER", "TEMP",
            "VELX", "VELY", "VELZ",
            "GAMC", "GAME",
            "DENA", "PRSA", "EINA", "ENRA",
            "VLXA", "VLYA", "VLZA"
        ]
    else:
        raise ValueError(f"Unknow index space {index_space}")

    return variable_order

def load_recipe(dimension, logger):
    """
    User-provided in recipe file passed to setup tool.

    .. todo::
        * What should argument list be for this function?
        * GC fills don't belong to an operation.  They need to be created so
          that a single GC fill will have all GC data needed by all operations
          involved in the subsequent invocation of the Orchestration unit.  How
          do you specify that in the recipe?  Should we specify them in the
          recipe or let the Flash-X code generation tool insert them into the
          CG-Kit graph?
    """
    recipe = TimeStepRecipe(logger)

    GC = recipe.add_guardcell_fill("Pre-Orch_1 GC Fill", invoke_after=recipe.root)

    hydro_begin     = recipe.begin_orchestration(BlockTypes.LEAVES, BlockLevels.ALL, use_tiling=False,
                                                                               invoke_after=GC), 
    sound_speed     = recipe.add_work("Hydro_op1::Hydro_computeSoundSpeedHll", invoke_after=hydro_begin,                    map_to=GPU)
    if dimension == 2:
        # In early testing on Summit, 2D blocks do *not* provide enough work to
        # see a benefit from running compute flux routines in parallel.
        comp_flx    = recipe.add_work("Hydro_op1::Hydro_computeFluxesHll_X",   invoke_after=sound_speed,                    map_to=GPU)
        comp_fly    = recipe.add_work("Hydro_op1::Hydro_computeFluxesHll_Y",   invoke_after=comp_flx,                       map_to=GPU)
        update_soln = recipe.add_work("Hydro_op1::Hydro_updateSolutionHll",    invoke_after=comp_fly,                       map_to=GPU)
    elif dimension == 3:
        # In early testing on Summit, 3D blocks *do* provide enough work to see
        # a benefit from running compute flux routines in parallel.
        comp_flx    = recipe.add_work("Hydro_op1::Hydro_computeFluxesHll_X",   invoke_after=sound_speed,                    map_to=GPU)
        comp_fly    = recipe.add_work("Hydro_op1::Hydro_computeFluxesHll_Y",   invoke_after=sound_speed,                    map_to=GPU)
        comp_flz    = recipe.add_work("Hydro_op1::Hydro_computeFluxesHll_Z",   invoke_after=sound_speed,                    map_to=GPU)
        update_soln = recipe.add_work("Hydro_op1::Hydro_updateSolutionHll",    invoke_after=[comp_flx, comp_fly, comp_flz], map_to=GPU)
    else:
        raise ValueError("Invalid dimension for recipe")
    do_eos          = recipe.add_work("Hydro_op1::Eos_wrapped",                invoke_after=update_soln,                    map_to=CPU)
    hydro_end       = recipe.end_orchestration(hydro_begin,                    invoke_after=do_eos)

    return recipe
