from flashx import (
    CPU, GPU
    BlockTypes, BlockLevels,
    TileStepRecipe
)

def get_variable_order(index_space, data_structure_index):
    """
    .. todo::
        * Do we need to specify variable order for all possible Grid data
          structures?  
        * Can this be the means for expressing to the setup tool what Grid data
          structures are needed?  Or at the very least can this be used to
          generate/override the commands needed to be specified to setup tool to
          determine what Grid data structures are needed?
    """
    variable_order = {}
    if index_space.lower() == "CENTER":
        assert data_structure_index == 1
        variable_order = {
            "DENS_VAR": 1,
            "PRES_VAR": 2,
            "VELX_VAR": 3,
            "VELY_VAR": 4,
            "VELZ_VAR": 5,
            "ENER_VAR": 6,
            "EINT_VAR": 7,
            "TEMP_VAR": 8,
            "GAMC_VAR": 9,
        }
    else:
        raise ValueError(f"Unknow index space {index_space}")

    return variable_order

def generate_simple_unsplit_recipe_2D(verbose=False):
    """
    .. todo::
        * Is it reasonable to assume that each of the Static Fortran Routines
          has a unique name in the sense that it will appear in one and only one
          operation spec?  Seems OK for the Hydro, but what about Eos_wrapper?
    """
    recipe = TimeStepRecipe(verbose=verbose)

    hydro_begin = recipe.begin_orchestration(BlockTypes.LEAVES, BlockLevels.ALL, use_tiling=False)
    sound_speed = recipe.add_work("Hydro_computeSoundSpeedHll", invoke_after=hydro_begin, map_to=GPU)
    comp_flx    = recipe.add_work("Hydro_computeFluxesHll_X",   invoke_after=sound_speed, map_to=GPU)
    comp_fly    = recipe.add_work("Hydro_computeFluxesHll_Y",   invoke_after=comp_flx,    map_to=GPU)
    update_soln = recipe.add_work("Hydro_updateSolutionHll",    invoke_after=comp_fly,    map_to=GPU)
    do_eos      = recipe.add_work("Eos_wrapped",                invoke_after=update_soln, map_to=CPU)
    hydro_end   = recipe.end_orchestration(hydro_begin,         invoke_after=do_eos)

    return recipe

def generate_simple_unsplit_recipe_3D(verbose=False):
    """
    .. todo::
        * Is it reasonable to assume that each of the Static Fortran Routines
          has a unique name in the sense that it will appear in one and only one
          operation spec?  Seems OK for the Hydro, but what about Eos_wrapper?
    """
    recipe = TimeStepRecipe(verbose=verbose)

    hydro_begin = recipe.begin_orchestration(BlockTypes.LEAVES, BlockLevels.ALL, use_tiling=False)
    sound_speed = recipe.add_work("Hydro_computeSoundSpeedHll", invoke_after=hydro_begin,     map_to=GPU)
    comp_flx    = recipe.add_work("Hydro_computeFluxesHll_X",   invoke_after=sound_speed,     map_to=GPU)
    comp_fly    = recipe.add_work("Hydro_computeFluxesHll_Y",   invoke_after=sound_speed,     map_to=GPU)
    comp_flz    = recipe.add_work("Hydro_computeFluxesHll_Z",   invoke_after=sound_speed,     map_to=GPU)
    update_soln = recipe.add_work("Hydro_updateSolutionHll",    invoke_after=[flx, fly, flz], map_to=GPU)
    do_eos      = recipe.add_work("Eos_wrapped",                invoke_after=update_soln,     map_to=CPU)
    hydro_end   = recipe.end_orchestration(hydro_begin,         invoke_after=do_eos)

    return recipe
