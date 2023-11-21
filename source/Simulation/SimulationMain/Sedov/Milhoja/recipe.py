# ----- CONTENTS OF RECIPE FILE WRITTEN BY RECIPE BUILDER
# I'm assuming for the moment that each recipe file will be a python file that
#
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
        * Is it reasonable to assume that each of the Static Fortran Routines
          has a unique name in the sense that it will appear in one and only one
          operation spec?  Seems OK for the Hydro, but what about Eos_wrapped?
    """
    recipe = TimeStepRecipe(logger)

    hydro_begin     = recipe.begin_orchestration(BlockTypes.LEAVES, BlockLevels.ALL, use_tiling=False)
    sound_speed     = recipe.add_work("Hydro_computeSoundSpeedHll", invoke_after=hydro_begin,     map_to=GPU)
    if dimension == 2:
        # In early testing on Summit, 2D blocks do *not* provide enough work to
        # see a benefit from running compute flux routines in parallel.
        comp_flx    = recipe.add_work("Hydro_computeFluxesHll_X",   invoke_after=sound_speed,     map_to=GPU)
        comp_fly    = recipe.add_work("Hydro_computeFluxesHll_Y",   invoke_after=comp_flx,        map_to=GPU)
        update_soln = recipe.add_work("Hydro_updateSolutionHll",    invoke_after=comp_fly,        map_to=GPU)
    elif dimension == 3:
        # In early testing on Summit, 3D blocks *do* provide enough work to see
        # a benefit from running compute flux routines in parallel.
        comp_flx    = recipe.add_work("Hydro_computeFluxesHll_X",   invoke_after=sound_speed,     map_to=GPU)
        comp_fly    = recipe.add_work("Hydro_computeFluxesHll_Y",   invoke_after=sound_speed,     map_to=GPU)
        comp_flz    = recipe.add_work("Hydro_computeFluxesHll_Z",   invoke_after=sound_speed,     map_to=GPU)
        update_soln = recipe.add_work("Hydro_updateSolutionHll",    invoke_after=[flx, fly, flz], map_to=GPU)
    else:
        raise ValueError("Invalid dimension for recipe")
    do_eos          = recipe.add_work("Eos_wrapped",                invoke_after=update_soln,     map_to=CPU)
    hydro_end       = recipe.end_orchestration(hydro_begin,         invoke_after=do_eos)

    return recipe

# ----- THIS REPRESENTS WORK THAT THE FLASH-X SETUP SYSTEM SHOULD EXECUTE
# Recipe builders do *not* need to know about this information and, therefore,
# it should *not* be included in recipe files.
if __init__ == "__main__":
    """
    .. todo::
        * Should we write all JSON/generated code to object folder or something
          like object/generated_code?
        * Keep track of all files that are being generated and where they are
          written to.
        * Use the TaskFunction objects to collect information such as
          constructor argument list that are needed to generated time advance
          code.  Perhaps the set of TaskFunction objects or their JSON files
          should be passed to the time advance generator code.
        * Once all code has been generated, write to some set of files
          information about all generated code that needs compiling.  We assume
          that the Flash-X build system will use that information to compile all
          files needed to build and link the simulation binary.
    """
    # Assume that the following values have been set by setup tool
    FLASHX_DIR = Path("/path/to/local/flashx/clone")
    OBJECT_DIR = FLASHX_DIR.joinpath("object")
    SIMULATION_H = OBJECT_DIR.joinpath("Simulation.h")
    DESTINATION = OBJECT_DIR.joinpath("generated_code")
    TIME_ADVANCE_F90 = DESTINATION.joinpath("FlashX_Name.F90")
    GRID_JSON = DESTINATION.joinpath("grid.json")
    MILHOJA_PATH = Path("/path/from/site/makefile")
    NO_OVERWRITE = False
    INDENT = 4

    # Assume that verbosity have been set by setup tool using values passed in
    # by user to setup command
    verbosity = LOG_LEVEL_NONE, LOG_LEVEL_BASIC, ..., or LOG_LEVEL_MAX
    logger = SetupLogger(verbosity)

    # Write grid specification to file
    with open(GRID_JSON, "w") as fptr:
        flashx.write_grid_json(fptr, NO_OVERWRITE)
        grid_spec = json.load(GRID_JSON)
    ndim = grid_spec["dimension"]

    # Write Simulation.h using variable ordering function in recipe file
    # provided to setup tool
    flashx.write_simulation_h(SIMULATION_H, 
                              get_variable_order,
                              NO_OVERWRITE,
                              logger)

    # Get recipe from recipe file provided to setup tool
    recipe = load_recipe(ndim, logger)

    # Compile recipe into intermediate representation
    recipe_IR = flashx.compile_recipe(recipe, DESTINATION, NO_OVERWRITE, logger)

    # Generate code for Orchestration unit
    for subgraph in recipe_IR:
        if subgraph.use_orchestration:
            # Construct task function & write its specification
            tf_name = subgraph.name
            tf_full_json = f"{tf_name}.json"
            tf_partial_json = f"{tf_name}_partial.json"

            # Write partial specification to file using metadata associated with
            # subgraph in IR by recipe compiler.
            data_item = subgraph.data_item_type
            tf_partial_spec = {
                "task_function": {
                    "language":       "Fortran",
                    "processor":      subgraph.processor,
                    "cpp_header":     f"{tf_name}.h",
                    "cpp_source":     f"{tf_name}.cxx",
                    "c2f_source":     f"{tf_name}_C2F.F90",
                    "fortran_source": f"{tf_name}.F90"
                },
                "data_item": {
                    "type":   data_item,
                    "header": f"{data_item}_{tf_name}.h",
                    "source": f"{data_item}_{tf_name}.cxx"
                }
            }
            with open(tf_partial_json, "w") as fptr:
                json.dump(tf_partial_spec, fptr) 

            # Find the group specification files that specify the subroutines to
            # be included in the task function.
            group_specs = []
            for subroutine in subgraph:
                group_specs.append(flashx.find_group_specification(subroutine))

            assembler = milhoja.TaskFunctionAssembler.from_milhoja_json(
                            tf_name, subgraph,
                            group_specs, grid_json, logger
                        )
            assembler.to_milhoja_json(tf_full_json, tf_partial_json, NO_OVERWRITE)

            # TODO: Run optimization methods over TF specification if users
            # specify that they want optimizations.  NO SUCH OPTIMIZATION
            # METHODS EXIST NOR WILL THEY BE CREATED ANYTIME SOON.

            # Write task function's code for use with Orchestration unit
            tf_spec = milhoja.TaskFunction.from_milhoja_json(tf_full_json, logger)
            milhoja.generate_data_item(
                tf_spec, DESTINATION, NO_OVERWRITE, MILHOJA_PATH, INDENT, logger
            )
            milhoja.generate_task_function(
                tf_spec, DESTINATION, NO_OVERWRITE, INDENT, logger
            )

    # Generate code for TimeAdvance
    flashx.generate_time_advance_code(
        recipe_IR, TIME_ADVANCE_F90, NO_OVERWRITE, INDENT, logger
    )
