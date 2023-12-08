# ----- THIS REPRESENTS WORK THAT THE FLASH-X SETUP SYSTEM SHOULD EXECUTE
# Recipe builders do *not* need to know about this information and, therefore,
# it should *not* be included in recipe files.
if __init__ == "__main__":
    """
    .. todo::
        * Should we write all JSON/generated code to object folder or something
          like object/generated_code?
        * Who is assembling all of the valid group specification files needed
          for this simulation, when?
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
        * Somehow the runtime parameters for each specific call to the
          Orchestration unit (e.g., N threads to use, N blocks/packet) need to
          get inserted into the Config system.
    """
    # Assume that the following values have been set by setup tool
    FLASHX_DIR = Path("/path/to/local/flashx/clone")
    OBJECT_DIR = FLASHX_DIR.joinpath("object")
    SIMULATION_H = OBJECT_DIR.joinpath("Simulation.h")
    DESTINATION = OBJECT_DIR.joinpath("generated_code")
    TIME_ADVANCE_F90 = OBJECT_DIR.joinpath("FlashX_Name.F90")
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
                # This assumes that the group specifications needed for this
                # simulation can be correctly filled out and located (in the
                # object folder) and that the name of the subroutine as it
                # appears in the subgraph can be used by the find code to
                # identify the one and only one group specification that
                # specifies the routine.
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
    #
    # This code will need to map the graph structure of each orchestration
    # invocation onto the graph structure of the corresponding thread team
    # configuration so that it knows which Orchestration_execute* subroutine to
    # use for the invocation.
    flashx.generate_time_advance_code(
        recipe_IR, TIME_ADVANCE_F90, NO_OVERWRITE, INDENT, logger
    )
