import FlashX_RecipeTools as fr
import milhoja

import json


# TODO:
# preprocessor vars needed for group specifications
# need to gather the following information from setup.py
PP = {
    "NDIM": 2,
    "NXB": 16,
    "NYB": 16,
    "NZB": 1,
    "NGUARD": 4,
    "K1D": 1,
    "K2D": 1,
    "K3D": 0,
    "NFLUXES": 5,
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


def simpleUnsplit(recipe, root, operation_spec):
    from _nodes import simpleUnsplit_nodes

    recipe.add_operation_spec(operation_spec)

    n = simpleUnsplit_nodes()

    _tileBegin = recipe.add_item(n.tileBegin, invoke_after=root, map_to="cpu")

    _soundSpd = recipe.add_item(n.soundSpd, invoke_after=_tileBegin, map_to="gpu")
    _flx = recipe.add_item(n.flx, invoke_after=_soundSpd, map_to="gpu")
    _fly = recipe.add_item(n.fly, invoke_after=_soundSpd, map_to="gpu")
    _flz = recipe.add_item(n.flz, invoke_after=_soundSpd, map_to="gpu")
    _updSoln = recipe.add_item(n.updSoln, invoke_after=[_flx, _fly, _flz], map_to="gpu")

    _doEos = recipe.add_item(n.doEos, invoke_after=_updSoln, map_to="cpu")

    _tileEnd = recipe.add_item(n.tileEnd, invoke_after=_doEos, map_to="cpu")

    return _tileEnd


def getGridSpec(pp):
    d = {
        "dimension": pp["NDIM"],
        "nxb": pp["NXB"],
        "nyb": pp["NYB"],
        "nzb": pp["NZB"],
        "nguardcells": pp["NGUARD"]
    }
    return d


def generate_tf_specs():
    # create empty recipe
    recipe = fr.Recipe()

    # gather grid spec
    grid_spec = getGridSpec(PP)
    with open("__grid.json", "w") as f:
        json.dump(grid_spec, f, indent=2)

    # preprocess, gather, and dump operation spec
    # this will write preprocessed file to `__Hydro_op1.json`
    hydro_spec = fr.opspec.load("Hydro_op1.json", PP)

    # build recipe
    _endNode = simpleUnsplit(recipe, recipe.root, operation_spec=hydro_spec)
    _endNode = recipe.add_item(fr.LeafNode(), invoke_after=_endNode)

    # gather argument list of each nodes, required for writing `Hydro.F90`
    # TODO: perhaps redundant for generating TF specs
    recipe.traverse(controllerNode=fr.opspec.Ctr_InitNodeFromOpspec(verbose=False))

    # transform into hierarchical graph
    recipe.traverse(controllerEdge=fr.Ctr_SetupEdge(verbose=False))
    h = recipe.extractHierarchicalGraph(controllerMarkEdge=fr.Ctr_MarkEdgeAsKeep(verbose=False),
                                        controllerInitSubgraph=fr.Ctr_InitSubgraph(verbose=False))

    # generate intermediate TF data
    ctrParseTFGraph = fr.opspec.Ctr_ParseTFGraph()
    ctrParseTFNode = fr.opspec.Ctr_ParseTFNode(ctrParseTFGraph)
    ctrParseTFMultiedge = fr.opspec.Ctr_ParseTFMultiEdge(ctrParseTFGraph)
    h.traverseHierarchy(controllerGraph=ctrParseTFGraph,
                        controllerNode=ctrParseTFNode,
                        controllerMultiEdge=ctrParseTFMultiedge)
    # dump TF data for debugging
    ctrParseTFGraph.dumpTFData(indent=2)

    # Milhoja
    tfData = ctrParseTFGraph.getTFData()
    logger = milhoja.BasicLogger(level=3)
    tfSpecTpl = "tf_spec_tpl.json"   # TODO: auto generate this
    for tf in tfData:
        name = tf["name"]
        call_graph = tf["subroutine_call_graph"]
        group_json_all = tf["operation_specs"]
        grid_json = "__grid.json"
        tfAssembler = milhoja.TaskFunctionAssembler.from_milhoja_json(
                    name, call_graph, group_json_all, grid_json, logger
                 )
        tfAssembler.to_milhoja_json(f"__tf_spec_{name}.json", tfSpecTpl, overwrite=True)

