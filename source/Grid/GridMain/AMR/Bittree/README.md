# What is Bittree?
When FLASH uses Paramesh for its underlying mesh structure, the user may find that their program incurs penalties in a massively parallel environment. This is because Paramesh distributes the tree metadata across nodes, and any particular node does not have access to the global structure. Bittree provides a memory efficient way of storing some global metadata on each processor, which can improve scaling behavior.

From John Bachan's paper: "The bittree is a memory efficient data structure for storing the entire octree structure of the mesh on each node. It also yields an efficient algorithm for mapping from a block's coordinates in the domain to its index along the space filling Morton curve which is used to distribute blocks among processors. Together, these features enable determination of the processor to which a given block belongs without requiring any off-node communication."

Currently, Bittree is primarily implemented to speed up the regridding process. The native Paramesh algorithm makes a huge number of MPI calls, including costly Allreduce's, throughout every step. With the global view of the tree provided by Bittree, this algorithm is heavily simplified and many rounds of communication are cut entirely. For details, see the alternative implementation of Paramesh files provided in this directory.
# How to use Bittree
In order to use Paramesh with Bittree enabled, include `Bittree=True` in your setup call. Note that Bittree does not currently support all capabilities of FLASH, notably 2d and 3d spherical coordinates.

Here is a useful setup configuration from the test suite that can be used to test that Bittree is running correctly. This configuration tests many aspects of Bittree, such as multiple root blocks, frequent derefinement while running, and the alternative Morton order.

`./setup DustCollapse -auto -2d +cylindrical +Mode1 +serialIO +uhd +newMpole -parfile=test_2dcyl.par AltMorton=True Bittree=True`

Note: Bittree's C++ source code requires a standard library that may not be automatically linked by your compiler. If you run into problems with the stanard libray not linking correctly, try putting the following line in your Makefile:
`LIB_STDCXX = -lc++ -lstdc++`
# Runtime Parameters
The user can control the behavior of Bittree with a runtime parameter. 
- `gr_btDistributedSort` - If using custom sort-by-work (i.e. gr\_btCustomWork=True), this controls behavior of the sorting algorithm. True = Uses less local storage but increases communication. False (default) = Uses less communication but high memory consumption.
# Implementation
First initialize the tree with the root block configuration. From there, at each refinement step update the Bittree and use it to calculate tree data, reducing MPI calls and hopefully speeding up the refinement process.

Outside of refinement, Bittree is accessed to locate a block (identify its proc ID and local block number) given its xyz coordinates and refinement level. Specifically, the gr\_xyzToBlock routine, (usually called from Grid\_getBlkIDFromPos), uses Bittree if possible and inter-processor communication if not. 
### How to refine/derefine the tree:
   1. `call bittree_refine_init`
   2. `call gr_btRefineMark`  for every local block thats going to be refined.
   3. `call bittree_refine_reduce`  once all processors have marked local blocks. This propagates local refinement patterns across all processors.
   4. `call bittree_refine_update`  once refinement pattern is settled. This creates a new tree alongside the old tree.
   5. `call bittree_refine_apply` once new tree data is calculated. This replaces the old tree with the new.
   6. `call gr_btSortMortonBittree`  to update localMortUB, allowing for blocks to be identified.

After that the tree corresponds to the new domain. (This is more complicated if any blocks need to be derefined. See amr\_refine\_derefine and amr\_check\_derefine for details.)

# Contents of Bittree package
### Wrapper routines
The follow routines each have their own file, and mostly serve as Fortran wrapper routines that take standard Fotran types (integer, logical) as input/output.
- gr\_btDerefineMark - marks a leaf block for derefinement (technially marks *parent* for nodetype change).
- gr\_btGetBitid - returns a block's bitid (its location in the Bittree's array of blocks).
- gr\_btGetDerefinee - returns whether a leaf block is marked for derefinement (technially whether *parent* is marked for nodetype change).
- gr\_btGetLocalBitids - given a range of global morton indices, returns an array of bitids blocks in that range
- gr\_btGetRefine - returns whether a leaf block is marked for refinement.
- gr\_btIsParent - returns whether a block is a parent or a leaf
- gr\_btLocate - given a block's bitid, return its location and refinement level
- gr\_btRefineMark - marks a leaf block for refinement
- gr\_btIdentify - routine that identifies a block's proc ID and local block number based on its location. Must be run after amr\_sort\_morton\_bittree updates morton index cutoffs.
### Management of Bittree
These are routines that create and verify Bittree. Need to be implemented by the Amr package.
- amr\_build\_bittree
- amr\_verify\_bittree
### Core Bittree functions (in c\_source/)
Most of the core functions that interact with Bittree's C++ source code are defined in `bittree_core.cxx` and interfaced in `bittree.F90`. `bittree_core.cxx` also defines all of the functions that actually create, modify, and access the Bittree-related arrays, bitmaps, variables, etc.
### C++ source code (in c\_source/)
The remaining files in the Bittree package are the C++ source files that actually manage the nitty-gritty details of the Bittree. See individual files for more documentation on their functionality and usage. For most purposes though, the FLASH user should consider Bittree a "block box" accessed and managed through the Fortran routines and wrappers. 
