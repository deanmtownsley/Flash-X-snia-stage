#include "bittree_bitarray.hxx"
#include "bittree_mortontree.hxx"
#include "bittree_ref.hxx"
#include "mpi.h"
#include <iostream>

using namespace BitTree;

typedef unsigned W;

namespace { // private globals
  struct TheTree_ {
    virtual unsigned block_count(bool updated) = 0;
    virtual unsigned leaf_count(bool updated) = 0;
    virtual unsigned delta_count() = 0;
    virtual bool check_refine_bit(int bitid) = 0;
    virtual bool is_parent(bool updated, int bitid) = 0;
    virtual void identify(bool updated, int *lev, int *ijk, int *mort, int *bitid) = 0;
    virtual void locate(bool updated, int bitid, int *lev, int *ijk, int *mort) = 0;
    virtual void get_bitid_list(bool updated, int mort_min, int mort_max, int *out) = 0;
    virtual void refine_init() = 0;
    virtual void refine_mark(int bitid, bool value) = 0;
    virtual void refine_reduce(MPI_Comm comm) = 0;
    virtual void refine_reduce_and(MPI_Comm comm) = 0;
    virtual void refine_update() = 0;
    virtual void refine_apply() = 0;
    virtual void print_2d(int datatype=0) = 0;
  };
  
  template<unsigned ndim>
  class TheTree: public TheTree_ {
    Ref<MortonTree<ndim,W> > tree;               //Actual Bittree
    Ref<MortonTree<ndim,W> > tree_updated;       //Updated Bittree, before refinement is applied
    Ref<BitArray<W> > refine_delta;              //(De)refinement flags for blocks
    bool is_reduced;                             //Flag to track whether refine_delta is up to date across processors
    bool is_updated;                             //Flag to track whether tree_updated matches latest refine_delta
    bool in_refine;                              //If in_refine=false, tree_updated and refine_delta should not exist 
  public:
    TheTree(const unsigned top[], const bool includes[]);
    unsigned block_count(bool updated);
    unsigned leaf_count(bool updated);
    unsigned delta_count();
    bool check_refine_bit(int bitid);
    bool is_parent(bool updated, int bitid);
    void identify(bool updated, int *lev, int *ijk, int *mort, int *bitid);
    void locate(bool updated, int bitid, int *lev, int *ijk, int *mort);
    void get_bitid_list(bool updated, int mort_min, int mort_max, int *out);
    void refine_init();
    void refine_mark(int bitid, bool value);
    void refine_reduce(MPI_Comm comm);
    void refine_reduce_and(MPI_Comm comm);
    void refine_update();
    void refine_apply();
    void print_2d(int datatype=0);
  };
  
  Ref<TheTree_> the_tree;
}

/** Constructor for TheTree */
template<unsigned ndim>
TheTree<ndim>::TheTree(const unsigned top[], const bool includes[]):
  tree(MortonTree<ndim,W>::make(top, includes)),
  is_reduced(false),
  is_updated(false),
  in_refine(false)  {
}

/** Check block count of Bittree. Wrapper for MortonTree's blocks() */
template<unsigned ndim>
unsigned TheTree<ndim>::block_count(bool updated) {
  if(updated && in_refine) {
    if (not is_updated) refine_update();
    return tree_updated->blocks();
  }
  else
    return tree->blocks();
}

/** Check block count of Bittree. Wrapper for MortonTree's leaves() */
template<unsigned ndim>
unsigned TheTree<ndim>::leaf_count(bool updated) {
  if(updated && in_refine) {
    if (not is_updated) refine_update();
    return tree_updated->leaves();
  }
  else
    return tree->leaves();
}
/** Check number of blocks marked for nodetype change */
template<unsigned ndim>
unsigned TheTree<ndim>::delta_count() {
  if(in_refine) return refine_delta->count();
  else return 0;
}

/** Check refinement bit. Wrapper for BitArray's get() */
template<unsigned ndim>
bool TheTree<ndim>::check_refine_bit(int bitid) {
  if(in_refine) return bool(refine_delta->get(bitid));
  else return 0;
}


/** Wrapper for MortonTree's block_is_parent */ 
template<unsigned ndim>
bool TheTree<ndim>::is_parent(bool updated, int bitid) {
  if (updated && in_refine){
    if (not is_updated) refine_update();
    return tree_updated->block_is_parent(bitid);
  }
  else
    return tree->block_is_parent(bitid);
}

/** Identify a block based on its integer coordinates and level of refinement. 
 *  If no block exists on the level specified, Return the block on the finest level
 *  possible. */
template<unsigned ndim>
void TheTree<ndim>::identify(
    bool updated,        //in
    int *lev,            //inout (0-based)
    int *ijk,            //inout 
    int *mort,           //out
    int *bitid           //out
    ) {
  unsigned coord[ndim];
  for(unsigned d=0; d < ndim; d++)
    coord[d] = ijk[d];

  if(updated && in_refine) {
    if (not is_updated) refine_update();
    if(tree_updated->inside(*lev, coord)) {
      typename MortonTree<ndim,W>::template Block<unsigned> b = tree_updated->identify(*lev, coord);
      *lev = b.level;
      for(unsigned d=0; d < ndim; d++)
        ijk[d] = b.coord[d];
      *mort = b.mort;
      *bitid = b.id;
    }
    else {
      *lev = -1;
      *mort = -1;
      *bitid = -1;
    }
  }
  else {
    if(tree->inside(*lev, coord)) {
      typename MortonTree<ndim,W>::template Block<unsigned> b = tree->identify(*lev, coord);
      *lev = b.level;
      for(unsigned d=0; d < ndim; d++)
        ijk[d] = b.coord[d];
      *mort = b.mort;
      *bitid = b.id;
    }
    else {
      *lev = -1;
      *mort = -1;
      *bitid = -1;
    }
  }
}

/** Return location information about a block based on its bitid number. (Location in Bittree) 
 *   */
template<unsigned ndim>
void TheTree<ndim>::locate(
    bool updated,        //in
    int bitid,           //in
    int *lev,            //out (0-based)
    int *ijk,            //out
    int *mort            //out
    ) {

  if(updated && in_refine) {
    if (not is_updated) refine_update();
    if(bitid < tree_updated->id_upper_bound() ) {
      typename MortonTree<ndim,W>::template Block<unsigned> b = tree_updated->locate(unsigned(bitid));
      *lev = b.level;
      for(unsigned d=0; d < ndim; d++)
        ijk[d] = b.coord[d];
      *mort = b.mort;
    }
    else {
      *lev = -1;
      for(unsigned d=0; d < ndim; d++)
        ijk[d] = -1;
      *mort = -1;
    }
  }
  else {
    if(bitid < tree->id_upper_bound() ) {
      typename MortonTree<ndim,W>::template Block<unsigned> b = tree->locate(unsigned(bitid));
      *lev = b.level;
      for(unsigned d=0; d < ndim; d++)
        ijk[d] = b.coord[d];
      *mort = b.mort;
    }
    else {
      *lev = -1;
      for(unsigned d=0; d < ndim; d++)
        ijk[d] = -1;
      *mort = -1;
    }
  }
}

/** Wrapper for MortonTree's bitid_list */
template<unsigned ndim>
void TheTree<ndim>::get_bitid_list(
    bool updated,       //in
    int mort_min,  //in
    int mort_max,  //in
    int *out      //out
    ) {
#ifndef BITTREE_SAFE
  if(updated && in_refine){
    if(not is_updated) refine_update();
    tree_updated->bitid_list(mort_min, mort_max, out);
  }
  else{
    tree->bitid_list(mort_min, mort_max, out);
  }
#else
  int outlist[(mort_max - mort_min)];
  if(updated && in_refine){
    if(not is_updated) refine_update();
    tree_updated->bitid_list(mort_min, mort_max, outlist);
  }
  else{
    tree->bitid_list(mort_min, mort_max, outlist);
  }
  for (int i=0;i<(mort_max - mort_min);i++)
    out[i] = outlist[i];
#endif
}


/** Creates refine_delta, and initializes all values to False.
  * First step of refinement. */
template<unsigned ndim>
void TheTree<ndim>::refine_init() {
  unsigned nbits = tree->id_upper_bound();
  refine_delta = BitArray<W>::make(nbits);
  refine_delta->fill(false);
  is_reduced = true;
  is_updated = false;
  in_refine = true;
}

/** Mark a bit on refine_delta */
template<unsigned ndim>
void TheTree<ndim>::refine_mark(
    int bitid,   // in
    bool value   // in
  ) {
  refine_delta->set(bitid, value);
  is_reduced = false;
  is_updated = false;
}


/** Reduce refine_delta across all processors by ORing. This means
 *  any blocks marked on one processor will be marked on all. */
template<unsigned ndim>
void TheTree<ndim>::refine_reduce(MPI_Comm comm) {
  MPI_Allreduce(
    MPI_IN_PLACE,
    refine_delta->word_buf(),
    refine_delta->word_count(),
    MPI_UNSIGNED,
    MPI_BOR,
    comm
  );
  is_reduced = true;
}

/** Reduce refine_delta across all processors by ANDing. This means
 *  any blocks unmarked on one processor will be unmarked on all. */
template<unsigned ndim>
void TheTree<ndim>::refine_reduce_and(MPI_Comm comm) {
  MPI_Allreduce(
    MPI_IN_PLACE,
    refine_delta->word_buf(),
    refine_delta->word_count(),
    MPI_UNSIGNED,
    MPI_BAND,
    comm
  );
  is_reduced = true;
}

/** Generates the updated tree from the original tree + refine_delta.
  * After call, both original and updated exist simultaneously. */
template<unsigned ndim>
void TheTree<ndim>::refine_update() {
  if (not is_reduced) {
    std::cout << "Bittree updating before reducing. Possible error." << endl;
  }
  tree_updated = tree->refine(refine_delta);
  is_updated = true;
}

/** Makes the updated tree the original tree. Final step of refinement. */
template<unsigned ndim>
void TheTree<ndim>::refine_apply() {
  if (not is_updated) {
    refine_update();
  }
  tree = tree_updated;
  refine_delta.nullify();
  tree_updated.nullify();
  in_refine = false;
}

/** Wrapper function to MortonTree's print_if_2d, which print a nice 
  * representation of the (2d) Bittree and refine_delta. 
  * If tree has been updated, print both original and updated version. 
  * Can be passed a datatype to change what number prints at each block loc. 
  * (0=bitid, 1=morton number, 2=parentage) */
template<unsigned ndim>
void TheTree<ndim>::print_2d(int datatype) {
  if (ndim==2) { 
    switch (datatype) {
      case 0: 
        std::cout << "printing original tree (datatype=bitid): " <<endl;
        break;
      case 1:
        std::cout << "printing original tree (datatype=mort): " <<endl;
        break;
      case 2:
        std::cout << "printing original tree (datatype=parent): " <<endl;
        break;
    }
    tree->print_if_2d(datatype);
    if(in_refine) {
    std::cout << "printing refine_delta (indexed by bitid):" <<endl;
    for (int j=0; j<refine_delta->length() ; j++){
      std::cout << j << ": " << refine_delta->get(j) << ";  " ;
    }
    }
    std::cout << endl;
    if (in_refine && is_updated) {
      switch (datatype) {
      case 0: 
        std::cout << "printing updated tree (datatype=bitid): " <<endl;
        break;
      case 1:
        std::cout << "printing updated tree (datatype=mort): " <<endl;
        break;
      case 2:
        std::cout << "printing updated tree (datatype=parent): " <<endl;
        break;
      }
      tree_updated->print_if_2d(datatype);
    }
  }
  else {
    std::cout << "Error: tried to print 2d Bittree but ndim is not 2!" <<endl;
  }
}


/** Checks if the_tree has been created */
extern "C" bool bittree_initialized() {
  return !!the_tree;
}

/** Essentially a wrapper for TheTree's constructor */
extern "C" void bittree_init(
    int *ndim,      // in
    int topsize[],  // in
    bool includes[] // in: includes[topsize[ndim-1]]...[topsize[0]]
  ) {
  unsigned top[3];
  for(unsigned d=0; d < *ndim; d++)
    top[d] = topsize[d];
  
  switch(*ndim) {
  case 1: {
    Ref<TheTree<1> > r; new(r.alloc()) TheTree<1>(top, includes);
    the_tree = r;
    } break;
  case 2: {
    Ref<TheTree<2> > r; new(r.alloc()) TheTree<2>(top, includes);
    the_tree = r;
    } break;
  case 3: {
    Ref<TheTree<3> > r; new(r.alloc()) TheTree<3>(top, includes);
    the_tree = r;
    } break;
  }
}

/** Wrapper function for block_count */
extern "C" void bittree_block_count(
    bool *updated,     //in: boolean
    int *count         //out
  ) {
  if(!!the_tree)
    *count = the_tree->block_count(*updated);
}


/** Wrapper function for leaf_count */
extern "C" void bittree_leaf_count(
    bool *updated,     //in: boolean
    int *count         //out
  ) {
  if(!!the_tree)
    *count = the_tree->leaf_count(*updated);
}

/** Wrapper function for delta_count */
extern "C" void bittree_delta_count(
    int *count         //out
  ) {
  if(!!the_tree)
    *count = the_tree->delta_count();
}


/** Wrapper function for check_refine_bit */
extern "C" void bittree_check_refine_bit(
    const int *bitid,   //in
    bool *bit_check     //out
  ) {
  if(!!the_tree)
    *bit_check = the_tree->check_refine_bit(*bitid);
}

/** Wrapper function for is_parent */
extern "C" void bittree_is_parent(
    bool *updated,      //in
    int *bitid,         //in
    bool *parent_check  //out
  ) {
  if(!!the_tree)
    *parent_check = the_tree->is_parent(*updated,*bitid);
}

/** Wrapper function for TheTree's identify, which 
  * itself wraps MortonTree's identify */
extern "C" void bittree_identify(
    bool *updated,      //in
    int *lev,           //inout (0-based)
    int *ijk,           //inout
    int *mort,          //out
    int *bitid          //out
  ) {
  if(!!the_tree)
    the_tree->identify(*updated, lev, ijk, mort, bitid);
}

/** Wrapper function for TheTree's locate, which 
  * itself wraps MortonTree's locate */
extern "C" void bittree_locate(
    bool *updated,      //in
    int *bitid,         //in
    int *lev,           //out (0-based)
    int *ijk,           //out
    int *mort          //out
  ) {
  if(!!the_tree)
    the_tree->locate(*updated, *bitid, lev, ijk, mort);
}

/** Wrapper function for TheTree's get_bitid_list, which 
  * itself wraps MortonTree's bitid_list */
extern "C" void bittree_get_bitid_list(
    bool *updated,      //in
    int *mort_min,      //in, 0-based
    int *mort_max,      //in, 0-based
    int *idout          //out
  ) {
  if(!!the_tree)
    the_tree->get_bitid_list(*updated, *mort_min, *mort_max, idout);
}

/** Wrapper function for refine_init */
extern "C" void bittree_refine_init() {
  if(!!the_tree)
    the_tree->refine_init();
}

/** Wrapper funciton for refine_mark */
extern "C" void bittree_refine_mark(
    int *bitid,        // in
    bool *value        // in
  ) {
  if(!!the_tree)
    the_tree->refine_mark(*bitid, *value);
}

/** Wrapper function for refine_reduce */
extern "C" void bittree_refine_reduce(int *comm_) {
  if(!!the_tree) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_);
    the_tree->refine_reduce(comm);
  }
}

/** Wrapper function for refine_reduce_and */
extern "C" void bittree_refine_reduce_and(int *comm_) {
  if(!!the_tree) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_);
    the_tree->refine_reduce_and(comm);
  }
}

/** Wrapper function for refine_update */
extern "C" void bittree_refine_update() {
  if(!!the_tree)
    the_tree->refine_update();
}

/** Wrapper function for refine_apply */
extern "C" void bittree_refine_apply() {
  if(!!the_tree)
    the_tree->refine_apply();
}

/** Wrapper function the print_2d */
extern "C" void bittree_print_2d(int *datatype=0)
{
  if(!!the_tree)
    the_tree->print_2d(*datatype) ;
}

#if 0
int main() {
  unsigned size[2] = {3,3};
  unsigned excl[1][2] = {1,1};
  Ref<MortonTree<2,unsigned> > tree = MortonTree<2,unsigned>::make(size, 1, excl);
  //for(int i=0; i < 10; i++) {
  for(int i=0; i < 3; i++) {
    Ref<BitArray<unsigned> > delta = BitArray<unsigned>::make(tree->level_id0(0) + tree->blocks());
    delta->fill(false, 0, tree->level_id0(tree->levels()-1));
    delta->fill(true, tree->level_id0(tree->levels()-1), tree->id_upper_bound());
    //delta->set(tree->level_id0(tree->levels()-1), true);
    //delta->set(tree->blocks()-1, true);
    tree = tree->refine(delta);
    cout << "i=" << i << " tree blocks=" << tree->blocks() << '\n';
  }
  if(false) {
    unsigned sum = 0;
    for(int i=0; i < 1<<14; i++) {
      unsigned x[2] = { i%(1<<11), i%(1<<11) };
      MortonTree<2,unsigned>::Block<unsigned> b0 = tree->identify(10, x);
      /*MortonTree<2,unsigned>::Block<unsigned> b1 = tree->locate<unsigned>(b0.id);
      DBG_ASSERT(b0.id == b1.id);
      DBG_ASSERT(b0.level == b1.level);
      DBG_ASSERT(b0.mort_ix == b1.mort_ix);
      DBG_ASSERT(b0.coord[0] == b1.coord[0] && b0.coord[1] == b1.coord[1]);
      */sum = 5*sum ^ b0.mort;
    }
    cout << "sum:" << sum << '\n';
  }
  else
    cout << (MortonTree<2,unsigned>*)tree;
}
#endif
