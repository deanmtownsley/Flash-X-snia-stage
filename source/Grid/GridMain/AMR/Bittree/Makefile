bittree_core.o: bittree_core.cxx bittree_bitarray.hxx \
 bittree_bitarray_defs.hxx bittree_prelude.hxx bittree_ref_defs.hxx \
 bittree_mem_defs.hxx bittree_bits.hxx bittree_bits_defs.hxx \
 bittree_mem.hxx bittree_ref.hxx bittree_mortontree.hxx \
 bittree_mortontree_defs.hxx Simulation.h

Grid += \
  bittree_core.o \
  bittree.o \
  gr_btIdentify.o \
  gr_btLocate.o \
  gr_btGetBitid.o \
  gr_btIsParent.o \
  gr_btRefineMark.o \
  gr_btDerefineMark.o \
  gr_btGetRefine.o \
  gr_btGetDerefine.o \
  gr_btGetLocalBitids.o


# The following GNU make special prevents that apparent dependencies
# on the file iso_c_binding.mod, which does usually not actually exist
# in the object directory but refers to the ISO_C_BINDING module known
# to FORTRAN compilers internally, trigger unnecessary recompilation
# of files that refer to the ISO_C_BINDING module.
.SECONDARY: iso_c_binding.mod
