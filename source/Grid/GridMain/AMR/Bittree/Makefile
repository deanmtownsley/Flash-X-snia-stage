bittree_core.o: bittree_core.cxx bittree_bitarray.hxx \
 bittree_bitarray_defs.hxx bittree_prelude.hxx bittree_ref_defs.hxx \
 bittree_mem_defs.hxx bittree_bits.hxx bittree_bits_defs.hxx \
 bittree_mem.hxx bittree_ref.hxx bittree_mortontree.hxx \
 bittree_mortontree_defs.hxx Simulation.h

Grid += \
  bittree_core.o \
  bittree.o \
  amr_identify_block.o \
  amr_bittree_locate.o \
  amr_bittree_get_bitid.o \
  amr_bittree_is_parent.o \
  amr_bittree_refine_mark.o \
  amr_bittree_derefine_mark.o \
  amr_bittree_get_refine.o \
  amr_bittree_get_derefine.o \
  amr_sort_morton_bittree.o \
  amr_bittree_get_local_bitids.o


# The following GNU make special prevents that apparent dependencies
# on the file iso_c_binding.mod, which does usually not actually exist
# in the object directory but refers to the ISO_C_BINDING module known
# to FORTRAN compilers internally, trigger unnecessary recompilation
# of files that refer to the ISO_C_BINDING module.
.SECONDARY: iso_c_binding.mod
