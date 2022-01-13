
#!
#! Modification history:
#!     Michael L. Rilee, November 2002, *dbz*
#!        Initial support for divergenceless prolongation
#!     Michael L. Rilee, December 2002, *clean_divb*
#!        Support for projecting field onto divergenceless field
#!

.SUFFIXES :  
.SUFFIXES : .o .c .f .F .F90 .f90 .fh .a

sources := \
  amr_1blk_cc_cp_remote.F90  \
  amr_1blk_cc_prol_gen_unk_fun.F90  \
  amr_1blk_cc_prol_inject.F90  \
  amr_1blk_cc_prol_linear.F90  \
  amr_1blk_cc_prol_genorder.F90 \
  amr_1blk_cc_prol_user.F90 \
  amr_1blk_cc_prol_gen_work_fun.F90 \
  amr_1blk_cc_prol_work_inject.F90 \
  amr_1blk_cc_prol_work_linear.F90 \
  amr_1blk_cc_prol_work_genorder.F90 \
  amr_1blk_cc_prol_work_user.F90 \
  amr_1blk_copy_soln.F90    \
  amr_1blk_ec_cp_remote.F90 \
  amr_1blk_ec_prol_gen_fun.F90 \
  amr_1blk_ec_prol_genorder.F90 \
  amr_1blk_ec_prol_user.F90 \
  amr_1blk_ec_prol_linear.F90 \
  amr_1blk_fc_prol_dbz.F90 \
  clean_field.F90 \
  poisson_sor.F90 \
  amr_1blk_fc_clean_divb.F90 \
  amr_1blk_fc_cp_remote.F90 \
  amr_1blk_fc_prol_gen_fun.F90 \
  amr_1blk_fc_prol_inject.F90 \
  amr_1blk_fc_prol_linear.F90 \
  amr_1blk_fc_prol_genorder.F90 \
  amr_1blk_fc_prol_user.F90 \
  amr_1blk_guardcell_reset.F90 \
  amr_1blk_guardcell_srl.F90 \
  set_f2c_indexes.F90 \
  amr_1blk_nc_cp_remote.F90 \
  amr_1blk_nc_prol_gen_fun.F90 \
  amr_1blk_nc_prol_genorder.F90 \
  amr_1blk_nc_prol_user.F90 \
  amr_1blk_nc_prol_linear.F90 \
  amr_1blk_save_soln.F90 \
  amr_1blk_t_to_perm.F90 \
  amr_1blk_to_perm.F90 \
  amr_bcset_init.F90 \
  amr_block_geometry.F90 \
  user_coord_transfm.F90 \
  amr_close.F90 \
  amr_initialize.F90 \
  amr_set_runtime_parameters.F90 \
  amr_mpi_find_blk_in_buffer.F90 \
  amr_perm_to_1blk.F90 \
  amr_prolong_cc_fun_init.F90 \
  amr_prolong_face_fun_init.F90 \
  amr_prolong_fun_init.F90 \
  amr_reorder_grid.F90 \
  amr_restrict_ec_fun.F90 \
  amr_restrict_ec_genorder.F90 \
  amr_restrict_ec_user.F90 \
  amr_restrict_edge.F90 \
  amr_restrict_fc_fun.F90 \
  amr_restrict_fc_genorder.F90 \
  amr_restrict_fc_user.F90 \
  amr_restrict_red.F90 \
  amr_restrict_unk_fun.F90 \
  amr_restrict_unk_genorder.F90 \
  amr_restrict_unk_user.F90 \
  amr_restrict_nc_fun.F90 \
  amr_restrict_nc_user.F90 \
  amr_restrict_nc_genorder.F90 \
  amr_restrict_work_fun.F90 \
  amr_restrict_work_genorder.F90 \
  amr_restrict_work_user.F90 \
  amr_restrict_work_fun_recip.F90 \
  amr_system_calls.F90 \
  amr_q_sort.F90 \
  amr_q_sort_real.F90 \
  mpi_amr_singular_line.F90 \
  mpi_amr_1blk_guardcell.F90 \
  mpi_amr_1blk_guardcell_c_to_f.F90 \
  mpi_amr_1blk_restrict.F90 \
  mpi_amr_comm_setup.F90 \
  mpi_amr_edge_average.F90 \
  mpi_amr_edge_average_udt.F90 \
  mpi_amr_edge_average_vdt.F90 \
  mpi_amr_edge_diagonal_check.F90 \
  mpi_amr_flux_conserve.F90 \
  mpi_amr_flux_conserve_udt.F90 \
  mpi_amr_flux_conserve_vdt.F90 \
  mpi_amr_get_remote_block.F90 \
  mpi_amr_get_remote_block_fvar.F90 \
  mpi_amr_global_domain_limits.F90 \
  mpi_amr_guardcell.F90 \
  mpi_amr_local_surr_blks_lkup.F90 \
  mpi_amr_prolong.F90 \
  mpi_amr_prolong_fc_divbconsist.F90 \
  mpi_amr_refine_derefine.F90 \
  amr_morton_process.F90 \
  local_tree.F90 \
  find_surrblks.F90 \
  local_tree_build.F90 \
  tree_search_for_surrblks.F90 \
  mpi_amr_restrict.F90 \
  mpi_amr_restrict_fulltree.F90 \
  mpi_amr_restrict_bnd_data_vdt.F90 \
  mpi_amr_restrict_edge_data_vdt.F90 \
  mpi_amr_store_comm_info.F90 \
  mpi_amr_timing_report.F90 \
  mpi_amr_tree_setup.F90 \
  mpi_get_buffer.F90 \
  mpi_get_edge_buffer.F90 \
  mpi_get_flux_buffer.F90 \
  mpi_morton_bnd.F90 \
  process_fetch_list.F90 \
  compress_fetch_list.F90 \
  mpi_morton_bnd_fluxcon.F90 \
  mpi_morton_bnd_prolong.F90 \
  mpi_morton_bnd_restrict.F90 \
  mpi_pack_blocks.F90 \
  mpi_unpack_blocks.F90 \
  mpi_put_buffer.F90 \
  mpi_pack_edges.F90 \
  mpi_pack_fluxes.F90 \
  mpi_put_edge_buffer.F90 \
  mpi_put_edge_buffer_1blk.F90 \
  mpi_put_flux_buffer.F90 \
  mpi_set_message_limits.F90 \
  mpi_set_message_sizes.F90 \
  mpi_unpack_edges.F90 \
  mpi_unpack_fluxes.F90 \
  rationalize_fetch_list.F90 \
  mpi_amr_checkpoint_wr.F90 \
  mpi_amr_checkpoint_re.F90 \
  mpi_amr_checkpoint_wr_default.F90 \
  mpi_amr_checkpoint_re_default.F90 \
  mpi_amr_checkpoint_wr_hdf5.F90 \
  mpi_amr_checkpoint_re_hdf5.F90 \
  mpi_amr_checkpoint_wr_mpiio.F90 \
  mpi_amr_checkpoint_re_mpiio.F90 \
  read_blocks_hdf5_r4.c \
  read_blocks_hdf5_r8.c \
  write_blocks_hdf5_r4.c \
  write_blocks_hdf5_r8.c \
  mpi_amr_plotfile_chombo.F90 \
  write_blocks_chombo_r4.c \
  write_blocks_chombo_r8.c \
  mpi_amr_derefine_blocks.F90 \
  mpi_amr_check_derefine.F90 \
  amr_morton_order.F90 \
  amr_compute_morton.F90 \
  amr_sort_morton.F90 \
  amr_sort_morton_reorder_grid.F90 \
  amr_sort_by_work.F90 \
  amr_migrate_tree_data.F90 \
  fill_old_loc.F90 \
  morton_sort.F90 \
  mpi_amr_redist_blk.F90 \
  send_block_data.F90 \
  mpi_amr_refine_blocks.F90 \
  amr_check_refine.F90 \
  mpi_amr_restrict_bnd_data.F90 \
  mpi_amr_restrict_edge_data.F90 \
  mpi_amr_boundary_block_info.F90 \
  mpi_amr_test_neigh_values.F90 \
  mpi_lib.F90 \
  mpi_pack_tree_info.F90 \
  mpi_unpack_tree_info.F90 \
  mpi_wrapper_int.F90 \
  mpi_wrapper_logical.F90 \
  mpi_wrapper_real.F90 \
  mpi_wrapper_dble.F90
 
%.o:%.F90
	$(FC) -c $(FFLAGS) $<

objects := \
	$(patsubst %.F90,%.o, \
	$(patsubst %.f90,%.o, \
	$(patsubst %.c,%.o, \
	$(patsubst %.F,%.o,$(sources)))))

vpath %.fh ../headers

# ar does not work on Macs for some reason (???)
#libparamesh.a: $(objects)
#	$(AR) $(ARFLAGS) $@ $^
#	cp -f libparamesh.a ../libs/.
libparamesh.a: $(objects)
	cp *.o ../libs/.

ifdef MY_CPP
GNUmakefile.include: $(sources)
	find . -name \*[.f90,.F90,.c,.F] | xargs $(MY_CPP) > $@
include GNUmakefile.include
else
$(objects): $(wildcard *.fh)
endif

.PHONY: clean
clean:
	$(RM) libparamesh.a *.o *~ GNUmakefile.include libmpi_paramesh.a *.i *~

#--------------------------------------------------------------------------

.SUFFIXES : .F90 .f90 .c .o

.F.o:
	$(FC) $(FFLAGS) -c -o $*.o $*.F

.f90.o:
	$(FC) $(FFLAGS) -c -o $*.o $*.f90

.F90.o:
	$(FC) $(FFLAGS) -c -o $*.o $*.F90

.c.o:
	$(CC) $(CFLAGS) -c -o $*.o $*.c

# .fh:;
# .mod:;

#--------------------------------------------------------------------------
