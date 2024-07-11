
include $(MESA_DIR)/utils/makefile_header


MESA_LIBS_NEEDED = -L$(MESA_LIB_DIR) -leos -lrates -lchem -lnet $(LOAD_MESA_NUMERICS)

TEST_OBJS = mesa_flash_nucreact.o test_burn_const_dens_temp.o

test_burn_const_dens_temp: $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o test_burn_const_dens_temp $(TEST_OBJS) $(MESA_LIBS_NEEDED)

mesa_flash_nucreact.o: ../../mesa_flash_nucreact.F90
	$(TEST_COMPILE) $(FCfree) $<

test_burn_const_dens_temp.o: test_burn_const_dens_temp.F90
	$(TEST_COMPILE) $(FCfree) $<

clean:
	-@rm -f *.o *.mod $(TEST_OBJS) test_burn_const_dens_temp
