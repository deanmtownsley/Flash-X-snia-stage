
include $(MESA_DIR)/utils/makefile_header


MESA_LIBS_NEEDED = -L$(MESA_LIB_DIR) -leos -lrates -lchem -lnet $(LOAD_MESA_NUMERICS)

GEN_OBJS = mesa_flash_nucreact.o generate_flash_species_config.o

generate_flash_species_config: $(GEN_OBJS)
	$(LOADER) $(FCopenmp) -o generate_flash_species_config $(GEN_OBJS) $(MESA_LIBS_NEEDED)

mesa_flash_nucreact.o: ../mesa_flash_nucreact.F90
	$(TEST_COMPILE) $(FCfree) $<

generate_flash_species_config.o: generate_flash_species_config.F90
	$(TEST_COMPILE) $(FCfree) $<

clean:
	-@rm -f *.o *.mod $(GEN_OBJS) generate_flash_species_config
