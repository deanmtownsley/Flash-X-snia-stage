Grid += gr_mpoleInterface.o gr_mpoleMoments.o gr_mpoleMomBins2Dspherical.o gr_mpoleMom2Dspherical.o gr_mpolePot2Dspherical.o gr_mpoleRad2Dspherical.o gr_mpoleCen2Dspherical.o #faceareaStore.o check_variable_sneo.o
Simulation += Simulation_data.o chimera_model.o model_interp.o quadrature.o hdf5_read_write.o

#Gravity_accelOneRow.o : Grid_interface.mod
#check_variable_sneo.o : faceareaStore.o
Simulation_init.o : Simulation_data.o chimera_model.o quadrature.o
Simulation_initBlock.o : Simulation_data.o chimera_model.o quadrature.o
chimera_model.o : Simulation_data.o model_interp.o quadrature.o hdf5_read_write.o
mesa_model.o : Simulation_data.o model_interp.o quadrature.o hdf5_read_write.o
