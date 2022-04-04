#	Makefile for unit Hydro/HydroMain/Spark
Hydro += Hydro_data.o Hydro_inhost_data.o Hydro_offload_data.o Hydro_init.o Hydro_computeDt.o \
				  hy_rk_getFaceFlux.o hy_rk_interface.o \
			   hy_rk_updateSoln.o Hydro_correctGravSrc.o \
				 hy_rk_getGravAccel.o hy_rk_correctFluxes.o 

hy_rk_getFaceFlux.o :  hydro_data.mod Hydro_data.o timers_interface.mod driver_interface.mod

