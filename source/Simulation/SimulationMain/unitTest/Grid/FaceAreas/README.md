#This test has been tested with gcc/10.2.0 compiler in debug or with -O2 optimization

#Below provided setup lines are for amr test. To perform a UG test, replace '+pm4dev' with '+ug -nofbs', and 'test_*.par with 'UG_test_*.par'

#Setup line for Cartesian test runs
#1D test
./setup unitTest/Grid/FaceAreas -auto -1d +cartesian +pm4dev -nxb=128 +noio -parfile=test_car.par
#2D test
./setup unitTest/Grid/FaceAreas -auto -2d +cartesian +pm4dev -nxb=128 -nyb=128 +noio -parfile=test_car.par 
#3D test
./setup unitTest/Grid/FaceAreas -auto -3d +cartesian +pm4dev -nxb=128 -nyb=128 -nzb=128 +noio -parfile=test_car.par

#Setup line for Cylindrical test runs
#1D test
./setup unitTest/Grid/FaceAreas -auto -1d +cylindrical +pm4dev -nxb=128 +noio -parfile=test_cyl.par
#2D test
./setup unitTest/Grid/FaceAreas -auto -2d +cylindrical +pm4dev -nxb=128 -nyb=128 +noio -parfile=test_cyl.par 
#3D test
./setup unitTest/Grid/FaceAreas -auto -3d +cylindrical +pm4dev -nxb=128 -nyb=128 -nzb=32 +noio -parfile=test_cyl.par

#Setup line for Spherical test runs
#1D test
./setup unitTest/Grid/FaceAreas -auto -1d +spherical +pm4dev -nxb=128 +noio -parfile=test_sph.par
#2D test
./setup unitTest/Grid/FaceAreas -auto -2d +spherical +pm4dev -nxb=128 -nyb=32 +noio -parfile=test_sph.par 
#3D test
./setup unitTest/Grid/FaceAreas -auto -3d +spherical +pm4dev -nxb=128 -nyb=32 -nzb=32 +noio -parfile=test_sph.par


