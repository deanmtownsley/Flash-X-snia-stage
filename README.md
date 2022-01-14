# Flash-X


## Git/Testing Workflow

The current rules for collaborating via git are as follows

1.  Base all feature branches off of the main branch.
2.  When development on a feature branch is finished, merge main branch on to the feature branch and run full test-suite on the feature branch. 
3.  If at any point during the previous steps errors are discovered that need to
    be fixed, fix code in the feature branch, then go back to point 2 above. Repeat until no test errors.
4. If a merge conflict occurs when merging main into the feature branch _do not_ attempt to resolve conflicts using the  GitHub web interface - such an attempt can results in an unintended merge.	
5.  Do not rebase a feature branch that has already been pushed to the GitHub
    repository.


## Libraries

The code of the PARMESH library is currently integrated into the GridMain implementation as a
subtree and not organized or built as a separate library.

Some applications and tests use external libraries that are expected to be already installed on the
system where Flash-X is being built and run. The directory locations of such library installations
should be made know to the Flash-X build system by a site-specific (or, as a fallback, OS-specific)
Makefile.h file. See the subdirectories under sites/ .

This applies in particular to the AMReX library. Separate library instances for 1D, 2D, and 3D
should be installed, and the appropriate locations mentioned in Makefile.h .

On the other hand, some applications and tests use INTERNAL libraries. Such libraries are built, as part of the Flash-X setup process, from source code that is located in subdirectories under lib/ . There are two cases for such libraries:

1. **Library source code is included as part of the Flash-X git respository.**

   An example is the sqrt3 library, whose source code is included in lib/sqrt/ .

 **Library source code must be retrieved from a separate repository.**

   Examples are the THORNADO and WEAKLIB libraries.
   Follow the instructions on submodules to automatically put the source code for these two
   in the right places in subdirectories under lib/.
   
## Git with Submodules

To prepare for building simulations that use libraries whose code must be retrieved
from a separate git repository, the following modified `git` commands can be used:

- `git pull --recurse-submodules=yes` (in place of the usual `git pull`)
- `git submodule update --init` (additionally, after `git pull`)

## Using MAPLE Container

MAPLE is a Python API and CLI that acts as a wrapper around ```docker```/```singularity``` to implement containerization of HPC applications and their developer environment. With MAPLE, one can run, compile, and develop ```Flash-X``` inside a container without having to install external libraries on a new machine. The only dependency for MAPLE is ```docker``` and ```singularity```.

We intend to develop MAPLE as a standalone application for multiple projects, and therefore provide it as a submodule within ```Flash-X```

If you want to test this new feature please follow these steps:

### Initialize ```maple``` submodule

  ```
  mkdir -p $HOME/.local/bin
  export PATH="$PATH:$HOME/.local/bin"

  git submodule update --init $Flash-X_HOME/tools/maple
  cd $Flash-X_HOME/tools/maple
  ./setup develop
  ./setup clean
  ```

### Writing a Maplefile

  ```Maplefile``` is used to define environment variables required by ```maple```. Following is a list of variables:
  
  	```maple_image```: Name of the image in remote registry 
  	
  	```maple_container```: Name of the local container
  	
	```maple_target```: Name of the target dir to mount src dir
	
	```maple_port```: Port ID for jupyter notebooks

 	```maple_backend```: Backend (docker/singularity)
       	
  ```maple``` passes these variables to its internal ```Dockerfile``` to build the images and containers.  

  Please refer to example ```Maplefile``` in ```Flash-X``` directory

### CLI use:

uitilty  - Build maple image from local image using ```maple build```

  - Getting shell access:

    ```maple pour```: to pour image into a container to enable shell access (only available with docker backend)

    ```maple shell```: provides shell access to the container

    ```maple commit```: Save changes from local container to local image (only available with docker backend)

    ```maple squash```: Prune redundant layers from a local container and save it to local image (do this to reduce size of an image, only available with docker backend)

    ```maple rinse```: this commands stops and deletes the local container (only available with docker backend)

  - Launch an ipython notebook inside the 

    ```maple notebook```: launches the notebook server

  - Execute commands inside the container

    ```maple execute "echo Hello World!"```: example to launch specific command inside the container

  - Cleanup

    ```maple clean```: deletes the local image, if you want to update remote image with changes to local image run ```maple push <remote_image_name:tag>``` before ```maple clean```

    ```maple remove```: deletes the instance of remote image on local machine, doing this means that ```maple build``` will have to perform the expensive task of pulling the remote image again if you decide to rebuild the local image.

  - To compile ```Flash-X``` inside the container use ```/home/site``` as the site directory

    ```
    ./setup incompFlow/PoolBoiling -auto -2d -site=/home/site +amrex -maxblocks=100
    ```

### API use:

  Python API provides a way to package ```Flash-X``` simulations for deployment within cloud computing environments and machine learning workflows. Example scripts are located in ```container``` sub-directory, which is designed to match ```Simulation/SimulationMain``` directory structure.

  Run ```python3 container/incompFlow/PoolBoiling.py``` and ```python3 container/incompFlow/RisingBubble.py``` to see how it works.
uitilty

## Tests 

The following are the setup commands of the tests that are currently included in the test-suite to confirm correctness of the functionality we consider functional in Flash-X.

#### Unit Tests

- unitTest/Eos/Helmholtz -auto +amrex -3d +noio
- unitTest/Eos/Helmholtz -auto +pm4dev -3d +noio
- unitTest/Gravity/Poisson3 -auto -2d +cylindrical +newmpole -debug -maxblocks=600 +noio +pm4dev -parfile=flash_2dcyl.par
- unitTest/Gravity/Poisson3 -auto -3d +newmpole +uhd -debug -maxblocks=550 -nxb=8 -nyb=8 -nzb=8 -gridinterpolation=monotonic
- unitTest/Grid/Amrex/TestFluxCorrection -auto -2d -nxb=8 -nyb=8 +noio +amrex
- unitTest/Grid/Amrex/TestFluxCorrection2 -auto -2d -nxb=8 -nyb=8 +noio +amrex
- unitTest/Grid/Amrex/TestInit -auto -2d -nxb=8 -nyb=4 +noio +amrex
- unitTest/Grid/Amrex/TestRefine -auto -2d -nxb=8 -nyb=8 +noio +amrex
- unitTest/Grid/Amrex/TestCyl2d -auto -2d -nxb=8 -nyb=4 +noio +amrex
- unitTest/Multigrid_Amrex -auto -3d +amrex +serialio -unit=IO/IOMain/hdf5/serial/AM -maxblocks=1000
- unitTest/Gravity/PointMass -auto -2d +pm4dev +cylindrical +serialIO -debug -parfile=test_2dcyl.par
- unitTest/Gravity/PointMass -auto -2d +amrex +cylindrical +noio -debug -parfile=test_2dcyl.par
- unitTest/Gravity/PointMass -auto -3d +pm4dev -maxblocks=550 -nxb=8 -nyb=8 -nzb=8
- unitTest/Gravity/PointMass -auto -3d +amrex +noio -maxblocks=10 -nxb=8 -nyb=8 -nzb=8
- unitTest/Grid/AnomalousRefine -auto -2d +spherical -nxb=8 -nyb=8
- unitTest/Grid/AnomalousRefine -auto -2d +spherical +amrex -nxb=8 -nyb=8 -unit=IO/IOMain/hdf5/serial/AM
- unitTest/IO -auto --index-reorder -3d +cube16 +parallelIO nVars=25 +hdf5AsyncIO

#### Regression Tests with Verified Baseline

- Sod -auto -2d -test +sHLL +ug +nofbs -parfile=test_ug_TBL_2d.par
- Sod -auto -2d -debug +uhd +ug +nofbs -parfile=test_pseudoug_2d.par
- Sod -auto -2d -debug +uhd +pm4dev +nolwf -gridinterpolation=monotonic -parfile=test\_{amr_unsplit,pseudoug}\_2d.par
- Sod -auto -2d -debug -nxb=8 -nyb=8 +sHLL +pm4dev -gridinterpolation=native -parfile=test\_{amr,pseudoug}\_2d.par
- Sod -auto -2d -nxb=8 -nyb=8 +sHLL +pm4dev Bittree=True -gridinterpolation=native -parfile=test_amr_2d.par
- Sod -auto -2d -debug +uhd +amrex +nolwf +serialio -unit=IO/IOMain/hdf5/serial/AM -parfile=test\_{amr_unsplit,pseudoug}\_2d.par
- Sod -auto -2d -debug +sHLL +amrex +nolwf +serialio -unit=IO/IOMain/hdf5/serial/AM -parfile=test\_{amr,pseudoug}\_2d.par
- Sod -auto -2d -debug +spark +ug +nofbs -parfile=test_pseudoug_2d.par
- Sod -auto -2d -nxb=12 -nyb=12 -debug +spark +pm4dev -gridinterpolation=monotonic -parfile=test_pseudoug_2d.par
- Sod -auto -2d -nxb=16 -nyb=16 -debug +spark +amrex +serialio -unit=IO/IOMain/hdf5/serial/AM -parfile=test_pseudoug.par
- Sod -auto -3d +cube16 +uhd +pm4dev +nolwf -parfile=test3d-1node_4lev.par
- Sedov -auto -2d -debug +uhd +ug +nofbs -parfile= test_pseudoug_2d.par
- Sedov -auto -3d -debug -nxb=8 -nyb=8 -nzb=8 +uhd +pm4dev -gridinterpolation=native -parfile=test\_{amr_unsplit,pseudoug}\_3d.par
- Sedov -auto -3d -nxb=8 -nyb=8 -nzb=8 +uhd +pm4dev Bittree=True -gridinterpolation=native -parfile=test_amr_unsplit_3d.par
- Sedov -auto -3d -debug +uhd +amrex +nolwf +serialio -unit=IO/IOMain/hdf5/serial/AM -parfile=test\_{amr_unsplit,pseudoug}\_3d.par
- Sedov -auto -2d +cylindrical -debug -nxb=16 -nyb=16 +uhd +pm4dev -gridinterpolation=monotonic DoAnalytical=True -parfile=test\_{amr,pseudoug}\_cyl_2d.par
- Sedov -auto -2d +cylindrical -debug -nxb=16 -nyb=16 +uhd +amrex +serialio -unit=IO/IOMain/hdf5/serial/AM DoAnalytical=True -parfile=test\_{amr,pseudoug}\_cyl_2d.par
- Sedov -auto -2d -debug -nxb=8 -nyb=8 +uhd +pm4dev -gridinterpolation=native -parfile=test\_{amr_unsplit,pseudoug}\_2d.par
- Sedov -auto -2d -debug -nxb=8 -nyb=8 +sHLL +pm4dev -gridinterpolation=native -parfile=test\_{amr,pseudoug}\_2d.par
- Sedov -auto -2d -debug +uhd +amrex +nolwf +serialio -unit=IO/IOMain/hdf5/serial/AM -parfile=test\_{amr_unsplit,pseudoug}\_2d.par
- Sedov -auto -2d -debug +sHLL +amrex +nolwf +serialio -unit=IO/IOMain/hdf5/serial/AM -parfile=test\_{amr,pseudoug}\_2d.par
- Sedov -auto -2d -nxb=12 -nyb=12 -debug +spark +ug +nofbs -parfile=test_pseudoug_2d.par
- Sedov -auto -3d -nxb=16 -nyb=16 -nzb=16 -debug +spark +amrex +serialio -unit=IO/IOMain/hdf5/serial/AM -parfile=test_pseudoug_3d.par
- Sedov -auto -2d +cylindrical -debug -nxb=16 -nyb=16 +spark +pm4dev -gridinterpolation=monotonic DoAnalytical=True -parfile=test_amr_cyl_2d.par
- Sedov -auto -2d -debug +spark +amrex +serialio -unit=IO/IOMain/hdf5/serial/AM -parfile=test_amr_unsplit_2d.par
- StreamingSineWave -auto -3d +cartesian -nxb=8 -nyb=8 -nzb=8 nE=2 nSpecies=1 nNodes=2 nMoments=4 momentClosure=MINERBO -parfile=test_paramesh_3d{,_restart}.par
- StreamingSineWave -auto -3d +cartesian nE=2 nSpecies=1 nNodes=2 nMoments=4 momentClosure=MINERBO +Mode3 +nolwf -parfile=test_amrex_3d.par
- Cellular -auto -2d -debug +a13 +uhd +pm4dev -gridinterpolation=monotonic -parfile=test_amr_2d.par
- SNIa_DoubleDetonation -auto -2d -test +uhd +nolwf +pm4dev +cylindrical -nxb=16 -nyb=16 +newMpole +xnet xnetData=Data_alpha Bittree=True AltMorton=True -parfile=testfaster_amr_unsplit_2d.par
- SNIa_DoubleDetonation -auto -2d -debug +cylindrical +Mode1 -nxb=16 -nyb=16 +newMpole +a13 -parfile=test_shellDet_2d.par
- SNIa_DoubleDetonation -auto -2d -test +spark +pm4dev +cylindrical -nxb=16 -nyb=16 +newMpole +xnet xnetData=Data_alpha Bittree=True AltMorton=True -parfile=testfaster_amr_unsplit_2d.par
- CCSN_WL -auto -1d +spherical -nxb=16 threadBlockList=False +pm4dev +hdf5 threadWithinBlock=False +newMpole +uhd -debug +mode1 -parfile=ccsn1d.par
- YahilLattimerCollapse -debug -auto -1d +spherical -maxblocks=16000 +hdf5 -without-unit=Grid/GridSolvers/Multipole -unit=Grid/GridSolvers/Multipole_new +uhd +nolwf -parfile=yahil_1d.par
- YahilLattimerCollapse -debug -auto -2d +cylindrical -maxblocks=16000 +hdf5 -without-unit=Grid/GridSolvers/Multipole -unit=Grid/GridSolvers/Multipole_new +uhd +nolwf -parfile=yahil_2d.par
- YahilLattimerCollapse -debug -auto -1d -nxb=12 +spherical -maxblocks=16000 +hdf5 -without-unit=Grid/GridSolvers/Multipole -unit=Grid/GridSolvers/Multipole_new +spark -parfile=yahil_1d.par
- YahilLattimerCollapse -debug -auto -2d -nxb=12 -nyb=12 +cylindrical -maxblocks=16000 +hdf5 -without-unit=Grid/GridSolvers/Multipole -unit=Grid/GridSolvers/Multipole_new +spark -parfile=yahil_2d.par
- HydroStatic -auto -2d -test +Mode1 +nolwf useFortran2003=True -parfile=flash.par

#### Regression Tests with Unverified Baseline

- Sod -auto -2d -debug -nxb=8 -nyb=8 +uhd +pm4dev AltMorton=True Bittree=True +threadBL -gridinterpolation=monotonic -parfile=test_amr_TBL_unsplit_2d.par
- DustCollapse -auto -3d +cartesian +Mode1 +serialIO +uhd +newMpole -debug -parfile=test_3dcar.par
- DustCollapse -auto -3d +cartesian +amrex +parallelIO +uhd +newMpole -debug -parfile=test_3dcar.debug.par
- DustCollapse -auto -2d +cylindrical +Mode3 +serialIO +uhd +newMpole -debug -parfile=test_2dcyl.debug.par
- DustCollapse -auto -2d +cylindrical +Mode1 AltMorton=True +serialIO +uhd +newMpole -parfile=test_2dcyl.par
- DustCollapse -auto -2d +cylindrical +Mode1 AltMorton=True Bittree=True +serialIO +uhd +newMpole -parfile=test_2dcyl.par
- DustCollapse -auto -1d +spherical +Mode1 +serialIO +uhd +newMpole -debug -parfile=test_1dsph.par
- DustCollapse -auto -1d +spherical +Mode3 +serialIO +uhd +newMpole -debug +nolwf -parfile=test_1dsph.debug.par
- DustCollapse -auto -1d -nxb=16 +spherical -unit=Grid/GridMain/AMR/Amrex +serialio -unit=IO/IOMain/hdf5/serial/AM --index-reorder +serialIO +spark +newMpole -debug -parfile=test_1dsph.debug.par
- IsentropicVortex -auto -2d -debug +uhd +amrex +nolwf +serialIO -unit=IO/IOMain/hdf5/serial/AM -unit=Particles
- IsentropicVortex -auto -2d -nxb=16 -nyb=16 -debug +spark +amrex +serialIO -unit=IO/IOMain/hdf5/serial/AM -unit=Particles
