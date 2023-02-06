## Milhoja Initialization Unit Test Design

This document is intended to be a living document that captures and communicates the __current__ design and use of the Milhoja `TestInit` unit test.  In this sense, users can assume that this document and all the contents of this folder are correct, up-to-date, and an encapsulation of the official Flash-X tests related to Milhoja initialization.

#### Keywords
Milhoja, Grid, Interface, Pseudo-UG, Cartesian, Unit Test, 1D, 2D, 3D

#### Motivation
This unit test was developed via Test-Driven Development (TDD) of the Milhoja grid backend's Fortran/C++ interoperability layer and the Flash-X Milhoja Grid unit implementation.  In particular, it was useful for building up functionality before the possibility of saving data to Flash-X format data files.  Since it is a simple and quick unit test and its checks are generally useful, there is no immediate need to remove the test.  Indeed, it might be invaluable for quick and correct TDD-based development of subsequent Milhoja grid backends.  

This test was designed largely to confirm basic Grid unit functionality (e.g., setting initial conditions and filling GCs) as well as to confirm that significant parts of the public and private Grid interfaces are correctly implemented by Milhoja's grid backend.  In addition, it also confirms that some private Grid variables were correctly initialized.  Therefore, including [123]D versions of this test in the testsuite seems good and useful.

One important detail of this test is that it confirms that data in arrays indexed by dimension have reasonable values above NDIM.  To the best of my knowledge there is no specification for what this data should be for many such arrays.  Therefore including [123]D versions of this test here will at least ensure that all Milhoja grid backends will set this data in the same way and that backend developers should discover quickly what values to use.

While blocks are typically specified in Flash-X to have the same number of cells along each direction, the Grid unit does not insist that all blocks must be specified this way.  Therefore, this test was explicitly written to confirm correct functionality when using blocks with a different number of cells along each direction.

#### Success vs. Failure
As this test is a unit test it indicates via the Flash-X-standard `unitTest_0000` file if all tests passed or if any test failed.  Note that some expected values used to assess correctness are hardcoded in the code.  As a result, the test can only be configured in specific ways via the par files and the setup command (See the specifications below).  It is unlikely that this test will function correctly if more than one MPI process is used.

#### Status
This test has been verified to function correctly on GCE/compute-12 with the Intel 20.4 compiler and can therefore be included carefully in testsuites.  Please note that more tests will likely be added as more Milhoja grid backend functionality is added.  Please refer to the Doxygen `todo` information contained in the various test-specific files to understand how this test will likely grow.

#### GCE Testsuite Specifications
```
# The following test used to be configured with -debug.
UnitTest/Grid/AMR/Milhoja/1D/TestInit:
  setupOptions: -auto  -1d -nxb=8 +noio --with-unofficial=Grid/GridMain/AMR/Milhoja
  parfiles: test_milhoja_grid.par

# The following test used to be configured with -debug.
UnitTest/Grid/AMR/Milhoja/2D/TestInit:
  setupOptions: -auto  -2d -nxb=8 -nyb=4 +noio --with-unofficial=Grid/GridMain/AMR/Milhoja
  parfiles: test_milhoja_grid.par

# The following test used to be configured with -debug.
UnitTest/Grid/AMR/Milhoja/3D/TestInit:
  setupOptions: -auto  -3d -nxb=8 -nyb=4 -nzb=2 +noio --with-unofficial=Grid/GridMain/AMR/Milhoja
  parfiles: test_milhoja_grid.par
```
