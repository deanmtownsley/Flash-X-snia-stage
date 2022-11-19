## DomainBBoxRadians Unit Test Design

This document is intended to capture and communicate the design and use of the `DomainBBoxRadians` unit test.

#### Keywords
Paramesh, AMReX, Milhoja, Uniform, Grid, Interface, Pseudo-UG, Cartesian, Unit Test, Geometry, Curvilinear, Angle Coordinate, 3D

#### Motivation
This test was designed largely to check whether angle coordinates are consistently converted from degrees (in Grid-owned runtime parameters) to radians (Flash-X code internal convention). All Grid implementations should have consistent behavior in this respect. Some other Grid features are also tested along the way, mostly those that plly to all, or nearly all, Grid implementations.

While blocks are typically specified in Flash-X to have the same number of cells along each direction, the Grid unit does not insist that all blocks must be specified this way.  Therefore, this test was explicitly written to confirm correct functionality when using blocks with a different number of cells along different directions.

#### Success vs. Failure
As this test is a unit test it indicates via the Flash-X-standard `unitTest_0000` file if all tests passed or if any test failed.  Note that some expected values used to assess correctness are hardcoded in the code.  As a result, the test can only be configured in specific ways via the par files and the setup command (See the specifications below).  It is unlikely that this test will function correctly if more than one MPI process is used.

#### Status
Initial adaptation from other Grid unitTest in progress; done for test_pm and test_ug.

#### Unofficial GCE Testsuite Specifications
```
[test_pm]
    testNode = "UnitTest/Grid/pm4/3dcyl/DomainBBoxRadians"
    setupOptions = "-auto -strict -3d -nxb=8 -nyb=8 -nzb=10 +noio -gridinterpolation=native"
    numProcs = 1
    parFile = "test_pm_grid.par"

[test_ug]
    testNode = "UnitTest/Grid/UG/3dcyl/DomainBBoxRadians"
    setupOptions = "-auto -strict -3d +cylindrical -nxb=8 -nyb=8 -nzb=10 --without-unit=physics/Eos +noio +ug +nofbs"
    numProcs = 1
    parFile = "test_ug_grid.par"

[test_amrex]
    testNode = "UnitTest/Grid/AMReX/3dcyl/DomainBBoxRadians"
    setupOptions = "-auto -3d -nxb=8 -nyb=8 -nzb=8 +noio +amrex"
    numProcs = 1
    parFile = "test_pm_grid.par"

[test_milhoja]
    testNode = "UnitTest/Grid/Milhoja/3dcyl/DomainBBoxRadians"
    setupOptions = "-auto -3d -nxb=8 -nyb=8 -nzb=8 +noio -with-unofficial=Grid/GridMain/AMR/Milhoja"
    numProcs = 1
    parFile = "test_pm_grid.par"
```
