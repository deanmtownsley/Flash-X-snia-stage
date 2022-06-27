## Milhoja Flux Data Access Unit Test Design

This document is intended to be a living document that captures and communicates
the __current__ design and use of the Milhoja `TestFluxes` unit test.  In this
sense, users can assume that this document and all the contents of this folder
are correct, up-to-date, and an encapsulation of the official Flash-X low-level
tests related to Milhoja flux data access.

#### Keywords
Milhoja, Grid, Interface, Flux, Cartesian, Unit Test, 1D, 2D, 3D

#### Motivation
This unit test was developed via Test-Driven Development (TDD) of the Milhoja grid backend's Fortran/C++ interoperability layer and the Flash-X Milhoja Grid unit implementation.  In particular, it was useful for building up functionality without having to develop a means to write flux data to an external file for review.  Since it is a simple and quick unit test and its checks are generally useful, there is no immediate need to remove the test.  Indeed, it might be invaluable for quick and correct TDD-based development of subsequent Milhoja grid backends.  

This test was designed largely to confirm basic read and write access of flux data using `Grid_iterator_t` and `Grid_tile_t`.  
This test also confirms the inability to access flux data for dimension above NDIM.  Therefore, including [123]D versions of this test in the testsuite seems good and useful.

#### Success vs. Failure
As this test is a unit test it indicates via the Flash-X-standard `unitTest_0000` file if all tests passed or if any test failed.  It is unlikely that this test will function correctly if more than one MPI process is used.

#### Status
This test has been verified to function correctly on GCE/compute-12 with the Intel 20.4 compiler and can therefore be included carefully in testsuites.

