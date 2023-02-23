__IMPORTANT: This document has been constructed as a prototype design document.  The contents have not been verified as correct.__

## AMReX Refinement Unit Test Design

This document is intended to be a living document that captures and communicates the __current__ design and use of the AMReX `TestRefine` unit test.  In this sense, users can assume that this document and all the contents of this folder are correct, up-to-date, and an encapsulation of the official Flash-X low-level test related to AMReX refinement.

#### Keywords
AMReX, Grid, Refinement, Cartesian, Unit Test, 2D

#### Motivation
This unit test was developed via Test-Driven Development (TDD) of the AMReX Grid unit implementation.  Rather than serve for automatic unit testing as part of continuous integration, the original intent was to allow for manually testing appropriate refinement with AMReX since AMReX refinement/derefinement was black box at the time that the test was created.  Specificallly, it allows users to control and therefore experiment with how AMReX mesh refinement and derefinement is working as implemented in Flash-X.  The code advances the data in UNK in time manually.  At each step, the code sets all data in UNK to zero except for possibly at a few points, whose non-zero values define what level of refinement must be achieved in the blocks that contain them.

#### Design
Unittest is configured such that
* the coarsest level has 2 blocks along both the X and Y dimensions and each block has 8 cells along both dimensions
* refinement/derefinement is done on every other time step
* only a single, cell-centered physical variable is managed.

Changes to the physical data are managed on a step-by-step basis and is done so in conjunction with a custom `gr_markRefineDerefineCallback` routine so that the non-zero data values in the physical data specify the refinement level to be achieved by AMReX for the block containing that point.

Time Stepping
* [Init](TestRefine_Init_Both.pdf) - One data point that refines down to level 3 and one to level 2.  After initial refinement, note that the lower-/upper-right blocks should be level 1.  However, they are promoted to level 2 due to periodic BC.
* [Steps 1/2](TestRefine_Step2_Both.pdf) - Set all data to zero to invoke full derefinement to level 1.
* [Steps 3/4](TestRefine_Step4_Both.pdf) - One data point near corner to invoke refinement to level 2 in its block only.
* [Steps 5/6](TestRefine_Step6_Both.pdf) - Same point to invoke refiment to level 5.  Should only achieve level 3 refinement under point.  Should see level 2 refinement elsewhere to maintain 1-level difference at refinement boundaries.
* __Steps 7/8__ - No change to data.  Let refinement achieve level 4 under point.  Due to periodic boundary conditions, refinement to level 3 on other three corners.
* [Steps 9/10](TestRefine_Step8_Both.pdf) - No change to data.  Refinement should be unchanged since the maximum refinement level is 4.
* __Steps 11/12__ - Add another point set to refine to level 4.  No refinement check here.  We just let this step advance toward full refinement.
* [Steps 13/14](TestRefine_Step14_Both.pdf) - No change to data.  Final refinement/derefinement should be achieved at this step.
* __Steps 15/16__ - Set all data to zero to invoke full derefinement to level 1.
* [Steps 17/18](TestRefine_Step18_Both.png) - Tag two neighboring cells that are divided by a block boundary to confirm that both blocks are refined in true octree fashion as opposed to having a single, refined block that is translated to contain both tagged cells.  This is an example of sacrificing grid efficiency for the organizational simplicity of octree.  A patch-based scheme might have refined a single block that is translated to include both of the tag cells.

#### Success vs. Failure
As this test is a unit test it indicates via the Flash-X-standard `unitTest_0000` file if all tests passed or if any test failed.  Note that some expected values used to assess correctness are hardcoded in the code.  As a result, the test can only be configured in specific ways via the par files and the setup command (See the specifications below).  It is unlikely that this test will function correctly if more than one MPI process is used.

#### Status
This test has been verified to function correctly on GCE/compute-12 with the Intel 20.4 compiler and can therefore be included carefully in testsuites.  Since the motivation for this test was as an initial tool for exploring AMReX through Flash-X, the need for inclusion in a testsuite is not obvious.  On one hand, the test required a fair amount of novel code and therefore could require more maintenance.  On the other hand, it might help us identify immediately if changes in AMReX violate some of our assumptions about AMReX's refinement.  Ideally, we should work with the AMReX team (if not already done so) to document the refinement so that it is not black box to Flash-X Grid developers or Flash-X users.

#### GCE Testsuite Specifications
```
UnitTest/Grid/AMR/AMReX/2d/Refine:
  setupOptions: -auto -2d -nxb=8 -nyb=8 +noio +amrex
  parfiles: test_amrex_grid.par
```

#### Conclusions Derived from Experimention
If a cell is tagged for refinement in a fine level, but not tagged in any lower levels, the cell's block is correctly refined.

When a block is refined, AMReX does not check the data at this new level for additional refinement.
