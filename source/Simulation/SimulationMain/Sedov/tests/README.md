## Sedov Comparison Test Designs
__TODO__: PUDs should write this.

#### Keywords
Grid, Milhoja, Pseudo-UG, Cartesian/2D, Cartesian/3D

#### Motivation
__TODO__: PUDs should write this.

Sedov comparison tests can be useful for testing and integrating Grid unit
implementations because it can be run with periodic boundary conditions for
short enough total simulation times.  This is helpful because running periodic
BCs can be implemented trivially for some AMR backends and, therefore,
full-featured BC functionality can effectively be ignored.

#### Success vs. Failure
Sedov comparison tests are expected to be deterministic.  Therefore, comparison
tests can be configured to demand that new results be bitwise identical to the
baseline.

It is expected that mass, linear momentums, and total energy are conserved.
The linear momentum for directions above NDIM should always be exactly zero.

__TODO__: Add in expected level of conservation for each conserved quantity?

