# Spacetime Static Unit Test

## Setup

Setup the test as follows

```bash
./setup unitTest/Spacetime/static -auto -1d +nofbs +noio --with-unofficial=physics/Spacetime
```
```bash
./setup unitTest/Spacetime/static -auto -2d +nofbs +noio --with-unofficial=physics/Spacetime
```
```bash
./setup unitTest/Spacetime/static -auto -3d +nofbs +noio --with-unofficial=physics/Spacetime
```

Other Grid-backends besides UG (which this defaults to) may be used, but additional refinement parameters will be necessary.  This functionality may be added here in the future along side tests of static non-vacuum backgrounds when relativistic hydro and/or radiation solvers are available.

## Test Description

This unit test will generate a Schwarzschild solution in Cartestion Kerr-Schild coordinates for a black hole with $M = 1 M_\odot$.  The ADM variables will be set according to (see Table 2.1 in [1])

### Lapse Function
$$
   \alpha = \left(1 + \frac{2M}{r}\right)^{-1/2}
$$

### Shift Vector
$$
   \beta^i = \frac{2M}{r}\alpha^2\ell^i
$$

### Spatial Metric
$$
   \gamma_{ij} = \eta_{ij} + \frac{2M}{r}\ell_i\ell_j
$$

### Extrinsic Curvature
$$
   K_{ij} = \frac{2M\alpha}{r}\left[\eta_{ij} -\left(2 + \frac{2M}{r}\right)\ell_i\ell_j\right]
$$

where $r$ is the radial coordinate, $\eta_{ij} = \delta_{ij}$ is a Minkowski metric with Cartesian coordinates $x^i$, and $\ell^\mu = (-1, x^i/r)$ is a null four-vector (with respect to the four-metric $g_{\mu\nu}$).

The provided `flash.par` file will setup the grid to avoid having the origin at $r = 0$ at a cell center.

After setting this initial data, the unit test will verify that these ADM variables were set and obtained correctly via the accessor procedures available in the Spacetime unit.


[1] Baumgarte, T. W., & Shapiro, S. L. 2010, _Numerical relativity: solving Einstein’s equations on the computer_ (Cambridge ; New York: Cambridge University Press)
