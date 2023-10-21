# Spacetime Unit

The Spacetime unit houses solvers responsible for evolving the dynamic geometry of spacetime as governed by the Einstein's field equations of general relativity (GR)

$$
   G_{\mu\nu} = 8\pi T_{\mu\nu}
$$

where $G_{\mu\nu}$ is Einstein tensor describing the curvature of spacetime, and $T_{\mu\nu}$ is the total stress-energy tensor describing the matter and energy composition of the spacetime. For simplicity, the `Spacetime` unit will assume units of $G = c = 1$.

While there is no requirement on how any particular Spacetime unit implementation solves these equations, this unit will provide access to a discretization of the full 4D spacetime into a 3+1 split of 3D spacelike hypersurfaces at constant time slices.  In particular, the Arnowitt, Deser and Misner (ADM) form of the invariant line element will be assumed

$$
   ds^2 = -\alpha^2 dt^2 + \gamma_{ij}\left(dx^i + \beta^i dt\right)\left(dx^j + \beta^j dt\right)
$$

where $\alpha$ is the lapse function which measures the proper time between spacelike hypersurfaces, $\beta^i$ is the shift vector representing the shift in coordinates from one spacelike hypersurface to the next, and $\gamma_{ij}$ is the metric on a spacelike hypersurface.  The extrinsic curvature, $K_{ij}$, of a spacelike hypersurface will also be provided.  See Chapter 2 of [1] for additional information on the 3+1 decomposition.

All implemented solvers will be expected to fill the provided grid variables with the above ADM quantities (the grid variables are defined at the level of `SpacetimeMain`).  By default, requiring use of the `Spacetime` unit without requesting a specific solver will only provide access to the ADM variables, and it will be the responsibility of the user to properly set the initial data for these variables, and to update the variables as needed during a simulation.


[1] Baumgarte, T. W., & Shapiro, S. L. 2010, _Numerical relativity: solving Einstein’s equations on the computer_ (Cambridge ; New York: Cambridge University Press)