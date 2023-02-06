# Method-of-Lines Implicit-Explicit Integrator Unit Test

- Solves the advection-reaction equation:
$$
\partial_t u + \alpha \partial_x u = -\beta u\\
u(x,t=0) = A \exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right),\ \ \ x\in [0,2]
$$
- The setup lines for this unit test is (plus any additional site-specific options):
```bash
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX --parfile=flash.par.fbe
```

## TODO

- Add method-specific order accuracy checks
- Debug IMEX-SSP methods (less than expected accuracy for explicit portion currently)
