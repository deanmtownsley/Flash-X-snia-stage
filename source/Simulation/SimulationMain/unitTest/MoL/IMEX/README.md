# Method-of-Lines Multi-Rate Integrator Unit Test

- Solves the KPR test problem presented in [Chinomona & Reynolds, 2021](https://arxiv.org/abs/2007.09776) (Section 5.1) using the first- through third-order IMEX integration schemes
- This differs from the multi-rate unit test by including the fast terms as implicit terms (by default IMEX will use any fast terms as part of its explicit terms, so the equations here were modified to include the fast terms in the implicit update)
- The setup line for this unit test is (plus any additional site-specific options):
```bash
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX --parfile=flash.par.fbe
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX --parfile=flash.par.ssp2-222
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX --parfile=flash.par.ssp2-322
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX --parfile=flash.par.ssp2-332
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX --parfile=flash.par.ssp3-332
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX --parfile=flash.par.ssp3-433
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX --parfile=flash.par.ark-111
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX --parfile=flash.par.ark-121
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX --parfile=flash.par.ark-122
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX --parfile=flash.par.ark-222
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX --parfile=flash.par.ark-232
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX --parfile=flash.par.ark-233
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX --parfile=flash.par.ark-343
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX --parfile=flash.par.ark-443
```
