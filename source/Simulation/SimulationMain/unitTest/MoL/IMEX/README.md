# Method-of-Lines Multi-Rate Integrator Unit Test

- Solves the KPR test problem presented in [Chinomona & Reynolds, 2021](https://arxiv.org/abs/2007.09776) (Section 5.1) using the first- through third-order IMEX integration schemes
- This differs from the multi-rate unit test by including the fast terms as implicit terms (by default IMEX will use any fast terms as part of its explicit terms, so the equations here were modified to include the fast terms in the implicit update)
- The setup line for this unit test is (plus any additional site-specific options):
```bash
./setup unitTest/MoL/IMEX -auto -1d +ug +nofbs +MoLIMEX
```
