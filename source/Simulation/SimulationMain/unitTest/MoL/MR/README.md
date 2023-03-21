# Method-of-Lines Multi-Rate Integrator Unit Test

- Solves the KPR test problem presented in [Chinomona & Reynolds, 2021](https://arxiv.org/abs/2007.09776) (Section 5.1) using the third-order IMEX-MRI-GARK3b integration scheme
- The IMEX-MRI-GARK4 scheme currently will not pass this unit test to a sufficient enough level of accuracy.  [Chinomona & Reynolds, 2021](https://arxiv.org/abs/2007.09776) does not recommend the use of this scheme due to an overly-restrictive joint-stability region.
- The setup line for this unit test is (plus any additional site-specific options):
```bash
./setup unitTest/MoL/MR -auto -1d +ug +nofbs +MoLMR
```
