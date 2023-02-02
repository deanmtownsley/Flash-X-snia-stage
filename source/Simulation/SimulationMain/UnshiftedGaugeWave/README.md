# Z4c Unshifted Gauge Wave

To setup this problem:
```bash
./setup UnshiftedGaugeWave -auto -3d +ug +nofbs +MoLERK --with-unofficial=physics/Spacetime
```

or

```bash
./setup UnshiftedGaugeWave -auto -3d +ug +nofbs +MoLMR --with-unofficial=physics/Spacetime
```

- Default `flash.par` is configured to run with four MPI ranks
- Included `plot.py` script can be ran to quickly check the solutions stability.  For the default `flash.par`:
```bash
python plot.py waveequation_hdf5_plt_cnt_000*
```

### TODO

- Add adjustable parameters for initial data amplitude/period/etc.
- Use as a Z4c unit test?
