# Helmholtz+WeakLib Hybrid EOS Unit Test

This unit tests the Helmholtz + WeakLib hybrid EOS implementation.  An grid of density, temperature, and electron fraction values provide the initial state for the EOS operating in `MODE_DENS_TEMP`.  The unit test will then cycle through additional modes (`MODE_DENS_EI` and `MODE_DENS_PRES`) and back around to `MODE_DENS_TEMP`, and check if the states obtained in the subsequent modes are within the provided tolerance of the initial "true" state.  By default, this unit test is restricted to only operate in the provided transition density range.

## Setting up the unit test

```bash
./setup unitTest/Eos/Helmholtz_Weaklib -auto +pm4dev -3d
```

If you are running the unit test outside of the object directory, you will need to also copy from/link to the files `helm_table.dat` and `wl-EOS-SFHo-15-25-50-noBCK.h5` located in the object directory.

## Notes

- A unit-test specific version of `Eos_putData.F90` changes some unused EOS variable maps to extract the inverted temperatures separately from Helmholtz and WeakLib
- Due to multiple incompatible versions of `Eos.h` existing in various subdirectories, the necessary one in `Hybrid` is included here
- `Helmholtz/Ye` contains a version of `Eos_wrapped` that does not call `Eos_getData`, so the necessary version in `EosMain` is included directly here.
