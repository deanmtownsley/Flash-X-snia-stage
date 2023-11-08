# Thornado Relaxation Benchmark

# uniform initial conditions, order 1
```
./setup Relaxation -auto -3d +cartesian +cube16 +thornado nE=16 swE=0 nSpecies=6 nNodes=2 nMoments=4 momentClosure=MINERBO thornadoOrder=ORDER_1 +weaklib -parfile=flash.par
```

# uniform initial conditions, order v/c
```
./setup Relaxation -auto -3d +cartesian +cube16 +thornado nE=16 swE=1 nSpecies=6 nNodes=2 nMoments=4 momentClosure=MINERBO thornadoOrder=ORDER_V +weaklib -parfile=flash.par
```

# get initial conditions from profile
```
./setup Relaxation -auto -1d +spherical -nxb=16 +thornado nE=5 swE=0 nSpecies=1 nNodes=2 nMoments=4 momentClosure=MINERBO thornadoOrder=ORDER_1 +weaklib -parfile=test_sph_1d.par
```
