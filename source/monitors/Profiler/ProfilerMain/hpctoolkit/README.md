## Using HPC ToolKit for granular profiling of Flash-X applications

This README only covers Flash-X specific API to enable granular
profiling of source code. To get a deeper understanding of best
practices on using sampling feature of HPC ToolKit we recommend 
reading its [documentation](http://hpctoolkit.org/).

To see a Flash-X specific example we recommend running experiments
in this [Lab Notebook](https://github.com/Lab-Notebooks/Flow-Boiling-Performance/blob/f1729bbfb5728904558a11277e167bbdb55c5063/simulation/FlowBoiling/Example2D/flashRun.sh#L1-L10).

### Setup and Compilation

This profiler is available by using the setup option, **+hpctoolkit**. 
This will require that you first install **hpctoolkit >= 2023.03.01** and define,

```make
LIB_HPCTOOLKIT = ${HPCTOOLKIT_PATH}/lib/hpctoolkit -lhpctoolkit
FFLAGS_HPCTOOLKIT = -g
```

in your site specific **Makefile.h**

### Usage

You can profile your Fortran source code in the following way,

```fortran
call Profiler_start("<profile-name>")
! .......
! code segment
! .......
call Profiler_stop("<profile-name>")
```

These subroutines interface with HPCToolKit API calls,

```cpp
void hpctoolkit_sampling_start(void);
void hpctoolkit_sampling_stop(void);
```

Flash-X specific API is designed to not allow multiple instances of
sampling to run concurrently. Therfore, applications will abort if a
new profiler is started without terminating the previous one. You can 
profile specific segements of code in the following way,

```fortran
call Profiler_start("<profile-segment-1>")
! .......
! first segment of code
! .......
call Profiler_stop("<profile-segement-1>")

! .......
! ignored segment of code
! .......

call Profiler_start("<profile-segment-2>")
! ......
! another segment of code
! ......
call Profiler_stop("<profile-segment-2>")
```

Please follow best practices to obtain realistic profiling results and
throughly read the [documentation](http://hpctoolkit.org/) before using 
this feature for granular profiling.

### Running Applications

When running this application we must tell hpcrun to initially turn
off sampling (it is on by default): use the -ds (or --delay-sampling)
option for hpcrun (dynamic) or set the `HPCRUN_DELAY_SAMPLING`
environment variable (static), along with additional options for trace
and event monitoring,

```bash
mpirun -n <procs> hpcrun -ds -t -e CYCLES -e CACHE-MISSES flashx
```

### Visualizing Output

If your application executes as expected and you do not face any permission
issues related to `pref_event_paranoid` then you should have a measurements
directory in your run directory. Next you need to execute,

```bash
hpcstruct <measurements-folder>
hpcprof <measurements-folder>
```

If these are succesfull a database directory will be created parallel to the 
measurements directory. To view results, you can execute,

```bash
hpcviewer <database-folder>
```
