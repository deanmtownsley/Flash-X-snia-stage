## Using HPC ToolKit for granular profiling of Flash-X applications

This README only covers Flash-X specific API to enable granular
profiling of source code. To get a deeper understanding of best
practices on basic usage and enabling sampling feature of HPC
ToolKit we recommend reading its documentation available at
(http://hpctoolkit.org/).

To see a Flash-X specific Makefile and Example we recommend following
updates on (https://github.com/Lab-Notebooks/Flow-Boiling-Performance)

### Setup and Compilation

This profiler is available by using the following setup shortcut `+hpctoolkit`

This will require that you first install hpctoolkit >= 2023.03.01
and define,

```
LIB_HPC_TOOLKIT = ${HPC_TOOLKIT_PATH}/lib/hpctoolkit -lhpctoolkit
```

and

```
FFLAGS_HPC_TOOLKIT = -g
```

in your site specific Makefile.h

### Usage

You can place calls to `Profiler_start(name)` and `Profiler_stop(name)`
to start and stop sampling in your FORTRAN source code. These subroutines
interface with HPC ToolKit API calls,

```
void hpctoolkit_sampling_start(void);
void hpctoolkit_sampling_stop(void);
```

Flash-X specific API is designed to not allow multiple instances of
sampling to run concurrently. Therfore, applications will abort if a
new  `Profiler_start` is called without terminating the previous one
using `Profiler_stop`.

See example under,

```
Simulation/SimulationMain/incompFlow/Driver_evolveAll
```

Please follow best practices to obtain realistic profiling results and
throughly read the documentation referenced above before using this
feature.

### Running Applications

When running this application we must tell hpcrun to initially turn
off sampling (it is on by default): use the -ds (or --delay-sampling)
option for hpcrun (dynamic) or set the `HPCRUN_DELAY_SAMPLING`
environment variable (static), along with additional options for trace
and event monitoring.

### Visualizing Output

If execution works out as expected and you do not face any permission
issues related to `pref_event_paranoid` then you should have a measurements
directory in your run directory. Next you need to execute,

```
hpcstruct <measurements-folder>
hpcprof <measurements-folder>
```

If these are succesfull a database directory will be created parallel to the 
measurements directory. You can execute,

```
hpcviewer <database-folder>
```

To view results

