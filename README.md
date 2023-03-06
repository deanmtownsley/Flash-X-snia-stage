# Flash-X, a Multiphysics Scientific Software System

Flash-X is a highly composable multiphysics software system that can
be used to simulate physical phenomena in several scientific
domains. It is derived from FLASH, which has a history of being a
community code for several communities. The Flash-X architecture has
been redesigned to be compatible with increasingly heterogeneous
hardware platforms. Part of the redesign is a newly designed
performance portability layer that is language agnostic.

## Libraries

The code of the PARAMESH library for Flash-X is provided as a git submodule by a separate code
repository named PARAMESH-for-Flash-X.  After successful initialization of the submodule
--- see further below ---
the code will appear as a subtree under `Grid/GridMain/AMR` in the `PM4_package` directory.
The Flash-X `setup` and build mechanisms will take care of compiling this code (if needed
for a simulation configuration); PARAMESH code is not organized or built as a separate library.

Some applications and tests use external libraries that are expected to be already installed on the
system where Flash-X is being built and run. The directory locations of such library installations
should be made known to the Flash-X build system by a site-specific (or, as a fallback, OS-specific)
Makefile.h file. See the subdirectories under sites/ .

This applies in particular to the AMReX library. Separate library instances for 1D, 2D, and 3D
should be installed, and the appropriate locations mentioned in Makefile.h .

On the other hand, some applications and tests use INTERNAL libraries. Such libraries are built, as part of the Flash-X setup process, from source code that is located in subdirectories under lib/ . There are two cases for such libraries:

1. **Library source code is included as part of the Flash-X git respository.**

   An example is the sqrt3 library, whose source code is included in lib/sqrt/ .

2. **Library source code must be retrieved from a separate repository.**

   Examples are the THORNADO and WEAKLIB libraries.
   Follow the instructions on submodules to automatically put the source code for these two
   in the right places in subdirectories under lib/.

## Git with Submodules

To prepare for building simulations that use libraries whose code must be retrieved
from a separate git repository, the following modified `git` commands can be used:

- `git pull --recurse-submodules=yes` (in place of the usual `git pull`)
- `git submodule update --init` (additionally, after `git pull`)
- `git submodule update --init source/Grid/GridMain/AMR/Paramesh4/PM4_package`
  (variant of the previous item, in case you want to only get the Paramesh package submodule)

## Git/Testing Workflow

The current rules for collaborating via GitHub are as follows:

Contributors with
read only permission to the Flash-X code repository should use the following
guidelines to create a pull request:

1. Create a fork.
2. Make your changes.
3. Create a PR to the **staged** branch whenever you wish.
   Give your PR a title that begins with the word "DRAFT".
   This will allow any discussion about the pull
   request to be conducted on github.
4. When you are ready for the pull request to be accepted, merge from **main**
   into your forked code, to ensure that your fork is not out of sync.
4. If a merge conflict occurs when merging **main** into the feature branch,
   _do not_ attempt to resolve conflicts using the  GitHub web interface - such an attempt can results in an unintended merge to **main**.
5. Run a local version of your test suite and make sure everything
   passes.
6. Make sure your latest commit has been pushed.
7. Remove "DRAFT" from your pull request name. If no further problems
   are found, this will cause the PR
   to be merged. The test suite is run at night if one of more
   PRs have been merged into the **staged** branch. PRs that come in
   before 6PM CST are more likely to be included in that night's test.
   Monitor the repo to see whether your PR was merged and the test suite passed.
   A comment will be added to your PR if the test suite failed.
8. If the test suite passes, a composite PR will be created from
   **staged** into **main**, and you won't have to do anything more. This will
   likely happen the day the test suite passed.
9. If the test suite fails, it is expected that you will prioritize resolving the
   failure. Note that the merged and colliding code will be available in the staged branch.
   You can copy that code into a local working copy to resolve the issue. **Please note that you
   should never make any commit into the staged branch.**
   If the test suite passes, you can reissue a PR and ask for a test suite run by leaving a comment in the PR.
   If the failures continue, we abandon the stage branch at the end of the day, and PRs have to be created again.
   If we determine that the interoperability is compromised, someone from the core team might have
   to get involved to help resolve.

Contributors with write permission should create a feature branch from the main branch
instead of a fork. The remainder of the workflow remains the same.

## Code Formatting
Use [fprettify](https://github.com/pseewald/fprettify) to format your source code and enforce consistency with indentation and whitespacing. Works for both `.F90` and `.F90-mc` files

## Special Syntax Highlighting

To enable syntax for **Flash-X** specific keywords implement following settings based on the editor:

### **VIM**

- Create `$HOME/.vim/after/syntax/fortran.vim` and add:

```
" Custom FORTRAN keywords
syn match customKeyword '\<*\zs@MACRO\>'
syn match customKeyword '\<*\zs@M\>'
syn match customKeyword '\<*\zs!!REORDER\>'
syn match customKeyword '\<*\zs!!NOVARIANTS\>'

" Definitions
hi def link customKeyword Keyword
```

- Create `$HOME/.vimrc` and add:
```
" Turn on syntax
syntax on

" Set file type
autocmd BufNewFile,BufRead *.F90-mc set filetype=fortran
```

## Containerization Workflows

[comment]: ![incompFlow](https://github.com/Flash-X/Flash-X/workflows/incompFlow/badge.svg)
[comment]: ![Sod](https://github.com/Flash-X/Flash-X/workflows/Sod/badge.svg)
[comment]: ![Sedov](https://github.com/Flash-X/Flash-X/workflows/Sedov/badge.svg)

These workflows are located in `.github/workflows` and are not part of default testing framework. Please to refer `.github/workflows/README.md` and `container/README.md` for details on containerization with **Flash-X**

## Tests
Test specifications for individual simulations are included under ``*/tests/tests.yaml`` files in each simulation directory under ``source/Simulation/SimulationMain``. New tests should be added as enteries in the prescribed YAML format before including it as a part of suites on different platforms.

Testing and maintainence of the code is implemented using command line tools available in Flash-X-Test repository: https://github.com/Flash-X/Flash-X-Test

Please refer to the instructions there to setup your own testing infrastructure. Also take a look at ``sites/ganon_jenkins/UnitTests.suite`` for an example of publicly available test suite which can be edited to enable code coverage for new modules.

Testing servers:

- Argonne, GCE:

  FlashTest server for running tests on `staged` branch - https://jenkins-gce.cels.anl.gov/job/Flash-X-staged

  FlashTestView - https://web.cels.anl.gov/projects/FLASH5/testsuite/home.py

- Ganon:

  FlashTest - http://ganon2.device.utk.edu:8080 
