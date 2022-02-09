# Flash-X

## Git/Testing Workflow

The current rules for collaborating via git are as follows

1.  Base all feature branches off of the main branch.
2.  When development on a feature branch is finished, merge the main branch onto the feature branch and run a test-suite on the feature branch.
3.  If at any point during the previous steps errors are discovered that need to
    be fixed, fix code in the feature branch, then go back to point 2 above. Repeat until no test errors.
4. If a merge conflict occurs when merging main into the feature branch _do not_ attempt to resolve conflicts using the  GitHub web interface - such an attempt can results in an unintended merge.
5.  Do not rebase a feature branch that has already been pushed to the GitHub
    repository.
6.  When you are ready, create a PR from the feature branch to the **staged** branch. If no conflicts occur in the **staged** branch, you are done. If conflicts occur, then depending on the extent and type of conflicts the resolution will be done on a case by case basis.

## Containerization Workflows

![incompFlow](https://github.com/Flash-X/Flash-X/workflows/incompFlow/badge.svg)
![Sod](https://github.com/Flash-X/Flash-X/workflows/Sod/badge.svg)
![Sedov](https://github.com/Flash-X/Flash-X/workflows/Sedov/badge.svg)

These workflows are located in `.github/workflows` and are not part of default testing framework. Please to refer `.github/workflows/README.md` and `container/README.md` for details on containerization with **Flash-X**

## Libraries

The code of the PARMESH library is currently integrated into the GridMain implementation as a
subtree and not organized or built as a separate library.

Some applications and tests use external libraries that are expected to be already installed on the
system where Flash-X is being built and run. The directory locations of such library installations
should be made know to the Flash-X build system by a site-specific (or, as a fallback, OS-specific)
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

## Tests 
You can obtain the source code for flashtest and a full set of tests from the 
Flash-X-test repository. The repository also has tools to help you setup your local testsuite.

