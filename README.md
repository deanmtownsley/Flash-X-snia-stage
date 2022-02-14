# Flash-X

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
   PRs have been merged into the **staged** branch.
8. If the test suite passes, a composite PR will be created from
   **staged** into **main**, and you won't have to do anything more.
9. If the test suite fails, it is expected that you will resolve the
   failures immediately. If the failures continue over several
   iterations, or if the conflicts prove to be non-trivial, the
   resolution will involve someone designated by the Council to work
   with you.

Contributors with write permission should create a feature branch
instead of a fork. The remainder of the workflow remains the same.

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

