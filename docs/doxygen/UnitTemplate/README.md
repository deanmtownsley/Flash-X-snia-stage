This folder contains a mock Flash-X unit that demonstrates and explains the
Flash-X, doxygen-based inline documentation scheme.  Its intended users are
Flash-X developers and maintainers.  The content will be maintained up-to-date
so that developers can use the appropriate example file as a template for
creating a new file.  In the future, we hope to provide a set of scripts to
generate new files with correct boilerplate doxygen headers.

The scheme is motivated by the desire to make it as easy as possible to create
and actively maintain inline documentation.  It is co-designed with the Flash-X
PR review process to ensure that all doxygen-based inline documentation that is
in `main` is correct.  As part of this, the scheme distinguishes between
general, universal documentation that specifies the implementation-agnostic
public interface of the unit and implementation-specific documentation.

The former is written in the inline documentation of the stubs; the latter, in
the inline documentation of concrete implementations.  By making this
separation, the content of the inline documentation across a stub and all of its
associated concrete implementations are disjoint so that there is no need for
that documentation to be maintained synchronized.  Also, the amount of inline
documentation in the concrete implementations should be minimal.

Together the implementation-agnostic inline documentation and the Fortran
programmatic subroutine/function declaration in each stub fully specify the
interface of that subroutine/function.  They form a contract between those who
use the code and those who implement the code so that code use is decoupled from
implementation.  If all parties abide by the contract (including the developers
and maintainers of all implementations of the unit), then implementation changes
won't result in users being forced to update their code.

Developers can view the inline documentation of this mock unit by building the
Flash-X doxygen documentation using the developers' Doxyfile.
