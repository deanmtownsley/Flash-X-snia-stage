__IMPORTANT__: The contents of this folder and its subfolders should *not* be
included in this repository.  However, this was written at a time when the
offline toolchain functionality was not available and was therefore written by
some poor soul by hand.  Once the offline toolchain can automatically generate
this content as part of official Flash-X offline processes, this content should
be removed from the repository.

The GPU version of the code requires that the cell-centered data be ordered in
memory in a way different from how the setup tool chooses to order them.
Therefore, after running setup but before compiling, users must copy the correct
Simulation_*D.h file from here to their object directory.

