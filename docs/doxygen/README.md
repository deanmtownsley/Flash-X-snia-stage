This folder contains
* two doxygen configurations and
* a mock Flash-X unit that demonstrates the Flash-X inline documentation scheme.

The configuration file `Doxyfile_API` will create doxygen content only for the
public interfaces of each unit, which is typically derived from the inline
documentation written in the stub files.  This content is what is made publicly
available through the Flash-X website and is intended for Flash-X users.

The configuration file `Doxyfile_developers` will create content for all
available doxygen inline documenation.  This includes the documentation in the
concrete implementation files, which is intended for developers and maintainers.

The official doxygen pages available through the Flash-X website are created
automatically with every new commit in the `main` branch.  Users can create
local doxygen pages with the command
```
doxygen <Doxyfile name>
```

Refer to the README.md in `UnitTemplate` for more information on the
documentation scheme.
