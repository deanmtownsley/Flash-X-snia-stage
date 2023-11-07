## Public testing for Flash-X

The current location of the jenkins server that is running tests is: http://ganon2.device.utk.edu:8080

## Adding new tests

Test configurations are recorded in `testsuites/*.suite` files which divide the tests based on their
respective coverage. You can add new tests in existing files or create a new one to increase the
coverage.

## Using benchmarks

Running tests that compare against an approved benchmark can be implemented by adding comparison
tests to `testsuites/*.suite` and adding the date stamp of the approved benchmark using `-cbase 0000-00-00`.

The benchmarks should be available as compressed tar files. See following files for example:

https://anl.box.com/s/3ic5yepg7ad2m98g6bf4grjynep4aj88

This tar archive should be unpacked into `/archive/ganon_jenkins/0000-00-00`. This is done to ensure that the size 
of benchmarks stays under control. If benchmarks are updated, a new tar file should be created along with addition 
of new lines in Dockerfile like the one below:

`RUN curl -s -L https://anl.box.com/shared/static/74yz29aufgc46v8t78r0pligdwm4cn6e | \
 tar xvz --strip-components=1 -C $ARCHIVEDIR`

The corresponding test for this archive are located in `testsuites/IncompNS.suite`
