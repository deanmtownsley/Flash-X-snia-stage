# python script to run a simulation inside flashx container using maple API - requires python3.8+
# Refere to README.md for details

# import maple (python API version of maple)
import maple.api as maple

if __name__ == "__main__":

    # create a flashx object 
    # base: remote image of flashx environment
    # container: name of the local container
    # source: basedir (Flash-X directory)
    # target: path of mount directory inside the container (mount source to target)
    flashx = maple.Maple(base='akashdhruv/flash:latest',container='rising_bubble',
                         target='/home/mount/Flash-X',backend="singularity")

    # build local image from remote image
    flashx.image.build('incompFlow')
    #
    flash.container.pour('incompFlow')
    #
    # execute commands inside the container
    # build and run amrex simulation
    flashx.container.execute("./setup incompFlow/RisingBubble -auto -2d -site=/home/site \
                              +amrex -maxblocks=100 && \
                              cd object && make && grep 'setup_flashRelease =' setup_flashRelease.F90 && \
                              mpirun -n 1 ./flashx && cat unitTest_0000")

    # clean local container
    flashx.container.rinse()
    #
    # remove (delete) instance of remote image from local machine
    flashx.image.delete('incompFlow')
