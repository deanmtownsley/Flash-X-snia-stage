# python script to run a simulation inside flashx container using maple API - requires python3.8+
# Refer to README.md for details

# import pymaple (python API version of maple)
import pymaple

if __name__ == "__main__":

    # create a flashx object 
    # image: remote image of flashx environment
    # container: name of the local container
    # source: basedir (Flash-X directory)
    # target: path of mount directory inside the container (mount source to target)
    flashx = pymaple.Maple(image='akashdhruv/flash:latest',container='pool_boiling',target='/home/mount/Flash-X')

    # build local image from remote image
    flashx.build()
    #
    # pour local image to local container
    flashx.pour()

    # execute commands inside the container
    flashx.execute("export OMPI_MCA_btl_vader_single_copy_mechanism=none")

    # build and run amrex simulation
    flashx.execute("./setup incompFlow/PoolBoiling -auto -2d -site=/home/site \
                    +amrex -maxblocks=100 && \
                    cd object && make && grep 'setup_flashRelease =' setup_flashRelease.F90 && \
                    mpirun -n 1 ./flash5 && cat unitTest_0000")

    # build an run paramesh simulation
    #flashx.execute("./setup incompFlow/PoolBoiling -auto -2d -site=/home/site \
    #                +pm4dev -gridinterpolation=native -maxblocks=100 && \
    #                cd object && make && grep 'setup_flashRelease =' setup_flashRelease.F90 && \
    #                mpirun -n 1 ./flash5 && cat unitTest_0000")

    # rinse (stop and delete) local container
    flashx.rinse()
    #
    # clean (delete) local image
    flashx.clean()
    #
    # remove (delete) instance of remote image from local machine
    flashx.remove()
