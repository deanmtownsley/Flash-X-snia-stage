# python script to run a simulation inside flashx container using maple API - requires python3.8+
# Refer to README.md for details

# import maple (python API version of maple)
import maple.api as maple

if __name__ == "__main__":
    # create an image object
    # name: name of the image
    # base: remote image of flashx environment
    # backend: docker/singularity
    image = maple.Image(name='incomp_flow',base='akashdhruv/hypre:2.22.0',backend='docker')
    
    # create a container object
    # name: name of the local container
    # source: basedir (Flash-X directory)
    # target: path of mount directory inside the container (mount source to target)
    container = maple.Container(name='pool_boiling',target='/home/mount/FlashX')

    # build local image
    image.build()

    # execute commands inside the container
    # build and run paramesh simulation
    container.run(image,"./setup incompFlow/PoolBoiling -auto -2d -site=container \
                         +pm4dev -gridinterpolation=native -maxblocks=100 && \
                         cd object && make && grep 'setup_flashRelease =' setup_flashRelease.F90 && \
                         mpirun -n 1 ./flashx && cat unitTest_0000")

    # delete image
    image.delete()
