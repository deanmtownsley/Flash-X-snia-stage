# python script to run a simulation inside flashx container using maple API - requires python3.8+
# Refer to README.md for details

# import maple (python API version of maple)
import maple.api as maple

if __name__ == "__main__":
    # create an image object
    # name: name of the image
    # base: remote image of flashx environment
    # backend: docker/singularity
    image = maple.Image(name='grid',base='akashdhruv/flashx:latest',backend='docker')
    
    # create a container object
    # name: name of the local container
    # source: basedir (Flash-X directory)
    # target: path of mount directory inside the container (mount source to target)
    container = maple.Container(name='grid',target='/home/mount/FlashX')

    # build local image
    image.build()

    # execute commands inside the container
    # build and run amrex simulation
    container.run(image,'FlashTest -z /home/mount/FlashX \
                                   -o /home/mount/FlashX/container/unitTest/TestResults \
                                   -i /home/mount/FlashX/container/unitTest/unitTest.xml \
                                   UnitTest/Grid/AMR/AMReX/2d/Init \
                                   UnitTest/Grid/AMR/AMReX/2d/Refine \
                                   UnitTest/Grid/AMR/AMReX/2d/FluxCorrection \
                                   UnitTest/Grid/AMR/AMReX/2d/FluxCorrection2 \
                                   UnitTest/Grid/AMR/AMReX/2d/TestCyl \
                                   UnitTest/GridAnomalousRefine/2d/Paramesh \
                                   UnitTest/GridAnomalousRefine/2d/AMReX')

    # delete image
    image.delete()
