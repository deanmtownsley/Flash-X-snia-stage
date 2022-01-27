# python script to run a simulation inside flashx container using maple API - requires python3.8+
# Refer to README.md for details.

# import maple (python API version of maple)
import maple.api as maple

if __name__ == "__main__":
    # create an image object
    # name: name of the image
    # base: remote image of flashx environment
    # backend: singularity/docker
    image = maple.Image(name='incomp_flow',base='docker://akashdhruv/flashx:latest',backend='singularity')    

    # create a container object
    # name: name of the local container
    # source: basedir (Flash-X directory)
    # target: path of mount directory inside the container (mount source to target)
    container = maple.Container(name='lidcav',target='/home/mount/FlashX')

    # build local image
    image.build()

    # Run FlashTest suite inside the container
    container.run(image,'FlashTestSuite -z /home/mount/FlashX \
                                        -o /home/mount/FlashX/TestResults \
                                        -s container \
                                        UnitTest/LidDrivenCavity/AMReX/2d')

    # delete image
    image.delete()
