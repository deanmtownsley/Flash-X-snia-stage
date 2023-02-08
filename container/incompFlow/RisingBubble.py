# python script to run a simulation inside flashx container using maple API - requires python3.8+
# Refer to README.md for details

import click

# import maple (python API version of maple)
import maple.api as maple


@click.command(name="RisingBubble")
@click.argument("base_image", default="akashdhruv/amrex:latest")
def RisingBubble(base_image):
    """source/Simulation/SimulationMain/incompFlow/RisingBubble"""

    # create an image object
    # name: name of the image
    # base: remote image of flashx environment
    # backend: docker/singularity
    image = maple.Image(name="flashx", base=base_image, backend="docker")

    # build local image
    image.build()

    # execute commands inside the container
    # build and run amrex simulation
    maple.Run(
        "rising_bubble",
        image,
        "./setup incompFlow/RisingBubble -auto -2d -site=container -makefile=amd64 \
                      +amrex -maxblocks=100 && \
                      cd object && make && grep 'setup_flashRelease =' setup_flashRelease.F90 && \
                      mpirun -n 1 ./flashx && cat unitTest_0000",
    )

    # delete image
    image.delete()


if __name__ == "__main__":
    # call the command
    RisingBubble()
