# python script to run a simulation inside flashx container using maple API - requires python3.8+
# Refer to README.md for details

import click

# import maple (python API version of maple)
import maple.api as maple


@click.command(name="DeformingBubble")
@click.argument("base", default="flashx")
def DeformingBubble(base):
    """source/Simulation/SimulationMain/incompFlow/DeformingBubble"""

    # create an image object
    # base: remote image of flashx environment
    image = maple.Image(name=base, base=base)

    if image.platform == "linux/amd64":
        makefile = "amd64"
    elif image.platform == "linux/ppc64le":
        makefile = "ppc64le"
    else:
        print(f"Cannot compile on platform: {image.platform}")
        raise NotImplementedError()

    # build local image
    maple.Publish(
        "deforming_bubble",
        image,
        f"./setup incompFlow/DeformingBubble -auto -3d -nxb=16 -nyb=16 -nzb=16 \
                         -site=container -makefile={makefile} --index-reorder +pm4dev -maxblocks=200 \
                         -gridinterpolation=native --without-unit=Grid/GridSolvers/HYPRE \
                         +noio +nolwf Bittree=True && mkdir -pv /home/run && cd object && \
                         make -j && cp flashx /home/run",
    )


if __name__ == "__main__":
    # call the command
    DeformingBubble()
