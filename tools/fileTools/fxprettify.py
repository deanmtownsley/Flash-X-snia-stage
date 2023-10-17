import sys
import click
import fprettify


@click.command("fxprettify")
@click.argument("filelist", nargs=-1, required=True)
def fxprettify(filelist):
    """
    \b
    Flash-X prettify

    \b
    This command applies fprettify with 3 white space
    indentation on Fortran files and then removes leading
    white spaces to set emacs style formatting
    """
    # Loop over files from filelist
    for filename in filelist:

        # Apply fprettify with 3 white space indentation
        fprettify.run([sys.argv[0], "-i", "3", filename])

        # Read lines from the formatted file
        with open(filename, "r") as ffile:
            lines = ffile.readlines()

        # Loop over lines and adjust leading white spaces
        # to satisfy emacs style formatting
        for index, line in enumerate(lines):
            if (line[0] == " ") and (line.strip()[:2] != "!!"):
                lines[index] = line[1:]

        # rewrite the file
        with open(filename, "w") as ffile:
            ffile.writelines(lines)


if __name__ == "__main__":
    fxprettify()
