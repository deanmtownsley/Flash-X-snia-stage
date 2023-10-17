import click


@click.command("fxprettify")
@click.argument("filelist", nargs=-1, required=True)
def fxprettify(filelist):
    """
    \b
    Flash-X prettify

    \b
    This simple command should be applied after using
    fprettify on Fortran files. The purpose of this
    command is to remove leading white space from each
    line to set emacs style formatting
    """
    # The logic is simple. Loop over files,
    # read lines, remove leading white space
    # if exists, and the rewrite the file
    for filename in filelist:

        with open(filename, "r") as ffile:
            lines = ffile.readlines()

        for index, line in enumerate(lines):
            if (line[0] == " ") and (line.strip()[:2] != "!!"):
                lines[index] = line[1:]

        with open(filename, "w") as ffile:
            ffile.writelines(lines)


if __name__ == "__main__":
    fxprettify()
