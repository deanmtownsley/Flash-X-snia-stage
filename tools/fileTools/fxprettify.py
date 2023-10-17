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
    # Store keywords for scope openers
    scope_keywords = ["program", "function", "subroutine", "module"]

    # Loop over files from filelist
    for filename in filelist:

        # Apply fprettify with 3 white space indentation
        fprettify.run([sys.argv[0], "-i", "3", filename])

        # Read lines from the formatted file
        with open(filename, "r") as ffile:
            fprettify_lines = ffile.readlines()

        # Create empty object for new lines and variables for storing
        # information related to previous line
        new_lines = []
        prev_line = ""
        prev_line_inscope = False

        # Loop over lines and adjust leading white spaces to satisfy emacs style formatting
        for index, line in enumerate(fprettify_lines):

            # Set scope flag to False for the current line
            curr_line_inscope = False

            # Check if previous line is attached to the opening scope
            if len(prev_line) > 0:

                # Test if line-continuation and scope keywords present in the pervious line
                if (prev_line[-1] == "&") and (
                    True in [prev_line[:6] == keyword[:6] for keyword in scope_keywords]
                ):
                    curr_line_inscope = True

                # Test if previous line is part of the opening scope and is continued
                elif (prev_line[-1] == "&") and (prev_line_inscope):
                    curr_line_inscope = True

            if (
                (line[0] == " ")
                and (line.strip()[:2] != "!!")
                and (not curr_line_inscope)
            ):
                line = line[1:]

            new_lines.append(line)
            prev_line = line.strip("\n")
            prev_line_inscope = curr_line_inscope

        # rewrite the file if lines have changed
        with open(filename, "w") as ffile:
            ffile.writelines(new_lines)


if __name__ == "__main__":
    fxprettify()
