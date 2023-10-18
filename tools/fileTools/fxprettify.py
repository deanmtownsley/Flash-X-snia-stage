from io import StringIO
import sys
import click
import fprettify


class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio  # free up some memory
        sys.stdout = self._stdout


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

        # Read lines from the original file
        with open(filename, "r") as ffile:
            original_lines = ffile.readlines()

        # Apply fprettify with 3 white space indentation
        with Capturing() as fprettify_lines:
            fprettify.run([sys.argv[0], "-i", "3", "-s", filename])

        # Create empty object for new lines and variables
        # for storing information related to previous line
        new_lines = []
        prev_line = ""
        prev_line_indent = True

        # Loop over lines and adjust leading white
        # spaces to satisfy emacs style formatting
        for line in fprettify_lines:

            line = line + "\n"

            # Set line indent to true and then perform
            # tests to see if indentation is needed
            line_indent = True

            # If previous line is > 0 perform indentation tests
            if len(prev_line) > 0:

                # Check if previous line is continued and does
                # not start with white space or comment symbol
                if (prev_line[-1] == "&") and (prev_line[0] not in [" ", "!"]):
                    line_indent = False

                # Check if previous line is continued and
                # does not have indentation set to False
                elif (prev_line[-1] == "&") and (not prev_line_indent):
                    line_indent = False

            # Now check if indentation needs to be applied to
            # the current line based on desired conditions
            if (line[0] == " ") and (line.strip()[:2] != "!!") and (line_indent):
                line = line[1:]

            new_lines.append(line)
            prev_line = line.strip("\n")
            prev_line_indent = line_indent

        # rewrite the file if lines have changed
        if new_lines != original_lines:
            with open(filename, "w") as ffile:
                ffile.writelines(new_lines)


if __name__ == "__main__":
    fxprettify()
