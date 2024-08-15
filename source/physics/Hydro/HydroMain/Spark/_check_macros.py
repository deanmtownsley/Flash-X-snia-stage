#!/usr/bin/env python3
import re
from pathlib import Path


def get_spark_fortran_mc_files(root_path):

    mclist = [f.as_posix() for f in root_path.rglob("*.F90-mc")]

    return mclist


def main():
    spark_root = Path(__file__).parent.absolute()
    flashx_root = spark_root.parents[4]
    flashx_bin = flashx_root / "bin"

    macro_files = [f.as_posix() for f in spark_root.glob("*.ini")]
    mc_files = [f.as_posix() for f in spark_root.rglob("*.F90-mc")]
    mc_files += [f.as_posix() for f in spark_root.rglob("*.ini")]

    import sys
    sys.path.append(flashx_bin.as_posix())
    from macroProcessor import macroProcessor, macro_regex, invocation_regex

    m = macroProcessor()

    m.loadDefsList(macro_files)

    all_macros = set(m.keylist)
    used_macros = set()

    for file in mc_files:
        # print(f"processing {file}...")
        with open(file, 'r') as f:
            for line in f:
                # scan line for comments and strip them
                line_stripped = line
                inString = False
                for i, ch in enumerate(line):
                    if ch == '"':
                        inString = not inString
                    if ch == "!" and not inString:
                        line_stripped = line[0:i]
                        break

                # get macro invocation lines
                invocation_list = re.findall(macro_regex, line_stripped)

                # find macro keys used in mc file
                for invocation in invocation_list:
                    invocation_parts = re.match(invocation_regex, invocation)
                    used_macros.add(invocation_parts.group("key"))


    unused = all_macros - used_macros

    print("Following macros are defined but never uesed:")
    for macro in unused:
        print(f"\t{macro}")

if __name__ == "__main__":
    main()
