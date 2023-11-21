#!/usr/bin/env python

import os
import sys
import shutil

from pathlib import Path

import milhoja

if __name__ == "__main__":
    # ----- HARDCODED VALUES
    DESTINATION = Path.cwd()
    INDENT = 4
    # Allow overwrite for now since code generation is not yet built into
    # Flash-X setup system and generated code is in version control system.
    NO_OVERWRITE = True

    # ----- COMMAND LINE ARGUMENTS
    if len(sys.argv) != 3:
        print()
        print("Please specify NDIM and location of Milhoja library installation")
        print()
        exit(1)

    ndim = int(sys.argv[1])
    milhoja_path = Path(sys.argv[2]).resolve()

    assert ndim in [2, 3]
    tf_json = f"gpu_tf_hydro_{ndim}D.json"
    src_json = milhoja_path.joinpath("include", f"sizes_{ndim}D.json")
    dst_json = milhoja_path.joinpath("include", "sizes.json")
    if dst_json.exists():
        os.remove(dst_json)
    shutil.copy(src_json, dst_json)

    # ----- GENERATE HYDRO ADVANCE DATA PACKET
    logger = milhoja.BasicLogger(milhoja.LOG_LEVEL_BASIC)    
    tf_spec = milhoja.TaskFunction.from_milhoja_json(tf_json)

    milhoja.generate_data_item(
        tf_spec, DESTINATION, NO_OVERWRITE, milhoja_path, INDENT, logger
    )

    os.remove(dst_json)
