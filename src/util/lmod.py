import os
import subprocess
import warnings


class LmodError(Exception):
    pass


def load_module(module_name):
    lmod = os.environ.get("LMOD_CMD")
    if lmod is None:
        raise LmodError('Environment variable "LMOD_CMD" not set. Is lmod available?')

    cmd = [lmod, "python", "load", module_name]
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, stderr = process.communicate()

    if process.returncode:
        raise LmodError(stderr.decode("utf-8"))
    if stderr:
        warnings.warn(stderr.decode("utf-8"))
    exec(stdout)
