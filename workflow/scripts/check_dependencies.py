import argparse
desc = "Checks PATH for versions of R and MEME indicated in config.yaml."

parser = argparse.ArgumentParser(
    description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
    "--MEME",
    type=str,
    required = True,
    help="Version number of installed version of MEME",
)
parser.add_argument(
    "--R",
    type=str,
    required = True,
    help="Version number of installed version of R",
)
args = parser.parse_args()
def check_dependency(dep, version):
    print_dep = dep.upper()
    if shutil.which(dep) != None:
        v = "--version"
        if dep == "meme":
            v = "-version"
       
        ro = subprocess.run([dep, v], capture_output = True)
        ver = version in str(ro.stdout) if ro.stdout != None else None
        #print(ver)
        if ver != None and version in str(ro.stdout):
            print(f"Successfully detected {print_dep} {version}.")
        elif  ver != None:
            print(f"{print_dep} {version}  not found in PATH. Please install {print_dep} {version} or change config to match installed version.")
            sys.exit()
        else:
            print(f"{print_dep} not found in PATH. Please install {print_dep} {version}.")
            sys.exit()
    else:
        print(f"{print_dep} not found in PATH. Please install {print_dep} {version}.")
        sys.exit()


#def main(args):
import shutil
import sys
import subprocess
from subprocess import CalledProcessError
deps = ['meme', 'R']
versions = dict(zip(deps, [args.MEME, args.R]))
for dep in deps:
    try:
        check_dependency(dep, versions[dep])
    except CalledProcessError:
        print(f"{print_dep} not found in PATH. Please install {print_dep} {version}.")
        sys.exit()
print('Dependency check completed.')