import shutil
import sys

versions = dict(zip(sys.argv[1].split(",") , sys.argv[2].split(",")))
def check_dependency(dep, version):
    if shutil.which(dep) != None:
        if not shutil.which(dep).contains(version):
            print(f"{dep} {version} not found in PATH. Please install this version or change config version to installed version.")
            sys.exit()
    else:
        print(f"{dep} not found in PATH. Please install {dep} {version}.")
        sys.exit()
        
for dep in versions.keys():
    check_dependency(dep, versions[dep])