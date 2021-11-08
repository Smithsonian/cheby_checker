# MJP: 2021-06-30
# This routine functions as a hack/quick way to execute the various
# commands needed to create an image and fire-up a container of that image
#
# DOCKER_BUILDKIT=0
# DOCKER_BUILDKIT=1

import os, sys

container_name  = 'chb_cont'
image_name      = 'chb_imag'

def kill():
    # commands to remove any extant running versions of the mpcremote image and/or container
    os.system(f"docker stop {container_name}")
    os.system(f"docker rm {container_name}")
    os.system(f"docker image rm {image_name}")

def build():
    ''' commands needed to create an image and fire-up a container of that image'''
    
    # commands to remove any extant running versions of the mpcremote image and/or container
    kill()
    
    # build the image and name it "{image_name}" (NB --no-cache will force a rebuild)
    os.system(f"docker build . -t {image_name}")# --no-cache")
    
    # get the path of the parent directory (i.e. for the entire repo)
    repo_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) ; print(f"repo_dir={repo_dir}")

    # run the image in a container and name the container as "{container_name}"
    # NB: Mounting local version of repo during development
    command = " ".join([
        f"docker run -d",                           # run the container using some options ...
        f"-v {repo_dir}:/cheby_checker",            # mount local version of repo
        f"--name {container_name} {image_name}"     # name the container
        ])
    print(f"running command... {command}")
    os.system(command)

    # execute interactively (so you get into the command-line)
    os.system(f"docker exec -it {container_name} bash")

if __name__ == '__main__':
    
    # If supply 'k' option, simply kills (doesn't build),
    if len(sys.argv) > 1 and sys.argv[1][0].lower() =='k':
        kill()
    # As default, complete rebuild
    else:
        build()
