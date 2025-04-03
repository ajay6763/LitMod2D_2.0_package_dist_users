#!/bin/bash
xhost +
sudo docker run --rm -ti -v /home/ajay/projects/LitMod2D/LitMod2d_projects:/home/work  --net=host -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix --env QT_X11_NO_MITSHM=1 8d08a094c6e6
