# create $HOME/feel directory
mkdir -p /home/$USER/feel
# start angiotk docker image
docker run --pid=host --ipc=host -it -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix -v /home/$USER/.Xauthority:/home/feelpp/.Xauthority -v /home/$USER/tp/:/home/feelpp/tp/ -v /home/$USER/feel:/feel feelpp/angiotk:master-ubuntu-16.10
