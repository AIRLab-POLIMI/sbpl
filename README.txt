This folder contains the sbpl library. To use it follow the instructions from root directory of sbpl:

mkdir build
cd build
cmake ..
make
sudo make install

Now the library is installed. To test it there are some environment configuration files in env_examples and in myenvfiles directories and there are some motion primitives files in matlab/mprim and in myprimitives directories.

Depending on which executable is launched there are different parameters to specify, but launching them with -h parameters an help appear.
Do not launch performance_test executable, it is only to do various sequential tests on planners and environments.
