This folder contains the C++ and MATLAB source code for the mex files. The existing mex files are built for x64 systems. If you're running an x86 system, you have two choices:
* Compile the C++ source using your own compiler, then place the resulting mex file in the "MEX Files" folder.
* Copy the MATLAB source code into the "Utils" folder and add that folder to the path.