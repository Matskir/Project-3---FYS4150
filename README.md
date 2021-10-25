# Project-3---FYS4150

Project3_FYS4150_2.cpp is one of the main C++ files for this project, and includes code to compute the motions of particles in the Penning trap. As it is set up in this repository it computes the motions of two particles without coulomb forces, as this is the fastest by far. It is Project3_FYS4150_2.cpp that is used to create all data files used in the python files to create the necessary plots. The .cpp file should be compiled and run in the following way:

g++ Project3_FYS4150_2.cpp -o project3_fys4150.exe -O3 -larmadillo

./project3_fys4150.exe

task_10.cpp contains the C++ code needed to simulate 100 particles inside the trap subjected to a time-dependent electric field. As it is set up in this repository it does this in the frequency interval [0.2, 0.6] MHz with a step size of 0.02 MHz with amplitude f = 0.4 without taking coulomb forces into account for speed's sake. This can easily be changed. Running the file creates a text file containing the number of particles still inside the Penning trap and the corresponding frequency, and can be read by read_resonance.py to create a plot. The file should be run in the following way:

g++ task_10.cpp -o task_10.exe -O3 -larmadillo

./task_10.exe

basic_read_plot.py reads data files created by Project3_FYS4150.exe and makes neat plots needed for the Result section. All the data files needed to run the .py file are included in the repository, but could also be created manually using Project3_FYS4150_2.cpp. The file should be run in the following way:

python basic_read_plot.py

read_calc_errors.py reads files from Project3_FYS4150_2.cpp and computes the relative errors from them, making nice plots. As some of the data files needed for this are relatively large they are not included in the repository, but could be created manually using the .cpp file if one is patient. The file should be run as:

python read_calc_errors.py

read_resonance.py is the last code file and reads data files from task_10.cpp and make plots of the fraction of particles still remaining within the Penning trap. As this took extremely long to run, the needed data files are here included. The file should be run as:

python read_resonance.py

