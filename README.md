# quantum_zeno_paradox

C++ program for simulating measurement of spin on X axis that is precessing due to field along Z axis. 
By reducing the interval between each measurement we simulate quantum Zeno paradox.
Project also includes comparison of simualtion where we have collapse of wave function on measurement and simulation where we don't.

For making plots, matplotlib_cpp is used. 
So in order for that to work, you'll need that library and compile the project with:
```
g++ -o QuantumZenoParadox QuantumZenoParadox.cpp -std=c++14 -I/usr/include/python3.8 -lpython3.8
```
where you need to set right python version that you have on your computer.
For more details about that library check https://github.com/lava/matplotlib-cpp.
Or you can just delete or comment out that part of the project.