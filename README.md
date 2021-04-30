# quantum_zeno_paradox

C++ program for simulating measurement of spin on X axis that is precessing due to field along Z axis. 
Part of it is class that demonstrates quantum Zeno paradox by reducing interval between each measurement.

For making plots, matplotlib_cpp is used. 
So in order for that to work, you'll need to compile with:
```
g++ -o QuantumZenoParadox QuantumZenoParadox.cpp -std=c++14 -I/usr/include/python3.8 -lpython3.8
```
where you need to set right python version that you have on your computer. 
For more details about that library check https://github.com/lava/matplotlib-cpp.

The code was written as part of the project for uni.