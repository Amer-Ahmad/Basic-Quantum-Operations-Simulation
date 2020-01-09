# Basic-Quantum-Operations-Simulation
this project is a simple simulation for handful of basic operations that could be applied on a quantum register or a linear operator.

### QuantumRegister
operations that are in the code: initialization, norm, innerProduct, TensorProduct, measurement, permuting a QuantumState.

### LinearOperator
operations that are in the code: initialization, initialization as outer product of two quantum states, addition of two linear operators and multiplying with a scalar, transform (applying a linear operator on a quantum state), partialTransform (applying a linear operator on certain indices of a quantum state).

### Tester
I initialized a bunch of random quantum states to test the methods on (keeping in mind that the norm of a quantum state is equal to 1).

this project was done as a part of the assignments for my intro to quantum computing class. Moreover, the educational value of it is to help students visualize how a quantum register function in a simplified manner, especially how a measurement is done on a quantum register.

#### P.S: the project was done on python 2.7, incase you are using another version you might face a problem in running the program mainly due to the random library. Thus, do the neccessary adjustments before running the program.
