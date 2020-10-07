1. Make the CX function able to work on any two qubits in the circuit.
a. thoughts: redo the inner for loop to be a helper function and then call this recursively.
The function would return an array of the matrices required to build the particular CX gate.
Lastly, this could then be multiplied out to create the desired CX gate.
2. From the remade CX function create the other core gates recquired to make the more interesting algorithms.
3. Make a true measurement function. This will determine the state of the given qubit. If the 
state was entangled prior to measurement, then the other qubits will be known at that time.
4. To assist with point 3, try to make a Tr over subsystem algorithm that is capable of removing
portions of the Hilbert space and generating the corresponding kets.