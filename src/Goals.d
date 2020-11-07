2. Make a true measurement function. This will determine the state of the given qubit. If the 
state was entangled prior to measurement, then the other qubits will be known at that time.
3. Try and make versions of the gates that operate directly on the involved qubits. Possibly copy the result of the control qubits to the targets.
If done properly, this will effectively remove the scaling issue the current version has as it will remove the need for the gate matrices to be large.
Should be able to remake the final statevector from the separate qubits in their particular superpositions. This would more closely emulate the actual devices.
4. Get QPE example working!!!!!!!!!!!!!!!!!!!