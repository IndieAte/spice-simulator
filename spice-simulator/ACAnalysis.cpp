#define _USE_MATH_DEFINES

#include "ACAnalysis.h"

using namespace Eigen;

// Function prototypes for functions that don't need to be externally visible
Vector3d voltageVectorToPolar(int outputNode, VectorXcd voltageVector, double frequency);
VectorXcd solveAtFrequency(std::vector<Component*> comps, std::vector<int> cSIndexes,
	std::vector<int> vSIndexes, std::vector<int> nSIndexes, int nNodes, double angFreq);
void nonSourceHandler(Component* comp, MatrixXcd& gMat, double angFreq);
void currentSourceHandler(Component* comp, VectorXcd& iVec);
void voltageSourceHandler(Component* comp, MatrixXcd& gMat, VectorXcd& iVec);

/* Function runACAnalysis
*		This function is responsible for running the AC analysis of the given circuit. It
*		should calculate the response of the circuit at frequencies spaced out logarithmically
*		over the specified range and return the results
*  
*  Inputs:
*		int outNode										- The ID of the node to be treated as the output
*		double startFreq							- The start frequency of the AC analysis
*		double stopFreq								- The stop frequency of the AC analysis
*		int pPD												- The points per decade of the AC analysis
*		std::vector<Component*> comps	- The vector of components in the circuit
*		int nNodes										- The number of nodes in the circuit, excluding ground
* 
*	 Outputs:
*		std::vector<Vector3d> output	- Vector of Eigen vectors, where each Eigen vector contains
*																	  the amplitude and phase of the voltage at the output
*																		node at a given frequency
*/
std::vector<Vector3d> runACAnalysis(int outNode, double startFreq, double stopFreq, int pPD,
	std::vector<Component*> comps, int nNodes) {

	// Check that the output node is actually a node in the netlist
	try {
		if (outNode > nNodes) {
			throw std::invalid_argument("Output node not in netlist");
		}
	} catch (std::invalid_argument& e) {
		std::cerr << e.what() << std::endl;
	}

	// Initialise vectors that will contain information about what kinds of components
	// are where in the component vector
	std::vector<int> vSIndexes, cSIndexes, nSIndexes;
	std::vector<int> vSTmp;

	// Loop over the component vector and populate the vectors with information
	// about where source and non-source components are
	for (int i = 0; i < comps.size(); i++) {
		Component* c = comps[i];

		if (typeid(*c) == typeid(ACCurrentSource) || typeid(*c) == typeid(DCCurrentSource)) {
			cSIndexes.push_back(i);
		// Voltage sources require a little more care - we want to handle DC or 0 value voltage
		// sources first, or otherwise there will be minor inaccuracies for voltage sources
		// with meaningful values at AC
		} else if (typeid(*c) == typeid(ACVoltageSource)) {
			std::vector<double> ppts = c->getProperties();
			if (ppts[0] != 0) {
				vSTmp.push_back(i);
			} else {
				vSIndexes.push_back(i);
			}
		} else if (typeid(*c) == typeid(DCVoltageSource)) {
			vSIndexes.push_back(i);
		} else {
			nSIndexes.push_back(i);
		}
	}

	// Add the meaningful value voltage sources onto the back of vSIndexes so they're
	// handled last
	for (int i = 0; i < vSTmp.size(); i++) {
		vSIndexes.push_back(vSTmp[i]);
	}

	// Initialise values for the logarithmic frequency sweep
	double currentFreq = startFreq;
	int n = 1;
	std::vector<Vector3d> output;

	while (currentFreq <= stopFreq) {
		// Convert frequency to angular frequency
		double currAngFreq = currentFreq * M_PI * 2;

		// Solve the circuit at the current frequency
		VectorXcd voltageVector = solveAtFrequency(comps, cSIndexes, vSIndexes, nSIndexes, nNodes, currAngFreq);
		output.push_back(voltageVectorToPolar(outNode, voltageVector, currentFreq));

		// Increment the frequency logarithmically to ensure there are pPD points per decade
		// and that the sweep runs over the correct frequencies
		double exponent = n / (double)pPD;
		currentFreq = pow(10, exponent) * startFreq;
		n++;
	}

	return output;
}

/* Function voltageVectorToPolar
*		This is a helper function to convert the outputs of solveAtFrequency
*		to polar form and package this information with the frequency at which the 
*		result was found
* 
*  Inputs:
*		int outNode					- The ID of the node to be treated as the output
*		VectorXcd voltVect  - The vector of node voltages in rectangular form
*		double freq					- The frequency at which this result was calculated
* 
*	 Outputs:
*		Vector3d output			- The node voltage of outNode, expressed in polar form
*/ 
Vector3d voltageVectorToPolar(int outNode, VectorXcd voltVect, double freq) {
	Vector3d output;
	int oNIndex = outNode - 1;

	double amplitude = abs(voltVect(oNIndex));
	double phase = arg(voltVect(oNIndex));

	output << amplitude, phase, freq;

	return output;
}

void convertToSmallSignal(std::vector<Component*>& comps) {
	std::vector<int> nlcIndexes;

	for (int i = 0; i < comps.size(); i++) {
		Component* c = comps[i];

		if (typeid(*c) == typeid(Diode)) nlcIndexes.push_back(i);
	}

	if (nlcIndexes.size() == 0) {
		return;
	}
}

/* Function solveAtFrequency
*		This function solves the node voltages of the circuit at the given angular frequency
*		and returns them in rectangular form
* 
*  Inputs:
* 	std::vector<Component*> comps - Vector of pointers to all components in the circuit
* 	std::vector<int> cSIndexes    - Vector of the indexes of current sources in comps
*   std::vector<int> vSIndexes    - Vector of the indexes of voltage sources in comps
* 	std::vector<int> nSIndexes		- Vector of the indexes of non source components in comps
* 	int nNodes										- The number of nodes in the circuit, excluding ground
* 	double angFreq								- The angular frequency at which to solve the circuit
* 
*  Outputs:
* 	VectorXcd soln								- A vector of complex numbers representing the voltages at 
*																	  each node
*/
VectorXcd solveAtFrequency(std::vector<Component*> comps, std::vector<int> cSIndexes, 
	std::vector<int> vSIndexes, std::vector<int> nSIndexes, int nNodes, double angFreq) {

	// Initialise the conductance matrix and the current vector
	MatrixXcd gMat = MatrixXcd::Zero(nNodes, nNodes);
	VectorXcd iVec = VectorXcd::Zero(nNodes);

	// Update the conductance matrix with the conductances of all non-source
	// components
	for (int i = 0; i < nSIndexes.size(); i++) {
		int j = nSIndexes[i];
		nonSourceHandler(comps[j], gMat, angFreq);
	}

	// Update the current vector with the currents from current sources
	for (int i = 0; i < cSIndexes.size(); i++) {
		int j = cSIndexes[i];
		currentSourceHandler(comps[j], iVec);
	}

	// Update the conductance matrix and current vector with the effects of
	// voltage sources
	for (int i = 0; i < vSIndexes.size(); i++) {
		int j = vSIndexes[i];
		voltageSourceHandler(comps[j], gMat, iVec);
	}

	// Initialise the solution vector
	VectorXcd soln;
	// Solve for the nodal voltages using an Eigen solver
	soln = gMat.colPivHouseholderQr().solve(iVec);
	return soln;
}

/* Function nonSourceHandler
*		This function updates the conductance matrix with the effects of
*		non-source components (ie R, C, L, D, etc) at the given angular frequency
* 
*  Inputs:
*		Component* comp	- The component being considered
*		MatrixXcd& gMat	- The conduction matrix as it currently is
*		double angFreq	- The angular frequency being considered
* 
*  Outputs (By Reference):
*		MatrixXcd& gMat	- The conduction matrix with the update
*/
void nonSourceHandler(Component* comp, MatrixXcd& gMat, double angFreq) {
	// Get a vector of the nodes connected to comp and thus see how many terminals
	// comp has
	std::vector<int> nodes = comp->getNodes();
	int nNodes = nodes.size();

	// Two terminal components are the simplest and most common case
	if (nNodes == 2) {
		// Shorter names for the nodes connected to comp
		int n0 = nodes[0];
		int n1 = nodes[1];

		// If comp connects a node to itself, then it will have no effect
		// as it is effectively short circuited
		if (n0 == n1) return;

		// Shorter names for the indexes within the conductance matrix
		// of the nodes connected to comp
		int n0i = n0 - 1;
		int n1i = n1 - 1;

		// Calculate the conductance of comp between it's terminals
		std::complex<double> g = comp->getConductance(n0, n1, angFreq);

		// Update the conduction matrix appropriately depending on whether either
		// node is ground
		if (n0 != 0 && n1 != 0) {
			gMat(n0i, n1i) -= g;
			gMat(n1i, n0i) -= g;
		}

		if (n0 != 0) gMat(n0i, n0i) += g;
		if (n1 != 0) gMat(n1i, n1i) += g;
	}
}

/* Function currentSourceHandler
*		This function updates the current vector with the effects of current
*		sources
* 
*	Inputs:
*		Component* comp	- The current source component
*		VectorXcd& iVec	- The current vector as it currently is
* 
* Outputs (By Reference):
*		VectorXcd& iVec	- The current vector after having been updated
*/
void currentSourceHandler(Component* comp, VectorXcd& iVec) {
	// Check if the current source is AC, as DC current sources will
	// have no effect on the AC analysis
	if (typeid(*comp) == typeid(ACCurrentSource)) {
		// Get the nodes connected to the source and name them as the input
		// and output for clarity
		std::vector<int> nodes = comp->getNodes();
		int nIn = nodes[0] - 1;
		int nOut = nodes[1] - 1;

		// If the current source is short circuited it will have no effect
		if (nIn == nOut) return;

		// Get the current of the source as a phasor
		std::vector<double> ppts = comp->getProperties();
		double ampl = ppts[0];
		double phase = ppts[1];
		std::complex<double> current = std::polar(ampl, phase);

		// Update the current vector appropriately depending on
		// whether either node is ground
		if (nIn != -1 && nOut != -1) {
			iVec(nIn) -= current;
			iVec(nOut) += current;
		} else if (nIn == -1) {
			iVec(nOut) += current;
		} else {
			iVec(nIn) -= current;
		}
	}
}

/* Function voltageSourceHandler
*		This function updates the conduction matrix and the current vector
*		with the effects of voltage sources
* 
*  Inputs:
*		Component* comp - The voltage source component
*		MatrixXcd& gMat - The conduction matrix as it currently is
*		VectorXcd& iVec - The current vector as it currently is
* 
*  Outputs (By Reference):
*		MatrixXcd& gMat	- The conduction matrix after having been updated
*		VectorXcd& iVec - The current vector after having been updated
*/
void voltageSourceHandler(Component* comp, MatrixXcd& gMat, VectorXcd& iVec) {
	// Get the nodes connected to the source and name them as positive
	// and negative for clarity
	std::vector<int> nodes = comp->getNodes();
	int nPos = nodes[0] - 1;
	int nNeg = nodes[1] - 1;

	// Check if the voltage source is short circuited and if it is, throw an
	// error as the node it's connected to will be undefined
	try {
		if (nPos == nNeg) throw std::invalid_argument("Voltage source shorted");
	} catch (std::invalid_argument& e) {
		std::cerr << e.what() << std::endl;
	}

	// Set the voltage to 0 by default, meaning no work has to be done for
	// DC sources
	std::complex<double> voltage = 0;

	// If we've got an AC source we need to find the voltage from it's amplitude 
	// and phase
	if (typeid(*comp) == typeid(ACVoltageSource)) {
		std::vector<double> ppts = comp->getProperties();
		double ampl = ppts[0];
		double phase = ppts[1];
		voltage = std::polar(ampl, phase);
	}

	// Handle the update appropriately depending on if either node is ground
	if (nPos != -1 && nNeg != -1) {
		// If the voltage source is floating, we need to add a new unknown, the current
		// through it, which means adding an extra colum to the conduction matrix. We
		// also need to add an extra equation to the conduction matrix representing the
		// voltage between nodes

		// Add an extra row and column to the conduction matrix
		int cols = gMat.cols() + 1;
		int rows = gMat.rows() + 1;
		gMat.conservativeResize(rows, cols);
		// Add an extra row to the current vector
		iVec.conservativeResize(rows);

		// Set the new rows to 0
		rows--;
		cols--;
		gMat.row(rows).setZero();
		gMat.col(cols).setZero();

		// Fill the new row of the conduction matrix with the voltage equation
		gMat(rows, nPos) = 1;
		gMat(rows, nNeg) = -1;
		// Attach the unknown current to the correct nodes
		gMat(nPos, cols) = -1;
		gMat(nNeg, cols) = 1;
		// Add the other side of the voltage equation to the current vector
		iVec(rows) = voltage;

	} else if (nNeg == -1) {
		// If the negative node is ground, the positive node is set to the voltage
		gMat.row(nPos).setZero();
		gMat(nPos, nPos) = 1;
		iVec(nPos) = voltage;
	} else {
		// If the positive node is ground, the negative node is set to negative the voltage
		gMat.row(nNeg).setZero();
		gMat(nNeg, nNeg) = -1;
		iVec(nNeg) = voltage;
	}
}