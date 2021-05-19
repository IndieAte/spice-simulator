#include "DCOpPoint.h"

using namespace Eigen;

VectorXd iterate(std::vector<Component*> comps, std::vector<int> cSIndexes, std::vector<int> vSIndexes,
	std::vector<int> lCIndexes, std::vector<int> nlCIndexes, int nNodes);
void linearComponentHandler(Component* comp, MatrixXd& gMat);
void nonlinearComponentHandler(Component* comp, MatrixXd& gMat, VectorXd& iVec);
void DCcurrentSourceHandler(Component* comp, VectorXd& iVec);
void DCvoltageSourceHandler(Component* comp, MatrixXd& gMat, VectorXd& iVec);
void updateNonlinearComponent(Component* comp, VectorXd vVec);

/* Function runDCOpPoint
*		This function is responsible for running a DC operating point analysis of
*		the circuit specified by the vector of components comps and returning the
*		node voltages as an Eigen vector
* 
*	 Inputs:
*		std::vector<Component*> comps - Vector of component objects specifiying the circuit
*		int nNodes                    - Number of nodes in the circuit
* 
*	 Outputs:
*		VectorXd currSoln             - The node voltages once the simulation has converged
*/
VectorXd runDCOpPoint(std::vector<Component*> comps, int nNodes) {
	// Initialise vectors containing indexes of various component types to
	// give each iteration quick access
	std::vector<int> cSIndexes, vSIndexes, lCIndexes, nlCIndexes;
	std::vector<int> vSTmp;

	// Iterate over comps to populate the index vectors
	for (int i = 0; i < comps.size(); i++) {
		Component* c = comps[i];

		if (typeid(*c) == typeid(ACCurrentSource) || typeid(*c) == typeid(DCCurrentSource)) {
			cSIndexes.push_back(i);
		// Note that we need to take care with voltage sources, as we want to process 0 valued voltage
		// sources first
		} else if (typeid(*c) == typeid(DCVoltageSource)) {
			std::vector<double> ppts = c->getProperties();
			if (ppts[0] != 0) {
				vSTmp.push_back(i);
			} else {
				vSIndexes.push_back(i);
			}
		// Also note that we treat inductors as 0 value voltage sources for DC operating point analysis
		} else if (typeid(*c) == typeid(Inductor) || typeid(*c) == typeid(ACVoltageSource)) {
			vSIndexes.push_back(i);
		} else if (typeid(*c) == typeid(Diode)) {
			nlCIndexes.push_back(i);
		} else {
			lCIndexes.push_back(i);
		}
	}

	// Make sure that meaningfully valued voltage sources are handled last so they overwrite 0 valued
	// sources, preventing minor errors in voltages
	for (int i = 0; i < vSTmp.size(); i++) {
		vSIndexes.push_back(vSTmp[i]);
	}

	// Create two voltage vectors and populate them with initial values
	VectorXd prevSoln, currSoln;
	prevSoln = iterate(comps, cSIndexes, vSIndexes, lCIndexes, nlCIndexes, nNodes);
	currSoln = iterate(comps, cSIndexes, vSIndexes, lCIndexes, nlCIndexes, nNodes);

	int n = 0;

	// Iterate the analysis until the current solution is approximately the
	// previous solution (ie convergence) or until an iteration cap is reached
	while (!currSoln.isApprox(prevSoln) && n < 40) {
		n++;
		prevSoln = currSoln;
		currSoln = iterate(comps, cSIndexes, vSIndexes, lCIndexes, nlCIndexes, nNodes);
	}

	// If we reach the iteration cap, alert the user so they know results may be
	// inaccurate
	if (n == 40) {
		std::cerr << "DC operating point iteration limit exceeded" << std::endl;
	}

	// Truncate the voltage vector to remove any unknown currents that were added to
	// handle floating voltage sources
	currSoln.conservativeResize(nNodes);

	return currSoln;
}

/* Function iterate
*		This function runs one iteration of the DC operating point analysis, which involves forming
*		the conductance matrix and current vector for the current iteration and solving it, which
*		performs a Newton-Raphson iteration so long as non-linear components are handled correctly
* 
*  Inputs:
*		std::vector<Component*> comps - Vector of components in the circuit
*		std::vector<int> cSIndexes    - Indexes of current sources within comps
*		std::vector<int> vSIndexes    - Indexes of voltage sources and inductors within comps
*		std::vector<int> lCIndexes    - Indexes of linear components within comps
*		std::vector<int> nlCIndexes   - Indexes of nonlinear components within comps
*		int nNodes                    - The number of nodes in the circuit, excluding ground
* 
*  Outputs:
*		VectorXd vVec                 - Vector of node voltages after the iteration
*/
VectorXd iterate(std::vector<Component*> comps, std::vector<int> cSIndexes, std::vector<int> vSIndexes,
	std::vector<int> lCIndexes, std::vector<int> nlCIndexes, int nNodes) {
	// Initialise the conductance matrix and current vector
	MatrixXd gMat = MatrixXd::Zero(nNodes, nNodes);
	VectorXd iVec = VectorXd::Zero(nNodes);

	// Loops over linear components and call the handler for them
	for (int i = 0; i < lCIndexes.size(); i++) {
		int j = lCIndexes[i];
		linearComponentHandler(comps[j], gMat);
	}

	// Loops over nonlinear components and call the handler for them
	for (int i = 0; i < nlCIndexes.size(); i++) {
		int j = nlCIndexes[i];
		nonlinearComponentHandler(comps[j], gMat, iVec);
	}

	// Loops over current sources and call the handler for them
	for (int i = 0; i < cSIndexes.size(); i++) {
		int j = cSIndexes[i];
		DCcurrentSourceHandler(comps[j], iVec);
	}

	// Loops over voltage sources (and inductors) and call the handler for them
	for (int i = 0; i < vSIndexes.size(); i++) {
		int j = vSIndexes[i];
		DCvoltageSourceHandler(comps[j], gMat, iVec);
	}

	// Solve for the nodal voltages
	VectorXd vVec = gMat.colPivHouseholderQr().solve(iVec);

	// Loops over nonlinear components and updates their properties based on the new nodal voltages
	for (int i = 0; i < nlCIndexes.size(); i++) {
		int j = nlCIndexes[i];
		updateNonlinearComponent(comps[j], vVec);
	}

	return vVec;
}

/* Function linearComponentHandler
*		This function updates the conductance matrix with the effects of linear components
* 
*  Inputs:
*		Component* comp - The linear component being handled
*		MatrixXd& gMat  - The conductance matrix before being updated
* 
*  Outputs (By Reference):
*		MatrixXd& gMat  - The conductance matrix after being updated
*/
void linearComponentHandler(Component* comp, MatrixXd& gMat) {
	// Get the nodes connected to the component
	std::vector<int> nodes = comp->getNodes();

	int n0 = nodes[0];
	int n1 = nodes[1];

	// Handle resistors
	if (typeid(*comp) == typeid(Resistor)) {
		// Get indexes into gMat for the two nodes
		int n0i = n0 - 1;
		int n1i = n1 - 1;

		// Get the conductance of the resistor
		 double g = std::real(comp->getConductance(n0, n1, 0));

		 // Update the conductance matrix accordingly to whether any of the connected
		 // nodes are ground
		if (n0 != 0 && n1 != 0) {
			gMat(n0i, n1i) -= g;
			gMat(n1i, n0i) -= g;
		}

		if (n0 != 0) gMat(n0i, n0i) += g;
		if (n1 != 0) gMat(n1i, n1i) += g;
	}
}

/* Function nonlinearComponentHandler
*		This function updates the conductance matrix and current vector with the effects of nonlinear components
* 
*  Inputs:
*		Component* comp - The nonlinear component being considered
*		MatrixXd& gMat  - The conductance matrix before being updated
*		VectorXd& iVec  - The current vector before being updated
* 
*  Outputs (By Reference):
*		MatrixXd& gMat  - The updated conductance matrix
*		VectorXd& iVec  - The updated current vector
*/
void nonlinearComponentHandler(Component* comp, MatrixXd& gMat, VectorXd& iVec) {
	// Get the nodes connected to the component
	std::vector<int> nodes = comp->getNodes();

	// Handle diodes
	if (typeid(*comp) == typeid(Diode)) {
		// Name the connected nodes for legibility
		int nAnode = nodes[0];
		int nCathode = nodes[1];

		// Get the indexes into gMat and iVec of the nodes
		int nAi = nAnode - 1;
		int nCi = nCathode - 1;

		// Get the companion model conductance and current of the diode in its current
		// state
		double g = std::real(comp->getConductance(nAnode, nCathode, 0));
		std::vector<double> ppts = comp->getProperties();
		double I = ppts[0];
		
		// Update gMat and iVec according to whether either node connected to the diode is ground
		if (nAi != -1 && nCi != -1) {
			gMat(nAi, nCi) -= g;
			gMat(nCi, nAi) -= g;
		}

		if (nAi != -1) { 
			gMat(nAi, nAi) += g; 
			iVec(nAi) -= I;
		}

		if (nCi != -1) { 
			gMat(nCi, nCi) += g;
			iVec(nCi) += I;
		}
	}
}

/* Function DCcurrentSourceHandler
*		This function handles updating the current vector with the effects of current sources (note that the
*		DC doesn't refer to DC current sources and just exists to disambiguate from a similar function in
*		ACAnalysis.cpp)
* 
*  Inputs:
*		Component* comp - The component being conidered
*		VectorXd& iVec  - The current vector before being updated
* 
*  Outputs (By Reference):
*		VectorXd& iVec  - The updated current vector
*/
void DCcurrentSourceHandler(Component* comp, VectorXd& iVec) {
	// We only need to make changes for DC current sources in DC operating point analysis
	if (typeid(*comp) == typeid(DCCurrentSource)) {
		// Get the nodes connected to the source and give them names
		std::vector<int> nodes = comp->getNodes();
		int nIn = nodes[0] - 1;
		int nOut = nodes[1] - 1;

		// Get the current of the source
		std::vector<double> ppts = comp->getProperties();
		double current = ppts[0];

		// Update the current vector depending on whether any nodes are ground
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

/* Function DCvoltageSourceHandler
*		This function handles updating the current vector and conductance matrix with the effects of 
*   voltage sources (note that the DC doesn't refer to DC voltage sources and just exists to 
*   disambiguate from a similar function in ACAnalysis.cpp)
*  
*  Inputs:
*		Component* comp - The component being considered
*		MatrixXd& gMat  - The conduction matrix before being updated
*		VectorXd& iVec  - The current vector before being updated
* 
*  Outputs (By Reference):
*		MatrixXd& gMat  - The updated conductance matrix
*		VectorXd& iVec  - The updated current vector
*/
void DCvoltageSourceHandler(Component* comp, MatrixXd& gMat, VectorXd& iVec) {
	// Get and name the nodes connected to the voltage source
	std::vector<int> nodes = comp->getNodes();
	int nPos = nodes[0] - 1;
	int nNeg = nodes[1] - 1;

	// By default assume the DC voltage is zero
	double voltage = 0;

	// Update the DC voltage if comp is a DC voltage source
	if (typeid(*comp) == typeid(DCVoltageSource)) {
		std::vector<double> ppts = comp->getProperties();
		voltage = ppts[0];
	}

	// Update gMat and iVec appropriately depending on whether either node
	// is ground
	if (nPos != -1 && nNeg != -1) {
		// Floating sources need special treatment, and need an extra unknown in the voltage vector
		// representing the current through the source and another equation describing the source.
		// This is implemented by adding an extra row and column to the conductance matrix and an extra
		// entry to the current vector
		int cols = gMat.cols() + 1;
		int rows = gMat.rows() + 1;
		gMat.conservativeResize(rows, cols);

		iVec.conservativeResize(rows);

		rows--;
		cols--;
		// If we don't call setZero, positions we don't explicitly set may contain garbage
		gMat.row(rows).setZero();
		gMat.col(cols).setZero();

		// Add the equation describing the source
		gMat(rows, nPos) = 1;
		gMat(rows, nNeg) = -1;

		// Add the unknown current through the source
		gMat(nPos, cols) = -1;
		gMat(nNeg, cols) = 1;

		iVec(rows) = voltage;
	} else if (nNeg == -1) {
		gMat.row(nPos).setZero();
		gMat(nPos, nPos) = 1;
		iVec(nPos) = voltage;
	} else {
		gMat.row(nNeg).setZero();
		gMat(nNeg, nNeg) = -1;
		iVec(nNeg) = voltage;
	}
}

/* Function updateNonlinearComponent
*		Updates the properties of a nonlinear component with the effect of the last
*		Newton-Raphson iteration
* 
*  Inputs:
*		Component* comp - The component to be updated
*		VectorXd vVec   - The results of the last iteration
*/
void updateNonlinearComponent(Component* comp, VectorXd vVec) {
	// Get the nodes connected to the component
	std::vector<int> nodes = comp->getNodes();

	// Update diodes
	if (typeid(*comp) == typeid(Diode)) {
		// Get the index into vVec of the nodes
		int nAi = nodes[0] - 1;
		int nCi = nodes[1] - 1;

		double Vd;

		// Calculate the voltage across the diode accordingly based on whether either node
		// is ground
		if (nAi != -1 && nCi != -1) Vd = vVec(nAi) - vVec(nCi);
		else if (nAi == -1) Vd = -vVec(nCi);
		else Vd = vVec(nAi);

		std::vector<double> ppts;
		ppts.push_back(Vd);

		// Update the voltage across the diode
		comp->setProperties(ppts);
	}
}