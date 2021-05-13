#include "Component.h"

using namespace std::complex_literals;

// ==================== COMPONENT (BASE CLASS) ====================

std::vector<int> Component::getNodes() {
	std::vector<int> nodes;

	return nodes;
}

std::complex<double> Component::getConductance(int p_node1, int p_node2, double p_angularFrequency) {
	return 0;
}

std::vector<double> Component::getProperties() {
	std::vector<double> properties;

	return properties;
}

// ======================= AC CURRENT SOURCE ======================

std::vector<int> ACCurrentSource::getNodes() {
	std::vector<int> nodes;

	// ACCurrentSource should return the input node first, output node
	// second to make creating the conduction matrix easier
	nodes.push_back(nodeIn);
	nodes.push_back(nodeOut);

	return nodes;
}

std::complex<double> ACCurrentSource::getConductance(int p_node1, int p_node2, double p_angularFrequency) {
	return 0;
}

std::vector<double> ACCurrentSource::getProperties() {
	std::vector<double> properties;

	properties.push_back(amplitude);
	properties.push_back(phase);

	return properties;
}

// ======================= DC CURRENT SOURCE ======================

std::vector<int> DCCurrentSource::getNodes() {
	std::vector<int> nodes;
	
	// DCCurrentSource should return the input node first, output node
	// second to make creating the conduction matrix easier
	nodes.push_back(nodeIn);
	nodes.push_back(nodeOut);

	return nodes;
}

std::complex<double> DCCurrentSource::getConductance(int p_node1, int p_node2, double p_angularFrequency) {
	return 0;
}

std::vector<double> DCCurrentSource::getProperties() {
	std::vector<double> properties;

	properties.push_back(current);

	return properties;
}

// ======================= AC VOLTAGE SOURCE ======================

std::vector<int> ACVoltageSource::getNodes() {
	std::vector<int> nodes;

	// ACVoltageSource should return the positive node first, then
	// the negative node
	nodes.push_back(nodePlus);
	nodes.push_back(nodeMinus);

	return nodes;
}

std::complex<double> ACVoltageSource::getConductance(int node1, int node2, double p_angularFrequency) {
	return 0;
}

std::vector<double> ACVoltageSource::getProperties() {
	std::vector<double> properties;

	properties.push_back(amplitude);
	properties.push_back(phase);

	return properties;
}

// ======================= DC VOLTAGE SOURCE ======================

std::vector<int> DCVoltageSource::getNodes() {
	std::vector<int> nodes;

	// DCVoltageSource should return the positive node first, then
	// the negative node
	nodes.push_back(nodePlus);
	nodes.push_back(nodeMinus);

	return nodes;
}

std::complex<double> DCVoltageSource::getConductance(int node1, int node2, double p_angularFrequency) {
	return 0;
}

std::vector<double> DCVoltageSource::getProperties() {
	std::vector<double> properties;

	properties.push_back(voltage);

	return properties;
}

// =========================== RESISTOR ===========================

std::vector<int> Resistor::getNodes() {
	std::vector<int> nodes;

	nodes.push_back(node1);
	nodes.push_back(node2);

	return nodes;
}

std::complex<double> Resistor::getConductance(int p_node1, int p_node2, double p_angularFrequency) {
	if (p_node1 != node1 && p_node1 != node2
		|| p_node2 != node1 && p_node2 != node2) {
		
		return 0;
	} else {
		return 1 / resistance;
	}
}

std::vector<double> Resistor::getProperties() {
	std::vector<double> properties;

	properties.push_back(resistance);

	return properties;
}

// =========================== CAPACITOR ===========================

std::vector<int> Capacitor::getNodes() {
	std::vector<int> nodes;

	nodes.push_back(node1);
	nodes.push_back(node2);

	return nodes;
}

std::complex<double> Capacitor::getConductance(int p_node1, int p_node2, double p_angularFrequency) {
	if (p_node1 != node1 && p_node1 != node2
		|| p_node2 != node1 && p_node2 != node2) {
		
		return 0;
	} else {
		return capacitance * p_angularFrequency * 1i;
	}
}

std::vector<double> Capacitor::getProperties() {
	std::vector<double> properties;

	properties.push_back(capacitance);

	return properties;
}

// =========================== INDUCTOR ===========================

std::vector<int> Inductor::getNodes() {
	std::vector<int> nodes;

	nodes.push_back(node1);
	nodes.push_back(node2);

	return nodes;
}

std::complex<double> Inductor::getConductance(int p_node1, int p_node2, double p_angularFrequency) {
	if (p_node1 != node1 && p_node1 != node2
		|| p_node2 != node1 && p_node2 != node2) {
		
		return 0;
	} else {
		return 1.0 / (inductance * p_angularFrequency * 1i);
	}
}

std::vector<double> Inductor::getProperties() {
	std::vector<double> properties;

	properties.push_back(inductance);

	return properties;
}

// =========================== DIODE ==============================

std::vector<int> Diode::getNodes() {
	std::vector<int> nodes;

	nodes.push_back(nodeAnode);
	nodes.push_back(nodeCathode);

	return nodes;
}

std::complex<double> Diode::getConductance(int p_node1, int p_node2, double p_angularFrequency) {
	if (p_node1 != nodeAnode && p_node1 != nodeCathode
		|| p_node2 != nodeAnode && p_node2 != nodeCathode) {
		
		return 0;
	} else {
		return 0; // NEED CONDUCTANCE HERE
	}
}

std::vector<double> Diode::getProperties() {
	std::vector<double> properties;

	properties.push_back(modelName);

	return properties;
}

// =========================== BJT ================================

std::vector<int> BJT::getNodes() {
	std::vector<int> nodes;

	nodes.push_back(nodeCollector);
	nodes.push_back(nodeBase);
	nodes.push_back(nodeEmitter);

	return nodes;
}

std::complex<double> Diode::getConductance(int p_node1, int p_node2, double p_angularFrequency) {
	if (!((std::count(getNodes().begin(), getNodes().end(), p_node1) || 
	std::count(getNodes().begin(), getNodes().end(), p_node2)))) {
		return 0;
	} else {
		return 0; // NEED CONDUCTANCE HERE
	}
}

std::vector<double> Diode::getProperties() {
	std::vector<double> properties;

	properties.push_back(modelName);

	return properties;
}