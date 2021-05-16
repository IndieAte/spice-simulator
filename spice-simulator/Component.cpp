#include "Component.h"

#define _VT 0.025851997

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

void Component::setProperties(std::vector<double> properties) {

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

void ACCurrentSource::setProperties(std::vector<double> properties) {

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

void DCCurrentSource::setProperties(std::vector<double> properties) {

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

void ACVoltageSource::setProperties(std::vector<double> properties) {

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

void DCVoltageSource::setProperties(std::vector<double> properties) {

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

void Resistor::setProperties(std::vector<double> properties) {

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

void Capacitor::setProperties(std::vector<double> properties) {

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

void Inductor::setProperties(std::vector<double> properties) {

}

// =========================== DIODE ==============================

Diode::Diode(std::string p_name, std::string p_modelName, int p_nodeAnode, int p_nodeCathode) : 
	Component{ p_name }, nodeAnode{ p_nodeAnode }, nodeCathode{ p_nodeCathode } {
	
	Vd = 0.7;

	if (p_modelName == "D") Is = pow(10, -12);
	else Is = pow(10, -12);

	Gd = (Is / _VT) * exp(Vd / _VT);
	Id = (Is * (exp(Vd / _VT) - 1)) - (Gd * Vd);
}

std::vector<int> Diode::getNodes() {
	std::vector<int> nodes;

	// Diode should return the anode first, then the cathode
	nodes.push_back(nodeAnode);
	nodes.push_back(nodeCathode);

	return nodes;
}

std::complex<double> Diode::getConductance(int p_node1, int p_node2, double p_angularFrequency) {
	if (p_node1 != nodeAnode && p_node1 != nodeCathode
		|| p_node2 != nodeAnode && p_node2 != nodeCathode) {
		
		return 0;
	} else {
		return Gd;
	}
}

std::vector<double> Diode::getProperties() {
	std::vector<double> properties;

	properties.push_back(Id);

	return properties;
}

void Diode::setProperties(std::vector<double> properties) {
	// For diode, the set properties vector should have Vd 
	// at index 0
	Vd = properties[0];

	Gd = (Is / _VT) * exp(Vd / _VT);
	Id = (Is * (exp(Vd / _VT) - 1)) - (Gd * Vd);
}

// =========================== BJT ================================

std::vector<int> BJT::getNodes() {
	std::vector<int> nodes;

	nodes.push_back(nodeCollector);
	nodes.push_back(nodeBase);
	nodes.push_back(nodeEmitter);

	return nodes;
}

std::complex<double> BJT::getConductance(int p_node1, int p_node2, double p_angularFrequency) {
	if (!((std::count(getNodes().begin(), getNodes().end(), p_node1) || 
	std::count(getNodes().begin(), getNodes().end(), p_node2)))) {
		return 0;
	} else {
		return 0; // NEED CONDUCTANCE HERE
	}
}

std::vector<double> BJT::getProperties() {
	std::vector<double> properties;

	properties.push_back(modelName);

	return properties;
}

void BJT::setProperties(std::vector<double> properties) {

}