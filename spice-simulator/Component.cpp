#include "Component.h"

// ==================== COMPONENT (BASE CLASS) ====================

std::vector<int> Component::getNodes() {
	std::vector<int> nodes;

	return nodes;
}

double Component::getConductance(int p_node1, int p_node2) {
	return 0;
}

std::vector<double> Component::getProperties() {
	std::vector<double> properties;

	return properties;
}

// ======================== CURRENT SOURCE ========================

std::vector<int> CurrentSource::getNodes() {
	std::vector<int> nodes;
	
	// CurrentSource should return the input node first, output node
	// second to make creating the conduction matrix easier
	nodes.push_back(nodeIn);
	nodes.push_back(nodeOut);

	return nodes;
}

double CurrentSource::getConductance(int p_node1, int p_node2) {
	return 0;
}

std::vector<double> CurrentSource::getProperties() {
	std::vector<double> properties;

	properties.push_back(current);

	return properties;
}

// =========================== RESISTOR ===========================

std::vector<int> Resistor::getNodes() {
	std::vector<int> nodes;

	nodes.push_back(node1);
	nodes.push_back(node2);

	return nodes;
}

double Resistor::getConductance(int p_node1, int p_node2) {
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