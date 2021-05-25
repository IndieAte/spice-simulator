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

// ============== VOLTAGE CONTROLLED CURRENT SOURCE ===============

std::vector<int> VoltageControlledCurrentSource::getNodes() {
	std::vector<int> nodes;
	
	// DCCurrentSource should return the input node first, output node
	// second to make creating the conduction matrix easier
	nodes.push_back(nodeIn);
	nodes.push_back(nodeOut);
	nodes.push_back(control_nodeIn);
	nodes.push_back(control_nodeOut);

	return nodes;
}

std::complex<double> VoltageControlledCurrentSource::getConductance(int p_node1, int p_node2, double p_angularFrequency) {
	return 0;
}

std::vector<double> VoltageControlledCurrentSource::getProperties() {
	std::vector<double> properties;

	properties.push_back(transconductance);

	return properties;
}

void VoltageControlledCurrentSource::setProperties(std::vector<double> properties) {

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

Diode::Diode(std::string p_name, int p_nodeAnode, int p_nodeCathode, Model* p_model) : 
	Component{ p_name }, nodeAnode{ p_nodeAnode }, nodeCathode{ p_nodeCathode }, model { p_model } {
	
	std::vector<double> model_values =  model->getDoubles();
	Is = model_values[0];

	Vd = 0.7;
	
	// try {
	// 	if (p_modelName == "D") Is = pow(10, -12);
	// 	else throw std::invalid_argument("Invalide diode model: " + p_modelName);
	// } catch (std::invalid_argument& e) {
	// 	std::cerr << e.what() << std::endl;
	// }

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
	properties.push_back(Is);

	return properties;
}

void Diode::setProperties(std::vector<double> properties) {
	// For diode, the set properties vector should have Vd 
	// at index 0
	Vd = properties[0];

	if (Vd > 1) Vd = 1;

	Gd = (Is / _VT) * exp(Vd / _VT);
	Id = (Is * (exp(Vd / _VT) - 1)) - (Gd * Vd);
}

// =========================== BJT ================================

BJT::BJT(std::string p_name, int p_nodeCollector, int p_nodeBase, int p_nodeEmitter, Model* p_model) :
	Component{ p_name }, nodeCollector{ p_nodeCollector }, nodeBase{ p_nodeBase },
	nodeEmitter{ p_nodeEmitter }, model { p_model } {
	
	std::vector<double> model_values = model->getDoubles();
	Vbe = 0.7;
	Vbc = 0.7;
	Is = model_values[0];
	bf = model_values[1];
	br = model_values[2];
	vaf = model_values[3];
	var = model_values[4];
	npn = model_values[5];

	// try {
	// 	if (modelName == "NPN") {
	// 		Vbe = 0.7;
	// 		Vbc = 0.7;
	// 		Is = pow(10, -12);
	// 		bf = 100;
	// 		br = 1;
	// 		npn = 1;
	// 	} else if (modelName == "PNP") {
	// 		Vbe = 0.7;
	// 		Vbc = 0.7;
	// 		Is = pow(10, -12);
	// 		bf = 100;
	// 		br = 1;
	// 		npn = 0;
	// 	} else {
	// 		throw std::invalid_argument("Invalid BJT model: " + modelName);
	// 	}
	// } catch (std::invalid_argument& e) {
	// 	std::cerr << e.what() << std::endl;
	// }

	double zeta = exp(Vbe / _VT);
	double xi = exp(Vbc / _VT);

	Gcc = (Is / _VT) * (1 + 1 / br) * xi;
	Gcb = (Is / _VT) * (zeta - (1 + 1 / br) * xi);
	Gce = -(Is / _VT) * zeta;

	Gbc = -(Is / (_VT * br)) * xi;
	Gbb = (Is / _VT) * ((zeta / bf) + (xi / br));
	Gbe = -(Is / (_VT * bf)) * zeta;

	Gec = -(Is / _VT) * xi;
	Geb = -(Is / _VT) * ((1 + 1 / bf) * zeta - xi);
	Gee = (Is / _VT) * (1 + 1 / bf) * zeta;

	if (npn == 1) {
		Ic = (Vbe * Is / _VT) * zeta - (Vbc * Is / _VT) * (1 + 1 / br) * xi - Is * (zeta - xi - (xi - 1) / br);
		Ib = (Vbe * Is / (bf * _VT)) * zeta + (Vbc * Is / (br * _VT)) * xi - Is * ((zeta - 1) / bf + (xi - 1) / br);
		Ie = (Vbc * Is / _VT) * xi - (Vbe * Is / _VT) * (1 + 1 / bf) * zeta + Is * (zeta - xi + (zeta - 1) / bf);
	} else {
		Ic = -((Vbe * Is / _VT) * zeta - (Vbc * Is / _VT) * (1 + 1 / br) * xi - Is * (zeta - xi - (xi - 1) / br));
		Ib = -((Vbe * Is / (bf * _VT)) * zeta + (Vbc * Is / (br * _VT)) * xi - Is * ((zeta - 1) / bf + (xi - 1) / br));
		Ie = -((Vbc * Is / _VT) * xi - (Vbe * Is / _VT) * (1 + 1 / bf) * zeta + Is * (zeta - xi + (zeta - 1) / bf));
	}
}

std::vector<int> BJT::getNodes() {
	std::vector<int> nodes;

	nodes.push_back(nodeCollector);
	nodes.push_back(nodeBase);
	nodes.push_back(nodeEmitter);

	return nodes;
}

std::complex<double> BJT::getConductance(int p_node1, int p_node2, double p_angularFrequency) {
	if (p_node1 == nodeCollector && p_node2 == nodeCollector) return Gcc;
	else if (p_node1 == nodeCollector && p_node2 == nodeBase) return Gcb;
	else if (p_node1 == nodeCollector && p_node2 == nodeEmitter) return Gce;
	else if (p_node1 == nodeBase && p_node2 == nodeCollector) return Gbc;
	else if (p_node1 == nodeBase && p_node2 == nodeBase) return Gbb;
	else if (p_node1 == nodeBase && p_node2 == nodeEmitter) return Gbe;
	else if (p_node1 == nodeEmitter && p_node2 == nodeCollector) return Gec;
	else if (p_node1 == nodeEmitter && p_node2 == nodeBase) return Geb;
	else if (p_node1 == nodeEmitter && p_node2 == nodeEmitter) return Gee;
	else return 0;
}

std::vector<double> BJT::getProperties() {
	std::vector<double> properties;

	properties.push_back(Ic);
	properties.push_back(Ib);
	properties.push_back(Ie);
	properties.push_back(npn);
	properties.push_back(Is);
	properties.push_back(bf);
	properties.push_back(br);

	return properties;
}

void BJT::setProperties(std::vector<double> properties) {
	if (npn == 1) {
		Vbe = properties[0];
		Vbc = properties[1];
	} else {
		Vbe = -properties[0];
		Vbc = -properties[1];
	}

	if (Vbe > 1) Vbe = 1;
	if (Vbc > 1) Vbc = 1;

	double zeta = exp(Vbe / _VT);
	double xi = exp(Vbc / _VT);

	Gcc = (Is / _VT) * (1 + 1 / br) * xi;
	Gcb = (Is / _VT) * (zeta - (1 + 1 / br) * xi);
	Gce = -(Is / _VT) * zeta;

	Gbc = -(Is / (_VT * br)) * xi;
	Gbb = (Is / _VT) * ((zeta / bf) + (xi / br));
	Gbe = -(Is / (_VT * bf)) * zeta;

	Gec = -(Is / _VT) * xi;
	Geb = -(Is / _VT) * ((1 + 1 / bf) * zeta - xi);
	Gee = (Is / _VT) * (1 + 1 / bf) * zeta;

	if (npn == 1) {
		Ic = (Vbe * Is / _VT) * zeta - (Vbc * Is / _VT) * (1 + 1 / br) * xi - Is * (zeta - xi - (xi - 1) / br);
		Ib = (Vbe * Is / (bf * _VT)) * zeta + (Vbc * Is / (br * _VT)) * xi - Is * ((zeta - 1) / bf + (xi - 1) / br);
		Ie = (Vbc * Is / _VT) * xi - (Vbe * Is / _VT) * (1 + 1 / bf) * zeta + Is * (zeta - xi + (zeta - 1) / bf);
	} else {
		Ic = -((Vbe * Is / _VT) * zeta - (Vbc * Is / _VT) * (1 + 1 / br) * xi - Is * (zeta - xi - (xi - 1) / br));
		Ib = -((Vbe * Is / (bf * _VT)) * zeta + (Vbc * Is / (br * _VT)) * xi - Is * ((zeta - 1) / bf + (xi - 1) / br));
		Ie = -((Vbc * Is / _VT) * xi - (Vbe * Is / _VT) * (1 + 1 / bf) * zeta + Is * (zeta - xi + (zeta - 1) / bf));
	}
}

MOSFET::MOSFET(std::string p_name, int p_nodeDrain, int p_nodeGate, int p_nodeSource, Model* model) :
	Component { p_name }, nodeDrain { p_nodeDrain }, nodeGate { p_nodeGate }, nodeSource { p_nodeSource } {
		std::vector<double> model_values = model->getDoubles();
		vto = model_values[0];
		k = model_values[1];
		nmos = model_values[2];
	}

std::vector<int> MOSFET::getNodes() {
	std::vector<int> nodes;
	nodes.push_back(nodeDrain);
	nodes.push_back(nodeGate);
	nodes.push_back(nodeSource);
	return nodes;
}

std::complex<double> MOSFET::getConductance(int p_node1, int p_node2, double p_angularFrequency) {
	return 0;
}

std::vector<double> MOSFET::getProperties() {
	std::vector<double> properties;
	properties.push_back(Id);
	properties.push_back(Ig);
	properties.push_back(Is);
	properties.push_back(Vgs);
	properties.push_back(Vds);
	properties.push_back(vto);
	properties.push_back(k);
	properties.push_back(nmos);
	return properties;
}

void MOSFET::setProperties(std::vector<double>) {

}