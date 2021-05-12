#define _USE_MATH_DEFINES

#include "ACAnalysis.h"
#include <cmath>

Eigen::Vector3d voltageVectorToAmpAndPhase(int outputNode, Eigen::VectorXcd voltageVector, double frequency) {
	Eigen::Vector3d output;

	double amplitude = abs(voltageVector(outputNode - 1));
	double phase = arg(voltageVector(outputNode - 1));

	output << amplitude, phase, frequency;

	return output;
}

std::vector<Eigen::Vector3d> runACAnalysis(int outputNode, double startFreq, double stopFreq, int pointsPerDecade,
	std::vector<Component*> components, int highest_node) {

	if (outputNode > highest_node) {
		std::cout << "Error: Output node not in netlist" << std::endl;
	}

	double currentFreq = startFreq;
	int n = 1;
	std::vector<Eigen::Vector3d> output;

	while (currentFreq <= stopFreq) {
		double currAngFreq = currentFreq * M_PI * 2;
		Eigen::VectorXcd voltageVector = solveAtFrequency(components, highest_node, currAngFreq);
		output.push_back(voltageVectorToAmpAndPhase(outputNode, voltageVector, currentFreq));

		double exponent = n / (double)pointsPerDecade;
		currentFreq = pow(10, exponent) * startFreq;
		n++;
	}

	return output;
}