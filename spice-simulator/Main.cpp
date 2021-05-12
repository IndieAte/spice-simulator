#include "ConductanceMatrix.h"
#include "Component.h"
#include "ParseFile.h"
#include "ACAnalysis.h"

int main() {
	std::ifstream infile;
	infile.open("../testCircuit.cir"); 
	if(!infile.is_open()){
		return EXIT_FAILURE;
	}

	//Here the highest node number is initialised and parsed into decode_file.
	//A vector of Component pointers is set to the output of decode_file
	int highest_node = 0;
	std::vector<Component*> components = decode_file(infile, highest_node);
	infile.close();

	int outputNode = 1;

	std::vector<Eigen::Vector3d> results = runACAnalysis(outputNode, 10, 1000, 10, components, highest_node);

	for (int i = 0; i < results.size(); i++) {
		double amplitude = results[i](0);
		double phase = radians_to_degrees(results[i](1));
		double frequency = results[i](2);

		std::cout << "At " << frequency << "Hz:" << std::endl;
		std::cout << amplitude << "V" << std::endl;
		std::cout << phase << " degrees" << std::endl << std::endl;
	}
}