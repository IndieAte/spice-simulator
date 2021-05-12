#include "ConductanceMatrix.h"
#include "Component.h"
#include "ParseFile.h"

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

	Eigen::VectorXcd solution = solveAtFrequency(components, highest_node, 0);

	for (int i = 0; i < highest_node; i++) {
		std::cout << "Node " << i + 1 << ": ";
		std::cout << solution(i) << std::endl;
	}
}