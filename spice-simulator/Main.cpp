#include "Component.h"
#include "ParseFile.h"
#include "ACAnalysis.h"

int main(int argc, char** argv) {
	std::ifstream infile;
	std::string outfilePath;
	int outputNode = -1;
	
	try {
		if (argc == 1) {
			throw std::invalid_argument("Error: No input file provided");
		} else if (argc == 2) {
			if (std::string(argv[1]).compare("-h") != 0) {
				std::cout << argv[1] << std::endl;
				infile.open(argv[1]);
				outputNode = -1;
				outfilePath = "output.csv";
			} else {
				std::cout << std::endl;
				std::cout << "There are 4 valid syntaxes to call this program:" << std::endl;
				std::cout << "  circuit-sim INPUT_FILE_PATH" << std::endl;
				std::cout << "  circuit-sim INPUT_FILE_PATH OUTPUT_NODE" << std::endl;
				std::cout << "  circuit-sim INPUT_FILE_PATH OUTPUT_FILE_PATH" << std::endl;
				std::cout << "  circuit-sim INPUT_FILE_PATH OUTPUT_NODE OUTPUT_FILE_PATH" << std::endl << std::endl;
				std::cout << "INPUT_FILE_PATH should be the path to the input netlist (typically a .cir file)" << std::endl;
				std::cout << "OUTPUT_NODE should be the number of the node to be treated as the output for an AC analysis, if not given you will be prompted to enter it later" << std::endl;
				std::cout << "OUTPUT_FILE_PATH should be the path to the output .csv file, if not given the default is 'output.csv'" << std::endl;
				return 0;
			}
		} else if (argc == 3) {
			if (is_number(argv[2], false)) {
				outputNode = std::stoi(argv[2]);
				outfilePath = "output.csv";
			} else {
				outfilePath = argv[2];
			}
			infile.open(argv[1]);
		} else if (argc == 4) {
			if (is_number(argv[2], false)) {
				outputNode = std::stoi(argv[2]);
			} else {
				throw std::invalid_argument("Output Node provided is not an integer.");
			}
			infile.open(argv[1]);
			outfilePath = argv[3];
		}

		if(!infile.is_open()){
			throw std::invalid_argument("Failed to open Input File.");
		}
	} catch (std::invalid_argument& e) {
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;	
	}

	//Here the highest node number is initialised and parsed into decode_file.
	//A vector of Component pointers is set to the output of decode_file
	int nNodes = 0;
	Command* command;
	std::vector<int> ac_source_indexes;
	std::vector<Component*> components = decode_file(infile, nNodes, command, ac_source_indexes);
	infile.close();

	std::vector<double> command_values = command->getValues();

	std::ofstream outfile(outfilePath);

	if (outfile.is_open()) {
		if (command->type == "AC") {
			try {
				if (ac_source_indexes.size() == 0) throw std::invalid_argument("No AC Sources in circuit.");

				std::string message = "No Output Node was provided. Please enter one: ";
				if (outputNode > nNodes) {
					outputNode = -1;
					message = "Ouput node provided is outside range of nodes. Please enter a new one: ";
				}
				while (outputNode == -1) {
					std::string tmp;
					std::cout << message;
					std::cin >> tmp;
					if (is_number(tmp, false)) {
						outputNode = std::stoi(tmp);
						if (outputNode > nNodes) {
							outputNode = -1;
							message = "Ouput node provided is outside range of nodes. Please enter a new one: ";
						}
					} else {
						message = "Input provided was not an integer. Please enter an Output Node: ";
					}
				}

				int input_index;
				if (ac_source_indexes.size() > 1) {
					std::string sources = "";
					for (int i=0; i<ac_source_indexes.size(); i++) {
						Component* ac_component = components[ac_source_indexes[i]];
						sources += "\t" + std::to_string(i+1) + ": " + ac_component->getName() + "\n";
					}
					std::cout << "More than one AC Source provided. Please enter which one you want to take as input:\n" << sources;
					std::cin >> input_index;
				}
				else input_index = ac_source_indexes[0];

				std::vector<std::vector<Eigen::Vector3d>> results = runACAnalysis(outputNode, input_index, command_values[2], command_values[3], command_values[1], components, nNodes);

				outfile << "Frequency / Hz, Amplitude / dB, Phase / Degrees" << std::endl;
				for (int i = 0; i < results.size(); i++) {
					double amplitude = 20 * log10(results[i][0](0)/results[i][1](0));
					double phase = radians_to_degrees(results[i][0](1)-results[i][1](1));
					double frequency = results[i][0](2);

					outfile << frequency << ", ";
					outfile << amplitude << ", ";
					outfile << phase << ", " << std::endl;
				}
			} catch (std::invalid_argument& e) {
				std::cerr << e.what() << std::endl;
				return EXIT_FAILURE;
			}
		} else if (command->type == "OP") {
			Eigen::VectorXd results = runDCOpPoint(components, nNodes);
			for (int i=0; i<nNodes; i++) {
				outfile << "V(N" << std::to_string(i+1) << "): " << results(i) << std::endl;
			}
		}
	} else {
		std::cout << "Failed to open output.csv" << std::endl;
		return EXIT_FAILURE;
	}
}