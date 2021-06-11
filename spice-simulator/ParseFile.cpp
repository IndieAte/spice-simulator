#define _USE_MATH_DEFINES

#include "ParseFile.h"

// This function returns a vector of a string that has been split at a chosen character.
// It only splits once even when it comes to multiple of the character in a row.
std::vector<std::string> string_split(const std::string& s, char split) {
	int n = 0;
	std::vector<std::string> v;

	v.push_back("");
	for (int i = 0; i < s.length(); i++) {
		if (s[i] != split) {
			v[n] += s[i];
		} else if (s[i] == split && v[n] != "") {
			v.push_back("");
			n++;
		}
	}
	return v;
}

// This function checks if a string is entirely a number.
// It also can check if the string int or decimal
bool is_number(std::string s, bool dec_check) {
	int counter = 0;
	for (int i=0; i<s.length(); i++) if (isdigit(s[i]) || (s[i] == '-' && i==0) || (s[i] == '.' && dec_check)) counter++;
	return counter == s.length();
}

// This function takes a node name in the format: "N001" and returns the integer
// number that the node refers to.
// The function also takes a parameter n and checks if any of the nodes that are
// parsed through the function have a higher node number and sets n to it if it
// is higher.
// The function also takes a vector node_count and counts how many components are
// connected to each node.
int get_node_number(const std::string& s, int& nNodes, std::vector<int>& node_count, bool count_node) {
	try {
		if (s.length() > 1) {
			std::string node_number_string = s.substr(1, 3);
			if (is_number(node_number_string, false)) {
				int node_number = std::stoi(node_number_string);
				if (node_number > nNodes && count_node) {
					nNodes = node_number;
					for (int i=node_count.size(); i<node_number+1; i++) node_count.push_back(0);
				}
				if (count_node) {
					node_count[node_number]++;
				}
				return node_number;
			} else {
				throw std::invalid_argument("Invalid Node: " + s);
			}
		} else if (s == "0") {
			if (count_node) {
				if (node_count.size() == 0) node_count.push_back(0);
				node_count[0]++;
			}
			return 0;
		} else {
			throw std::invalid_argument("Invalid Node: " + s);
		}
	} catch (std::invalid_argument& e) {
		std::cerr << e.what() << std::endl;
	}
}

// This function gets a numerical value from a string containing a number and multiplier.
double decode_value(std::string s) {
	std::string end = "";
	int counter = s.length() - 1;

	while (!isdigit(s[counter])) {
		end = s[counter] + end;
		counter--;
	}

	std::string front = s.substr(0,counter+1);
	std::cout << front << std::endl;
	try {
		if (is_number(front, true)) {
			double d = std::stod(front);
			std::cout << d << std::endl << std::endl;
			if (end == "p") {
				return d * 0.000000000001;
			} else if (end == "n") {
				return d * 0.000000001;
			} else if (end == "u") {
				return d * 0.000001;
			} else if (end == "m") {
				return d * 0.001;
			} else if (end == "k") {
				return d * 1000;
			} else if (end == "Meg") {
				return d * 1000000;
			} else if (end == "G") {
				return d * 1000000000;
			} else if (end == "") {
				return d;
			} else {
				throw std::invalid_argument("Invalid Multiplier: " + end);
			}
		} else {
			throw std::invalid_argument("Invalid Value: " + s);
		}
	} catch (std::invalid_argument& e) {
		std::cerr << e.what() << std::endl;
	}
}

// This function removes the AC and brackets from the amplitude and phase inputs.
std::vector<double> decode_ac(std::string a, std::string p) {
	a.erase(0,3);
	p.erase(p.length()-1,1);
	std::vector<double> v2;
	v2.push_back(decode_value(a));
	v2.push_back(decode_value(p));
	return v2;
}

// This function converts degrees to radians.
double degrees_to_radians(double d) {
	return (d * M_PI)/180;
}

// This function converts radians to degrees.
double radians_to_degrees(double d) {
	return (d * 180) / M_PI;
}

// This function converts a sweep into its numerical value.
double decode_sweep(std::string sweep) {
	if (sweep == "dec") {
		return 10;
	} else if (sweep == "oct") {
		return 8;
	} else if (sweep == "lin") {
		return 1;
	} else {
		throw std::invalid_argument("Invalid Sweep: " + sweep);
	}
}

// This function returns a vector of what is inside a set of brackets split at the spaces.
std::vector<std::string> open_brackets(std::string s) {
	std::vector<std::string> v;
	v.push_back("");
	int s_counter = 0, v_counter = 0;

	while (s[s_counter] != '(') {
		if (s[s_counter] != ' ') v[v_counter] += s[s_counter];
		s_counter++;
	}
	v_counter++;
	s_counter++;

	std::string s2 = s.substr(s_counter,s.length()-s_counter-1);
	std::vector<std::string> v2 = string_split(s2,' ');

	for (int i=0; i<v2.size(); i++) {
		v.push_back(v2[i]);
	}
	return v;
}

// This function takes a vector and index and returns a string corrosponding to the concatenation
// of the index provided of the vector and the rest of the vector elements after the index.
std::string get_final_elements(int index, std::vector<std::string> v) {
	std::string s = "";
	for (int i = index; i < v.size(); i++) {
		if (i == 0) s += v[i];
		else s += " " + v[i];
	}
	return s;
}

// This function generates a model for a specified component from a vector of the command line .model split at spaces.
// From this it will return a model of default values for parameters that are not specified
// and values that are specified will be set to the value provided.
Model* create_model(std::vector<std::string> v) {
	std::vector<std::string> end_v = open_brackets(get_final_elements(1,v));

	if (end_v[0] == "D") {
		double Is = pow(10,-12);
		for (int i=1; i<end_v.size(); i++) {
			std::vector<std::string> values = string_split(end_v[i],'=');

			std::transform(values[0].begin(), values[0].end(), values[0].begin(),
				[](unsigned char c) { return std::tolower(c); });

			if (values[0] == "is") {
				Is = decode_value(values[1]);
			}
		}
		return new DModel(v[0], end_v[0], Is);
	} else if (end_v[0] == "NPN" || end_v[0] == "PNP") {
		double Is = pow(10, -12), bf = 100, br = 1, npn = 1;
		if (end_v[0] == "PNP") npn = 0;
		double vaf = 10000, var = 10000;
		double cjc = 0, vjc = 0.75, mjc = 0.33;
		double cje = 0, vje = 0.75, mje = 0.33;
		double fc = 0.5;

		for (int i = 1; i < end_v.size(); i++) {
			std::vector<std::string> values = string_split(end_v[i],'=');

			std::transform(values[0].begin(), values[0].end(), values[0].begin(),
				[](unsigned char c) { return std::tolower(c); });

			if (values[0] == "is") {
				Is = decode_value(values[1]);
			} else if (values[0] == "bf") {
				bf = decode_value(values[1]);
			} else if (values[0] == "br") {
				br = decode_value(values[1]);
			} else if (values[0] == "vaf") {
				vaf = decode_value(values[1]);
			} else if (values[0] == "var") {
				var = decode_value(values[1]);
			} else if (values[0] == "cjc") {
				cjc = decode_value(values[1]);
			} else if (values[0] == "vjc") {
				vjc = decode_value(values[1]);
			} else if (values[0] == "mjc") {
				mjc = decode_value(values[1]);
			} else if (values[0] == "cje") {
				cje = decode_value(values[1]);
			} else if (values[0] == "vje") {
				vje = decode_value(values[1]);
			} else if (values[0] == "mje") {
				mje = decode_value(values[1]);
			} else if (values[0] == "fc") {
				fc = decode_value(values[1]);
			}
		}
		return new QModel(v[0], "Q", Is, bf, br, vaf, var, npn, cjc, vjc, mjc, cje, vje, mje, fc);
	} else if (end_v[0] == "NMOS" || end_v[0] == "PMOS") {
		double vto = 2.9, k = 0.005, va = 100, nmos = 1;
		if (end_v[0] == "PMOS") nmos = 0;

		for (int i=1; i<end_v.size(); i++) {
			std::vector<std::string> values = string_split(end_v[i],'=');

			std::transform(values[0].begin(), values[0].end(), values[0].begin(),
				[](unsigned char c) { return std::tolower(c); });

			if (values[0] == "vto") {
				vto = decode_value(values[1]);
			} else if (values[0] == "k") {
				k = decode_value(values[1]);
			} else if (values[0] == "va") {
				va = decode_value(values[1]);
			}
		}
		return new MModel(v[0], "M", vto, k, va, nmos);
	}
}

// This function checks if a there is a model already specified with the same name and component type
// as were passed into it. If there is no model that had been specified from the .cir file, then it will
// return the default model for the component type.
Model* get_model(std::string model_name, std::string model_type, std::vector<Model*>& models) {
	for (int i=0; i<models.size(); i++) {
		if (model_name == models[i]->name && model_type == models[i]->component) return models[i];
	}

	std::cout << "No model matching: '" << model_name << "'. Creating a new default model." << std::endl;

	std::vector<std::string> tmp;
	tmp.push_back(model_name);

	std::string model_component;
	if (model_type == "D") model_component = "D";
	else if (model_type == "Q") {
		if (model_name == "PNP") model_component = "PNP";
		else model_component = "NPN";
	}
	else if (model_type == "M") {
		if (model_name == "PMOS") model_component = "PMOS";
		else model_component == "NMOS";
	}
	tmp.push_back(model_component);
	
	tmp.push_back("()");
	Model* new_model = create_model(tmp);

	models.push_back(new_model);
	return new_model;
}

// This function takes a file and returns a vector of Component pointers that are specified by
// the .cir file.
std::vector<Component*> decode_file(std::ifstream& infile, int& nNodes, Command*& command, std::vector<int>& ac_source_indexes) {

	// This creates a vector of a vector with each element being a vector of each line split at spaces.
	std::string tmp;
	std::vector<std::vector<std::string>> file_vector;
	while (std::getline(infile, tmp)) {
		file_vector.push_back(string_split(tmp,' '));	
	}

	// This creates a vector of models which have been specified by the .cir file.
	std::vector<Model*> models;
	for (int i=0; i<file_vector.size(); i++) {
		std::vector<std::string> line_vector = file_vector[i];
		if (line_vector[0] == ".model") {
			line_vector.erase(line_vector.begin());
			models.push_back(create_model(line_vector));
		}
	}


	// This adds each component to the vector of components.
	std::vector<Component*> comps;
	std::vector<int> node_count;
	for (int i=0; i<file_vector.size(); i++) {
		std::vector<std::string> line_vector = file_vector[i];
		// std::cout << line_vector[0] << std::endl;
		try {
			switch (toupper(line_vector[0][0])) {
				case 'R': {
					if (line_vector[1] != line_vector[2] && line_vector.size() == 4) {
						comps.push_back(new Resistor(line_vector[0], decode_value(line_vector[3]), 
							get_node_number(line_vector[1], nNodes, node_count, true), 
							get_node_number(line_vector[2], nNodes, node_count, true)));
					} else if (line_vector[1] != line_vector[2]) {
						throw std::invalid_argument("Invalid Formatting of Resistor: " + line_vector[0]);
					}
					break;
				}
				case 'I': {
					if (line_vector[1] != line_vector[2] && line_vector[3][0] == 'A' && line_vector.size() == 5) {
						ac_source_indexes.push_back(comps.size());
						std::vector<double> ac_values = decode_ac(line_vector[3],line_vector[4]);
						comps.push_back(new ACCurrentSource(line_vector[0], ac_values[0], 
							degrees_to_radians(ac_values[1]), 
							get_node_number(line_vector[1], nNodes, node_count, true), 
							get_node_number(line_vector[2], nNodes, node_count, true)));
					} else if (line_vector[1] != line_vector[2] && line_vector.size() == 4) {
						comps.push_back(new DCCurrentSource(line_vector[0], decode_value(line_vector[3]), 
							get_node_number(line_vector[1], nNodes, node_count, true), 
							get_node_number(line_vector[2], nNodes, node_count, true)));
					} else if (line_vector[1] != line_vector[2]) {
						throw std::invalid_argument("Invalid Formatting of Current Source: " + line_vector[0]);
					}
					break;
				}
				case 'C': {
					if (line_vector[1] != line_vector[2] && line_vector.size() == 4) {
						comps.push_back(new Capacitor(line_vector[0], decode_value(line_vector[3]), 
							get_node_number(line_vector[1], nNodes, node_count, true), 
							get_node_number(line_vector[2], nNodes, node_count, true)));
					} else if (line_vector[1] != line_vector[2]) {
						throw std::invalid_argument("Invalid Formatting of Capacitor: " + line_vector[0]);
					}
					break;
				}
				case 'L': {
					if (line_vector[1] != line_vector[2] && line_vector.size() == 4) {
						comps.push_back(new Inductor(line_vector[0], decode_value(line_vector[3]), 
							get_node_number(line_vector[1], nNodes, node_count, true), 
							get_node_number(line_vector[2], nNodes, node_count, true)));
					} else if (line_vector[1] != line_vector[2]) {
						throw std::invalid_argument("Invalid Formatting of Inductor: " + line_vector[0]);
					}
					break;
				}
				case 'V': {
					if (line_vector[1] != line_vector[2] && line_vector[3][0] == 'A' && line_vector.size() == 5) {
						ac_source_indexes.push_back(comps.size());
						std::vector<double> ac_values = decode_ac(line_vector[3],line_vector[4]);
						comps.push_back(new ACVoltageSource(line_vector[0], ac_values[0], degrees_to_radians(ac_values[1]), 
							get_node_number(line_vector[1], nNodes, node_count, true), 
							get_node_number(line_vector[2], nNodes, node_count, true)));
					} else if (line_vector[1] != line_vector[2] && line_vector.size() == 4) {
						comps.push_back(new DCVoltageSource(line_vector[0], decode_value(line_vector[3]), 
							get_node_number(line_vector[1], nNodes, node_count, true), 
							get_node_number(line_vector[2], nNodes, node_count, true)));
					} else if (line_vector[1] != line_vector[2]) {
						throw std::invalid_argument("Invalid Formatting of Voltage Source: " + line_vector[0]);
					}
					break;
				}
				case 'D': {
					if (line_vector[1] != line_vector[2] && line_vector.size() == 4) {
						comps.push_back(new Diode(line_vector[0], 
							get_node_number(line_vector[1], nNodes, node_count, true), 
							get_node_number(line_vector[2], nNodes, node_count, true), 
							get_model(line_vector[3], "D", models)));
					} else if (line_vector[1] != line_vector[2]) {
						throw std::invalid_argument("Invalid Formatting of Diode: " + line_vector[0]);
					}
					break;
				}
				case 'Q': {
					if (line_vector.size() == 5) {
						comps.push_back(new BJT(line_vector[0], 
							get_node_number(line_vector[1], nNodes, node_count, true), 
							get_node_number(line_vector[2], nNodes, node_count, true), 
							get_node_number(line_vector[3], nNodes, node_count, true), 
							get_model(line_vector[4], "Q", models)));
					} else {
						throw std::invalid_argument("Invalid Formatting of BJT: " + line_vector[0]);
					}
					break;
				}
				case 'G': {
					if (line_vector[1] != line_vector[2] && line_vector.size() == 6) {
						comps.push_back(new VoltageControlledCurrentSource(line_vector[0], 
							decode_value(line_vector[5]), 
							get_node_number(line_vector[1], nNodes, node_count, true), 
							get_node_number(line_vector[2], nNodes, node_count, true), 
							get_node_number(line_vector[3], nNodes, node_count, false), 
							get_node_number(line_vector[4], nNodes, node_count, false)));
					} else if (line_vector[1] != line_vector[2]) {
						throw std::invalid_argument("Invalid Formatting of Voltage Controlled Current Source: " + line_vector[0]);
					}
					break;
				}
				case 'M': {
					if (line_vector.size() == 5) {
						comps.push_back(new MOSFET(line_vector[0], 
							get_node_number(line_vector[1], nNodes, node_count, true), 
							get_node_number(line_vector[2], nNodes, node_count, true), 
							get_node_number(line_vector[3], nNodes, node_count, true), 
							get_model(line_vector[4], "M", models)));
					} else {
						throw std::invalid_argument("Invalid Formatting of MOSFET: " + line_vector[0]);
					}
					break;
				}
				case '.': {
					if (line_vector[0] == ".ac") {
						if (line_vector.size() == 5) {
							command = new ACCommand("AC", decode_sweep(line_vector[1]), 
								decode_value(line_vector[2]), decode_value(line_vector[3]), 
								decode_value(line_vector[4]));
						} else {
							throw std::invalid_argument("Invalid Formatting of AC Command");
						}
					} else if (line_vector[0] == ".op") {
						if (line_vector.size() == 1) {
							command = new OPCommand("OP");
						} else {
							throw std::invalid_argument("Invalid Formatting of OP Command");
						}
					} else if (line_vector[0] == ".end") {
						return comps;
					}
					break;
				}
			}
		} catch (std::invalid_argument& e) {
			std::cerr << e.what() << std::endl;
		}
	}

	// This finds any floating nodes;
	for (int i=0; i<node_count.size(); i++) {
		try {
			if (node_count[i] == 1) throw std::invalid_argument("Found Floating Node: Node "+std::to_string(i));
		} catch (std::invalid_argument& e) {
			std::cerr << e.what() << std::endl;
		}
	}

	return comps;
}