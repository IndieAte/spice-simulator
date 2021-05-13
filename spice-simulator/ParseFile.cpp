#define _USE_MATH_DEFINES

#include "ParseFile.h"
#include <cmath>

//This function returns a vector of a string that has been split at each of a chosen character.
//The function also detects the character inside brackets and skips the character
std::vector<std::string> string_split(const std::string& s, char c) {
	int n = 0;
	bool skip = false;
	std::vector<std::string> v;

	v.push_back("");
	for (int i = 0; i < s.length(); i++) {
		if (s[i] != c) {
			if (s[i] == '(') skip = true;
			v[n] += s[i];
		} else if (skip) {
			v[n] += s[i];
			skip = false;
		} else if (!skip) {
			v.push_back("");
			n++;
		}
	}
	return v;
}

//This function takes a node name in the format: "N001" and returns the integer
//number that the node refers to.
//The function also takes a parameter n and checks if any of the nodes that are
//parsed through the function have a higher node number and sets n to it if it
//is higher.
int get_node_number(const std::string& s, int& n) {
	if (s.length() > 1) {
		int m = std::stoi(s.substr(1, 3));
		if (m > n) n = m;
		return m;
	} else {
		return 0;
	}
}

//Get value from integer and multiplier
double decode_value(std::string s) {
	std::string end = "";
	int counter = s.length() - 1;

	while (!isdigit(s[counter])) {
		end = s[counter] + end;
		counter--;
	}

	double d = std::stod(s.substr(0,counter+1));
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
		throw std::invalid_argument("Invalid Multiplier");
	}
}

//Get Amplitude and phase from AC String
std::vector<double> decode_ac(std::string s) {
	s.erase(0,3);
	s.erase(s.length()-1,1);
	std::vector<std::string> v1 = string_split(s,' ');
	
	std::vector<double> v2;
	v2.push_back(decode_value(v1[0]));
	v2.push_back(decode_value(v1[1]));
	return v2;
}

double degrees_to_radians(double d) {
	return (d * M_PI)/180;
}

double radians_to_degrees(double d) {
	return (d * 180) / M_PI;
}

//This function takes a file and returns a vector of Component pointers.
std::vector<Component*> decode_file(std::ifstream& infile, int& n, Command*& command) {
	std::vector<Component*> v1;
	std::string tmp;

	//Here the function iterates over each line and makes a component if necessary
	//with the data provided in the line.
	while (std::getline(infile, tmp)) {
		std::vector<std::string> v2 = string_split(tmp,' ');

		switch (toupper(tmp[0])) {
			case 'R': {
				if (v2[1] != v2[2]) v1.push_back(new Resistor(v2[0], decode_value(v2[3]), get_node_number(v2[1], n), get_node_number(v2[2], n)));
				break;
			}
			case 'I': {
				if (v2[1] != v2[2] && v2[3][0] == 'A') {
					std::vector<double> v3 = decode_ac(v2[3]);
					v1.push_back(new ACCurrentSource(v2[0], v3[0], degrees_to_radians(v3[1]), get_node_number(v2[1], n), get_node_number(v2[2], n)));
				} else if (v2[1] != v2[2]) {
					v1.push_back(new DCCurrentSource(v2[0], decode_value(v2[3]), get_node_number(v2[1], n), get_node_number(v2[2], n)));
				}
				break;
			}
			case 'C': {
				if (v2[1] != v2[2]) v1.push_back(new Capacitor(v2[0], decode_value(v2[3]), get_node_number(v2[1], n), get_node_number(v2[2], n)));
				break;
			}
			case 'L': {
				if (v2[1] != v2[2]) v1.push_back(new Inductor(v2[0], decode_value(v2[3]), get_node_number(v2[1], n), get_node_number(v2[2], n)));
				break;
			}
			case 'V': {
				if (v2[1] != v2[2] && v2[3][0] == 'A') {
					std::vector<double> v3 = decode_ac(v2[3]);
					v1.push_back(new ACVoltageSource(v2[0], v3[0], degrees_to_radians(v3[1]), get_node_number(v2[1], n), get_node_number(v2[2], n)));
				} else if (v2[1] != v2[2]) {
					v1.push_back(new DCVoltageSource(v2[0], decode_value(v2[3]), get_node_number(v2[1], n), get_node_number(v2[2], n)));
				}
				break;
			}
			case '.': {
				if (v2[0] == ".ac") {
					command = new ACCommand("AC", v2[1], decode_value(v2[2]), decode_value(v2[3]), decode_value(v2[4]));
				}
			}
		}
	}

	return v1;
}