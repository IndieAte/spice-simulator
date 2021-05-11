#include "ParseFile.h"

//This function returns a vector of a string that has been split at each space.
std::vector<std::string> string_split(const std::string& s) {
	int c = 0;
	std::vector<std::string> v;

	v.push_back("");
	for (int i = 0; i < s.length(); i++) {
		if (s[i] != ' ') {
			v[c] += s[i];
		} else {
			v.push_back("");
			c++;
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

//This function takes a file and returns a vector of Component pointers.
std::vector<Component*> decode_file(std::ifstream& infile, int& n) {
	std::vector<Component*> v1;
	std::string tmp;

	//Here the function iterates over each line and makes a component if necessary
	//with the data provided in the line.
	while (std::getline(infile, tmp)) {
		std::vector<std::string> v2 = string_split(tmp);
		switch (tmp[0]) {
		case 'R': {
			if (v2[1] != v2[2]) v1.push_back(new Resistor(v2[0], std::stod(v2[3]), get_node_number(v2[1], n), get_node_number(v2[2], n)));
			break;
		}
		case 'I': {
			if (v2[1] != v2[2]) v1.push_back(new DCCurrentSource(v2[0], std::stod(v2[3]), get_node_number(v2[1], n), get_node_number(v2[2], n)));
			break;
		}
		}
	}

	return v1;
}