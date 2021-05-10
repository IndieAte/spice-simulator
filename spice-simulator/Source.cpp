#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Component.h"

std::vector<std::string> string_split(const std::string& s) {
	int c = 0;
	std::vector<std::string> v;

	v.push_back("");
	for (int i=0; i<s.length(); i++) {
		if (s[i] !=  ' ') {
			v[c] += s[i];
		} else {
			v.push_back("");
			c++;
		}
	}
	return v;
}
int get_node_number(const std::string& s) {
	if (s.length() > 1) {
		return std::stoi(s.substr(1,4));
	} else {
		return 0;
	}
}

Component decode_line(std::string line) {
	std::vector<std::string> v = string_split(line);

	// Component comp;


	switch(line[0]) {
		case 'R': {
			std::cout << "R" << get_node_number(v[1]) << get_node_number(v[2]) << std::endl;

			Resistor r1(v[0], std::stod(v[3]), get_node_number(v[1]), get_node_number(v[2]));
			return r1;
			break;
		}
		case 'I': {
			std::cout << "I" << std::endl;

			Resistor r2(v[0], std::stod(v[3]), get_node_number(v[1]), get_node_number(v[2]));
			return r2;

			break;
		}
	}	

	Resistor r3(v[0], std::stod(v[3]), get_node_number(v[1]), get_node_number(v[2]));
	return r3;

}


int main() {

	// std::cout << "Test" << std::endl;


	// std::ifstream infile;
	// infile.open("../testCircuit.cir");
	// std::vector<Component> components;
 
	// if(!infile.is_open()){
	// 	return EXIT_FAILURE;
	// }

	// std::string tmp;
	// while(std::getline(infile,tmp)){
	// 	components.push_back(decode_line(tmp));
	// }

	// infile.close();

	Resistor r3("hi", 12, 1, 2);
	r3.getConductance(1,2);
}