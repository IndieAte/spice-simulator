#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Component.h"
#include "Node.h"


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

Component* decode_line(std::string line) {
	std::vector<std::string> v = string_split(line);

	switch(line[0]) {
		case 'R':
			std::cout << "R" << std::endl;
			break;
		case 'I':
			std::cout << "I" << std::endl;
			break;
	}
}


int main() {

	std::cout << "Test" << std::endl;


	std::ifstream infile;
	infile.open("../testCircuit.cir");
	std::vector<Component*> components;
 
	if(!infile.is_open()){ 
		return EXIT_FAILURE;
	}

	std::string tmp;
	while(std::getline(infile,tmp)){
		components.push_back(decode_line(tmp));
	}

	infile.close();

}