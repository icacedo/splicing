#include <iostream>
#include <fstream>
#include <string>
#include <vector>
// some guy on reddit said std::list has poor performance, never use
// some guy on stack overflow said malloc + array is not safe

std::string desc;
std::string seq;
std::string line;
std::vector<int> dons;
std::vector<int> accs;

void readFasta(std::string faFile);
void gtag(std::string seq);

int main(int argc, char** argv) 
{		
	readFasta(argv[1]);	
	std::cout << desc << "\n";
	std::cout << seq << "\n";
	gtag(seq);	
	// range based for loop
	// const makes variable immutable
	// & declares i as a reference
	for (const int& i : dons) {
		std::cout << i << "\n";
	} 
	std::cout << "#####\n";

	for (const int& i : accs) {
		std::cout << i << "\n";
	}
	return 0;
}

void readFasta(std::string faFile){
	std::ifstream MyReadFile(faFile);
	while (getline (MyReadFile, line)) {
		if (line.rfind(">", 0) == 0) {
			desc += line;
		}
		else {
			seq += line;
		}
	}
}

void gtag(std::string seq){
	for (int i = 0; i < seq.length(); i++) {
		if (seq.substr(i, 2) == "GT") {
			dons.push_back(i);
		};
		if (seq.substr(i, 2) == "AG") {
			accs.push_back(i);
		};
	}
}

