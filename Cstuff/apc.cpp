#include <iostream>
#include <fstream>
#include <string>
#include <vector>
// some guy on reddit said std::list has poor performance, never use
// some guy on stack overflow said malloc + array is not safe
// i need a way to handle problems like forgetting a input file

std::string desc;
std::string seq;
std::string line;
//int flank = 99;
//int emin = 25;
//int imin = 35;
int flank = 3;
int emin = 3;
int imin = 3;
std::vector<int> dons;
std::vector<int> accs;
//std::vector<int> introns;
std::vector<std::vector<int>> introns;

void readFasta(std::string faFile);
void gtag(std::string seq);
void apc(const std::vector<int>& dons, const std::vector<int>& accs, 
		//std::vector<int> introns);
		std::vector<std::vector<int>> introns);
int main(int argc, char** argv) 
{		
	readFasta(argv[1]);		
	gtag(seq);	
	apc(dons, accs, introns);	
	
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
	for (int i = emin+flank; i < seq.length()-flank-emin; i++) {
		if (seq.substr(i, 2) == "GT") {
			dons.push_back(i);
		};
		if (seq.substr(i, 2) == "AG") {
			accs.push_back(i);
		};
	}
}

// the tuple plus vector thing is very hard to do in C++
// look into using std::pair for tuples
// not really sure when const and & are needed
void apc(const std::vector<int>& dons, const std::vector<int>& accs, 
		//std::vector<int> introns){
		std::vector<std::vector<int>> introns){
	int don = dons[0];
	for (int aix = 0; aix < accs.size(); aix++) {
		if (accs[aix] - don + 1 < imin) {
			continue;
		}
		std::vector<int> intron = {don, accs[aix]};
		introns.push_back(intron);
	}
	// fun with vectors
	std::vector<std::vector<int>> vec1;
	std::vector<int> vec2;
	vec1.push_back(vec2);
}

