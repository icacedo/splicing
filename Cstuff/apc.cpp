#include <iostream>
#include <fstream>

std::string desc;
std::string seq;
std::string line;

void readFasta(std::string faFile);
int gtag() {
	int num = 42;
	return num;
} 

int main(int argc, char** argv) 
{		
	readFasta(argv[1]);	
	std::cout << desc << "\n";
	std::cout << seq << "\n";
	int result = gtag();
	std::cout << result << "\n";
	
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





