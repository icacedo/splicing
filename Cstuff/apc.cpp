#include <iostream>
#include <fstream>

std::string desc;
std::string line;

void readFasta(std::string faFile);

// argc is a count of arguments on the line
// argv is a character array
// char* is a pointer to a char
// char** is a point to a pointer of char
// char** is a memory address with char*
int main(int argc, char** argv) 
{		
	readFasta(argv[1]);

	return 0;
}

void readFasta(std::string faFile){
	std::ifstream MyReadFile(faFile);
	while (getline (MyReadFile, line)) {
		if (line.rfind(">", 0) == 0) {
			desc += line;
			std::cout << desc << "\n";
		}
		else {
			std::cout << line << "\n";
		}
	}
}






