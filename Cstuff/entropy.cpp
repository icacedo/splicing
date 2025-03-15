#include <iostream>
#include <fstream>
#include <vector>
#include <typeinfo>
// uisng std includes hundreds of functions, names may clash
using namespace std;

// argc counts the number of command line arguments
// argv is an array of char* with the command line arguments
int main(int argc, char* argv[]) {

    // create text string to output text file
    string line;
    ifstream MyReadFile(argv[1]);

    vector<char> sequence;

    

    while (getline (MyReadFile, line)) {
        if (line[0] == '>') {
            cout << line << '\n';
        } else {
            cout << line.substr(0, 5);
        }
    }

    cout << "wow";

    MyReadFile.close();

    return 0;
}

