#include <iostream>

// need to go to preferences > settings > find code runner
// check run in terminal
int main() {
    int age;

    std::string name;

    std::cout << "What's your age?: ";
    std::cin >> age; // \n is inserted in input buffer

    std::cout << "What's your name?: ";
    // std::ws elimantes newline characters or white spaces before using input
    std::getline(std::cin >> std::ws, name);
    // cannot read in lines with spaces
    //std::cin >> name;

    // move to top
    //std::cout << "What's your age?: ";
    //std::cin >> age;

    std::cout << "Hello " << name;
    std::cout << " You are " << age << " years old";

    return 0;
}