#include <iostream>
#include <vector>


//typedef and using section
/*
// pairlist is the alias for the new data type
//typedef std::vector<std::pair<std::string, int>> pairlist_t;
//instead of typedef, it is better to use the using keyword
//typedef std::string text_t;
//typedef int number_t;
using text_t = std::string;
using number_t = int;


int main(){

    using std::cout;

    // typedef is a reserved keyword used to create alias for another data type
    // used to increase readability

    //pairlist_t pailist;

    text_t firstname = "Snake";
    number_t age = 234;
    
    cout << firstname << '\n';
    cout << age;

  

    return 0;
}
*/

//arithmetic operators

//order of operations:
//parenthesis, multiplication and division, then addition and subtraction

int main(){

    // need int to use modulus operator
    int students = 20;
    // if division creates decimals, use double
    //double students = 20;

    //students = students + 1;
    students+=2;
    //increment operator, adds 1
    students++;

    //decrement operator, subtracts 1
    students--;

    //students = students * 2;
    students*=2;

    students/=3;

    //modulus
    int remainder = students % 3;

    std::cout << students << '\n';
    std::cout << remainder << '\n';

    //type conversion, from one data type to another

    //int x = 3.14;
    // converts 3.14 to an integer and stores it in the double
    double x = (int) 3.14;

    std::cout << x << '\n';

    int correct = 8;
    int questions = 10;
    double score = correct/(double)questions * 100;
    //double score = correct/questions * 100;

    std::cout << score << "%";

    return 0;
}