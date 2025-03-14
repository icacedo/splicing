#include <iostream>

namespace first{
    int x = 1;
}
namespace second{
    int x = 2;
}

int main()
{
    // Hello World
    // standard character output
    // semi colon ends statements
    /*
    std::cout << "I like donuts" << std::endl;
    std::cout << "They are full of sugar..." << '\n';
    std::cout << "wow";
    */
    
    // Variables
    /*
    int x; //declaration
    x = 5; //assignment
    int y = 6;
    int sum = x + y;

    std::cout << y << '\n';
    std::cout << x << '\n';
    std::cout << sum;
    */

    /*
    //integer
    int age = 21;
    int year = 2025;
    // if declared as int, will be displayed as 7
    int days = 7.5;

    //doubles include decimals
    double price = 10.99;
    double gpa = 2.5;
    double temperature = 25.1;

    // char data type stores a single character type
    char grade = 'A';
    char initial = 'B';
    // can only store a single character, not 'BC'
    char currency = '$';

    //boolean
    bool student = true;
    bool power = true;
    bool forSale = false;

    //strings
    std::string name = "Snake";

    std::cout << price << '\n';
    std::cout << days << '\n';
    std::cout << initial << '\n';
    std::cout << forSale << '\n';
    std::cout << "Hello " << name << '\n';
    std::cout << "You are " << age << "years old";
    */

    // const keyword, specifies that a variable's value is constant
    // compiler will not modify that variable, it is read only
    // show as all uppercase
    // use constants as often as possible, unless you know the variable will be changed

    /*
    const double PI = 3.14159;
    double radius = 10;
    double circumference = 2 * PI * radius;
    const int WIDTH = 1920;
    const int HEIGHT = 1080;

    std::cout << circumference << "cm";
    */

    //namespaces
    //prevents name conflicts, can allow for identically named entities so long as 
    //namespace is different
    //needs to be created outside of the main function
    //:: is the scope resolution operator
    //can refer to a certain version of x

    std::cout << first::x << '\n';

    using namespace second;
    // if using namespace std, don't need to include std in prefixes
    // but std namespace has hundreds of entities
    // can specify exactly what you want to use
    // but usually don't want to uses namespace std
    using std::cout;

    cout << x << '\n';
    cout << first::x;


    // if 0 is not returned, program did not execute properly
    return 0;
}