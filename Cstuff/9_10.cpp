// 9. Useful math related functions

#include <iostream>
#include <cmath>
/*
int main() 
{
    double x = 3;
    double y = 4;
    double z;

    //z = std::min(x, y);
    //z = std::max(x, y);
    //z = pow(2, 3); 
    //z = sqrt(9);
    z = abs(-7);

    // there are also round, ceil (round up), floor (round down) functions

    std::cout << z << '\n';

    return 0;
}
*/

// 10.Hypotenuse calculator 
/*
int main()
{
    double a;
    double b;
    double c;

    std::cout << "Enter side A: ";
    std::cin >> a;

    std::cout << "Enter side B: ";
    std::cin >> b;

    a = pow(a, 2);
    b = pow(b, 2);
    c = sqrt(a + b);

    std::cout << "side C: " << c << '\n';

    return 0;

}
*/

// 11.If statements
/*
int main()
{
    int age;

    std::cout << "Enter your age: ";
    std::cin >> age;

    if(age >= 18){
        std::cout << "Welcome";
    }
    else if(age < 0){
        std::cout << "Impossible";
    }
    else if(age >= 100){
        std::cout << "Very";
    }
    else{
        std::cout << "No";
    }

    return 0;
}
*/

// 12.Switches
// alternative to using many elif statements
// can compare one value against matching cases

/*
int main()
{
    int month;
    std::cout << "Enter the month (1-3): ";
    std::cin >> month;

    switch(month){
        case 1:
            std::cout << "It is January";
            break;
        case 2:
            std::cout << "It is February";
            break;
        case 3:
            std::cout << "It is March";
            break;
        default:
            std::cout << "Please enter in only numbers (1-3)";
    }
}

int main()
{
    char grade;

    std::cout << "What letter grade?: ";
    std::cin >> grade;

    switch(grade){
        case 'A':
            std::cout << "Nice";
            break;
        case 'C':
            std::cout << "Meh";
            break;
    }

    return 0;
}
*/

// 13.Console calculator program
/*
int main()
{
    char op;
    double num1;
    double num2;
    double result;

    std::cout << "*****CALC*****\n";

    std::cout << "Enter either (+ - * /): ";
    std::cin >> op;

    // too lazy to finish
    // but uses swith to compare cin and pick correct operator

    std::cout << "**************\n";

    return 0;
}
*/

// 14.ternary operator ?: = replacement to an if/else statement
// condition ? expression1 : expression2;

/*

if grad is greater than 60, do thing, else don't do 
grade >=60 ? std::cout << "Do thing" : std::cout << "Don't do";

*/

// 15.Logical operators
// && (and) check if two or more conditions are true
// if(temp > 0 && temp < 30){} both conditions must be true to execute
// || (or) check if at least one of two conditions is true
// if(temp <=0 || temp >= 30){} at least one of the condistions must be true
// ! (not) reverses the logical state of its operand
// bool sunny = true;
// if(!sunny){} checks if not sunny

// 17.Useful string methods



















