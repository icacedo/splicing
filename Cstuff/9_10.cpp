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
/*
int main()
{
    std::string name;

    std::cout << "Enter your name: ";
    std::getline(std::cin, name);

    if(name.empty()){
        std::cout << "You didn't enter your name\n";
    }

    if(name.length() > 12){
        std::cout << "Your name can't be over 12 characters\n";
    }
    else{
        std::cout << "Welcome " << name << "\n";
    }

    // use name.clear() to clear the name

    name.append("@gmail.com");

    std::cout << "Your username is now " << name << "\n";

    // get character at a particular position
    std::cout << name.at(0);

    // insert string at certain position
    name.insert(0, "@@@");

    std::cout << name;

    // get positon of a character
    // name.find(' ') where is there a white space?
    // name.erase(0, 3) deletes the first 3 characters, not inclusive

    return 0;
}
*/

// 18.While loops
/*
int main()
{
    std::string name;

    // can only continue with program if 
    // name is not empty
    while(name.empty()){
        std::cout << "Enter your name: ";
        std::getline(std::cin, name);
    }

    std::cout << "Hello" << name;
}
*/

// 19.Do while loops
/*
int main()
{
    // do while loop = do some block of code first
    // then repeat again if condition is true

    int number;

    // with a do while loop, don't need these here
    //std::cout << "Enter a positive #: " << std::endl;
    //std::cin >> number;

    do{
        std::cout << "Enter a positive #: " << std::endl;
        std::cin >> number;
    }while(number < 0);

    std::cout << "The # is: " << number << std::endl;

    return 0;

}
*/

// 20.For loops
/*
int main()
{
    // first part is index, then a stopping condition
    // need to increment your index
    // can increment by different amounts
    // i += 2 increments by 2
    // can also decrment i-- or i-=2
    for(int i = 1; i <= 3; i++){
        std::cout << i << std::endl;
    }

    std::cout << "WOW\n";

    return 0;
}
*/

// 21.Break and continue
// beak out of a loop
// skip current iteration of loop
// use an if statement and continue to skip or break

// 22.Nested loops
/*
int main()
{
    int rows;
    int columns;
    char symbol;

    std::cout << "How many rows?: \n";
    std::cin >> rows;

    std::cout << "How many columns?: \n";
    std::cin >> columns;

    std::cout << "Enter a symbol to use: \n";
    std::cin >> symbol;

    for(int i = 1; i <= rows; i++){
        for(int j = 1; j <= columns; j++){
            std::cout << symbol << ' ';
        }
        std::cout << '\n';
    }
}
*/

// 23.Random number generator

/*
#include <ctime>

int main()
{
    srand(time(NULL));

    int num1 = (rand() % 6) +1 ;
    int num2 = (rand() % 6) +1 ;
    int num3 = (rand() % 6) +1 ;

    std::cout << num1 << '\n';
    std::cout << num2 << '\n';
    std::cout << num3 << '\n';

    return 0;
}
*/

// 24.Random event generator
/*
#include <ctime>

int main()
{
    srand(time(0));
    int randNum = rand() % 5 + 1;

    // execute code depending on different case
    switch(randNum){
        case 1: std::cout << "You win a bumper sticker!\n";
            // exit out of switch after one statement
            break;
        case 2: std::cout << "You win a t-shirt!\n";
    }
}
*/

// 25. Number guessing game
// do-while loop
/*
do{
i++ can increment in do part
}while();
*/
/*
#include <ctime>

int main()
{
    int randNum = rand() % 100 + 1;

    std::cout << randNum << "\n";

    int num;
    int tries;

    do{
        std::cout << "Guess number!\n";
        std::cin >> num;
        if(num < randNum){
            std::cout << "Too low!\n";
        }
        else if(num > randNum){
            std::cout << "Too high!\n";
        }
        else{
            std::cout << "Correct! # of tries: " << tries << "\n";
        }
        tries++;
    }while (num != randNum); 

    std::cout << num << "\n";



    return 0;
}
*/

// 26.User defined functions

#include <iostream>

void happyBirthday(){
    std::cout << "HBDTY";
}

int main()
{

    happyBirthday();

    return 0;
}

// functions should be declared BEFORE the main function
// functions are usually defined AFTER the main function













