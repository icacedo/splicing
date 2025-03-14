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