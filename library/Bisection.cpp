#include <iostream>
#include <cmath>
#include <iomanip>

double myfunc(double x)
{
    return 8 * std::pow(x, 3) - 4 * std::pow(x, 2) - 3 * x + 1;
}

double bisection(double a, double b, double epsilon)
{
    double fa{myfunc(a)};
    double fb{myfunc(b)};
    if (std::abs(a - b) < epsilon)
    {
        return a;
    }
    if (fa * fb < 0)
    {
        double x0 = 0.5 * (a + b);
        if (fa * myfunc(x0) < 0)
        {
            return bisection(a, x0, epsilon);
        }
        else if (fb * myfunc(x0) < 0)
        {
            return bisection(b, x0, epsilon);
        }
    }
    // else if ((fa * fb > 0))
    // {
    //     std::cout << "No root in this range.\n"
    // }
    // else if ((fa == 0))
    // {
    //     std::cout << "Root is at a = " << a << '\n';
    // }
    // else if ((fb == 0))
    // {
    //     std::cout << "Root is at b = " << b << '\n';
    // }
}

int main()
{
    // Error eps = 1e-12
    double epsilon{1e-10};
    double x0 = bisection(0, 0.5, epsilon);

    // In với độ chính xác cao
    std::cout
        << std::fixed << std::setprecision(15);
    std::cout << "Result: x0 = " << x0 << std::endl;

    return 0;
}