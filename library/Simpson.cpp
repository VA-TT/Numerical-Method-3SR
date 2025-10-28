#include <iostream>
#include <cmath>
#include <iomanip>
#include <functional>

double myfunc(double x)
{
    return (std::pow(x, 7));
}

double simpson_quad(std::function<double(double)> f, double a, double b, int n)
{
    if (n % 2 != 0)
        throw std::invalid_argument("n must be even");
    double delta = (b - a) / n;
    double S{0.0};
    double fa{}, fb{}, fm{};
    for (int i{1}; i <= n; ++i)
    {
        fa = f(a + (i - 1) * delta);
        fb = f(a + i * delta);
        fm = f(a + i * delta - delta / 2);
        S += fa + fb + 4 * fm;
    }
    return S * delta / 6;
}

int main()
{

    // Tich phan xac dinh theo phuong phap Simpson cua ham x^4 trong doan [0,1], chia nho 10 lan
    std::cout << std::fixed << std::setprecision(15);
    double I = simpson_quad(myfunc, 0, 2, 10);
    std::cout << "Result: x0 = " << I << std::endl;

    return 0;
}