#include <iostream>
#include <cmath>
#include <iomanip>

double myfunc(double x)
{
    return (8 * std::pow(x, 3) - 4 * std::pow(x, 2) - 3 * x + 1);
}

double derivative(double x)
{
    return (24 * std::pow(x, 2) - 8 * x - 3);
}

double newtonRaphson(double x0, double epsilon, int max_iter = 1000)
{
    double x_n = x0;
    const double tiny = 1e-16; // ngưỡng để tránh chia cho 0
    for (int i = 0; i < max_iter; ++i)
    {
        double f_n = myfunc(x_n);
        double f_prime = derivative(x_n);

        if (std::abs(f_prime) < tiny)
        {
            // derivative too small: stop or fallback
            break;
        }

        double x_next = x_n - f_n / f_prime;

        // dừng theo bước (x) hoặc theo giá trị hàm (tuỳ bạn)
        if (std::abs(x_next - x_n) <= epsilon)
            return x_next;
        // or: if (std::abs(f_n) <= epsilon) return x_n;

        x_n = x_next;
    }
    return x_n; // trả x gần đúng sau max_iter
}

int main()
{
    // Error eps = 1e-12
    double epsilon{1e-10};
    double x0{0.5};
    double x_star{newtonRaphson(x0, epsilon)};

    // In với độ chính xác cao
    std::cout
        << std::fixed << std::setprecision(15);
    std::cout << "Result: x0 = " << x_star << std::endl;

    return 0;
}