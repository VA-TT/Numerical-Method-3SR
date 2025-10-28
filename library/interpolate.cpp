#include <iostream>
#include <array>
#include <cmath>
#include <iomanip>

// Hàm tính đa thức cơ sở Lagrange l_i(x)
template <std::size_t N>
double basisPolynomial(std::size_t i,
                       const std::array<double, N>& data_x,
                       double x)
{
    double li = 1.0;
    double xi = data_x[i];

    for (std::size_t j = 0; j < N; ++j)
    {
        if (j != i)
        {
            double xj = data_x[j];
            li *= (x - xj) / (xi - xj);
        }
    }
    return li;
}

// Hàm nội suy Lagrange
template <std::size_t N>
double interpolatePolynomial(const std::array<double, N>& data_x,
                             const std::array<double, N>& data_y,
                             double x)
{
    double fx = 0.0;
    for (std::size_t i = 0; i < N; ++i)
    {
        fx += data_y[i] * basisPolynomial(i, data_x, x);
    }
    return fx;
}

int main()
{
    constexpr std::size_t numberPoints = 7;

    std::array<double, numberPoints> data_x
        {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    std::array<double, numberPoints> data_y
        {0.0, 0.8415, 0.9093, 0.1411, -0.7568, -0.9589, -0.2794};

    std::cout << std::fixed << std::setprecision(6);

    // Tính giá trị nội suy tại các điểm cách nhau 0.1
    for (double x = 0.0; x <= 6.0; x += 0.1)
    {
        double fx = interpolatePolynomial(data_x, data_y, x);
        std::cout << "P(" << x << ") = " << fx << '\n';
    }

    return 0;
}
