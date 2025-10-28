#pragma once
#include <iostream>
#include <numeric> // for std::gcd
#include <cmath>

class Fraction
{
private:
    int numerator;   // Tử số
    int denominator; // Mẫu số

    void simplify()
    {
        if (denominator < 0)
        {
            numerator = -numerator;
            denominator = -denominator;
        }
        int gcd = std::gcd(std::abs(numerator), std::abs(denominator));
        numerator /= gcd;
        denominator /= gcd;
    }

public:
    Fraction(int num = 0, int den = 1) : numerator(num), denominator(den)
    {
        if (denominator == 0)
            throw std::runtime_error("Denominator cannot be zero");
        simplify();
    }

    Fraction operator+(const Fraction &other) const
    {
        return Fraction(numerator * other.denominator + other.numerator * denominator,
                        denominator * other.denominator);
    }

    Fraction operator-(const Fraction &other) const
    {
        return Fraction(numerator * other.denominator - other.numerator * denominator,
                        denominator * other.denominator);
    }

    Fraction operator*(const Fraction &other) const
    {
        return Fraction(numerator * other.numerator,
                        denominator * other.denominator);
    }

    Fraction operator/(const Fraction &other) const
    {
        if (other.numerator == 0)
            throw std::runtime_error("Division by zero");
        return Fraction(numerator * other.denominator,
                        denominator * other.numerator);
    }

    Fraction operator-() const
    {
        return Fraction(-numerator, denominator);
    }

    bool operator==(const Fraction &other) const
    {
        return numerator * other.denominator == other.numerator * denominator;
    }

    bool operator!=(const Fraction &other) const
    {
        return !(*this == other);
    }

    friend std::ostream &operator<<(std::ostream &out, const Fraction &f)
    {
        if (f.denominator == 1)
        {
            out << f.numerator;
        }
        else
        {
            out << f.numerator << "/" << f.denominator;
        }
        return out;
    }
};