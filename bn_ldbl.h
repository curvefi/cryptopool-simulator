#pragma once
#define _CRT_SECURE_NO_WARNINGS 1

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cassert>
#include <memory.h>
#include <ctime>
#include <cmath>
#include <map>
#include <set>
#include <tuple>

using namespace std;
#define NDEBUG 1

// Файл bn.h
class BN {
public:
    BN(BN const &orig) {
        _b = orig._b;
    }

    explicit BN(long long init_int) {
        _b = init_int;
    }

    explicit BN(double init_double) {
        _b = init_double;
    }

    BN() {
        _b = 0;
    }


    explicit BN(string const &in) {
        _b = atof(in.c_str());
    }

    BN(const char *init_string) {
        _b = atof(init_string);
    }


    ~BN() {
    }

    bool is_zero() const {
        return _b == 0;
    }


    BN &abs() { // Взять модуль
        _b = fabs(_b);
        return *this;
    }

    // Вернуть <0, если отрицательное, >0, если положительное, 0.
    int sign() const {
        if (_b < 0) return -1;
        if (_b == 0) return 0;
        return 1;
    }

    // Выдать представление BN в системе счисления radix в виде строки
    string to_string(int radix = 10) const {
        char buf[100];
        sprintf(buf, "%.0Lf", _b);
        return buf;
    }

    // Операции, аналогичные +=, -=, *=, /=, %=
    BN &add_to(BN const &r) {
        _b += r._b;
        return *this;
    }

    BN &sub_to(BN const &r) {
        _b -= r._b;
        return *this;
    }

    BN &mul_to(BN const &r) {
        _b *= r._b;
        return *this;
    }

    BN &div_to(BN const &r) {
        _b /= r._b;
        return *this;
    }


    // Извлечь корень степени reciprocal из BN (бонусная функция)
    BN &root_to(int root = 2) {
        if (root == 2) {
            _b = sqrt(_b);
        }
        return *this;
    }

    BN &operator=(BN const &orig) { // Создать копию существующего BN
        _b = orig._b;
        return *this;
    }


    bool operator<(BN const &right) const {
        return _b < right._b;
    }

    bool operator>(BN const &right) const {
        return _b > right._b;
    }

    bool operator<=(BN const &right) const {
        return _b <= right._b;
    }

    bool operator>=(BN const &right) const {
        return _b >= right._b;
    }

    bool operator==(BN const &right) const {
        return _b == right._b;
    }

    bool operator!=(BN const &right) const {
        return _b != right._b;
    }

    BN &operator+=(const BN &oth) {
        add_to(oth);
        return *this;
    }

    BN &operator*=(const BN &oth) {
        mul_to(oth);
        return *this;
    }

    BN &operator-=(const BN &oth) {
        sub_to(oth);
        return *this;
    }

    BN &operator/=(const BN &oth) {
        div_to(oth);
        return *this;
    }

    BN operator+(const BN &r) const {
        BN ret(*this);
        return ret += r;
    }

    BN operator-(const BN &r) const {
        BN ret(*this);
        return ret -= r;
    }

    BN operator*(const BN &r) const {
        BN ret(*this);
        return ret *= r;
    }

    BN operator/(const BN &r) const {
        BN ret(*this);
        return ret /= r;
    }

    double get_double() const { // Veeeeeeeryyyy ugly
        return _b;
    }

private:
    long double _b = 0;
};



