#pragma GCC optimize(2)

#ifndef BigInteger
#define BigInteger

#include <algorithm>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <vector>
#include <math.h>
using namespace std;

class Poly_div
{
public:
    struct complex
    {
        double a, b;
        complex() : a(0), b(0) {}
        complex(double x, double y) : a(x), b(y) {}
        void operator+=(const complex& z) { a += z.a, b += z.b; }
        complex operator-(const complex& z) { return { a - z.a, b - z.b }; }
        complex operator*(const complex& z) { return { a * z.a - b * z.b, a * z.b + b * z.a }; }
    };
    const long long Mod = 998244353;
    long long L;

    complex* W[2];
    int* a, * b, * bi, * x, * y, * rev;
    complex* t1, * t2;
    long long* t3;
    int n, m, d;

    // 构造函数
    Poly_div(long long num)
    {
        L = num;
        rev = new int[L << 1]();
        W[0] = new complex[L]();
        W[1] = new complex[L]();
        a = new int[L]();
        b = new int[L]();
        bi = new int[L]();
        x = new int[L]();
        y = new int[L]();
        t1 = new complex[L]();
        t2 = new complex[L]();
        t3 = new long long[L + 1]();
        n = 0;
        m = 0;
        d = 0;
    }
    ~Poly_div()
    {
        delete[] rev;
        delete[] W[0];
        delete[] W[1];
        delete[] a;
        delete[] b;
        delete[] bi;
        delete[] x;
        delete[] y;
        delete[] t1;
        delete[] t2;
        delete[] t3;
    }
    // 得到函数
    void get(std::vector<int>& A, std::vector<int>& B)
    { // A B都是规范化之后的
        n = A.size();
        for (int i = 0; i < A.size(); i++)
        {
            a[A.size() - 1 - i] = A[i];
        }
        m = B.size();
        for (int i = 0; i < B.size(); i++)
        {
            b[B.size() - 1 - i] = B[i];
        }
    }
    // 除法
    void div(std::vector<int>& A)
    {
        FFT_Init();
        d = divide(n, m);
        lack(n, m, d);
        for (int i = 0; i < n + m; i++)
        {
            A.push_back(x[i]);
        }
        return;
    }

    // 初始化
    void FFT_Init()
    {

        int i, j;
        const int l = L >> 1;
        const double pi = acos(-1);
        for (i = j = 1; j < L; j <<= 1)
            for (; i < j << 1; i++)
                rev[i << 1] = rev[i], rev[i << 1 | 1] = rev[i] + j;
        for (i = 0; i < l; i++)
            W[0][i + l] = { cos(pi * i / l), sin(pi * i / l) };
        for (i = l - 1; i; i--)
            W[0][i] = W[0][i << 1];

        memcpy(W[1], W[0], L * sizeof(complex));
        for (i = 1; i < L; i++)
            W[1][i].b *= -1;
    }
    // 过程
    void FFT(complex* f, int len, int sign)
    {
        int i = 1, j, k;
        complex* p, * q, * w, t;
        t.a = 0;
        t.b = 0;
        sign = (int)(sign < 0);
        for (int* r = rev + len + 1; i < len; i++, r++)
            if (i < (k = *r))
                t = f[i], f[i] = f[k], f[k] = t;

        for (i = 1; i < len; i <<= 1)
            for (j = 0; j < len; j += i, j += i)
                for (q = (p = f + j) + i, w = W[sign] + i, k = 0; k < i; k++, p++, q++)
                    t = *q * *w++, * q = *p - t, * p += t;

        if (sign)
        {
            double p = 1. / len;
            for (i = 0; i < len; i++, f++)
                f->a *= p, f->b *= p;
        }
    }

    void majutsu(int len) // 计算b的倒数，有效位数为len，结果存进bi里
    {

        bool g = 1;
        int l = 16, l2 = 32, l4 = 64, i;
        long double d = 0, e = 1;
        for (i = 0; i < 20; i++)
        {
            d = d + e * b[i], e *= 0.1;
        }
        d = 10. / d;

        if (d < 10)
        {
            for (i = 0; i <= l; i++)
                bi[i] = d, d = (d - bi[i]) * 10;
            bi[l - 1] += (bi[l] > 4);
        }
        else
            bi[0] = 10;

        while (l < len)
        {
        p:
            memset(t1, 0, sizeof(complex) * l4);
            memset(t2, 0, sizeof(complex) * l4);
            for (i = 0; i < l2; i++)
                t1[i].a = b[i];
            for (i = 0; i < l; i++)
                t2[i].a = bi[i];
            FFT(t1, l4, 1), FFT(t2, l4, 1);

            for (i = 0; i < l4; i++)
            {
                t1[i] = t1[i] * t2[i];
                t1[i].a = 20 - t1[i].a, t1[i].b = -t1[i].b;
                t2[i] = t1[i] * t2[i];
            }

            FFT(t2, l4, -1);
            t3[l4] = 0;

            for (i = l4 - 1; i >= 0; i--)
            {
                t3[i] = (long long)(floor(t2[i].a + 0.5)) + t3[i + 1] / 10, t3[i + 1] %= 10;
                if (t3[i + 1] < 0)
                    t3[i + 1] += 10, t3[i]--;
            }
            if (t3[0] > 9)
            {
                bi[0] = t3[0] / 10, t3[0] %= 10;
                for (i = 1; i < l2; i++)
                    bi[i] = t3[i - 1];
                bi[l2 - 1] += (t3[l2 - 1] > 4);
            }
            else
            {
                for (i = 0; i < l2; i++)
                    bi[i] = t3[i];
                bi[l2 - 1] += (t3[l2] > 4);
            }
            l <<= 1, l2 <<= 1, l4 <<= 1;
        }

        if (g)
        {
            g = 0, l >>= 1, l2 >>= 1, l4 >>= 1;
            goto p;
        } // 迭代过程结束后末几位仍会有偏差，于是再单独迭代一次
    }
    int divide(int n, int m) // 计算a/b，整数部分存进x里
    {
        int p = n - m + 16, l, i;
        majutsu(p);
        for (l = 1; l < n + p; l <<= 1)
            ;
        memset(t1, 0, sizeof(complex) * l);
        memset(t2, 0, sizeof(complex) * l);
        for (i = 0; i < n; i++)
            t1[i].a = a[i];
        for (i = 0; i < p; i++)
            t2[i].a = bi[i];
        FFT(t1, l, 1), FFT(t2, l, 1);
        for (i = 0; i < l; i++)
            t1[i] = t1[i] * t2[i];
        FFT(t1, l, -1);
        t3[l] = 0;
        for (i = l - 1; i >= 0; i--)
            t3[i] = (long long)(t1[i].a + 0.5) + t3[i + 1] / 10, t3[i + 1] %= 10;
        if (t3[0] > 9)
        {
            x[0] = t3[0] / 10, t3[0] %= 10, l = n - m + 1;
            for (i = 0; i < n - m; i++)
                x[i + 1] = t3[i];
        }
        else
            for (l = n - m, i = 0; i < n - m; i++)
                x[i] = t3[i];
        return l;
    }
    void lack(int n, int m, int& d) // 微调
    {
        int l, i, j;
        char t;
        long long tl = 0;
        for (i = 0, j = n - 1; i < j; i++, j--)
            t = a[i], a[i] = a[j], a[j] = t;
        for (i = 0, j = m - 1; i < j; i++, j--)
            t = b[i], b[i] = b[j], b[j] = t;
        for (i = 0, j = d - 1; i < j; i++, j--)
            t = x[i], x[i] = x[j], x[j] = t;
        for (l = 1; l < n; l <<= 1)
            ;
        memset(t1, 0, sizeof(complex) * l);
        memset(t2, 0, sizeof(complex) * l);
        for (i = 0; i < m; i++)
            t1[i].a = b[i];
        for (i = 0; i < d; i++)
            t2[i].a = x[i];
        t2[0].a += 1.;
        FFT(t1, l, 1), FFT(t2, l, 1);
        for (i = 0; i < l; i++)
            t1[i] = t1[i] * t2[i];
        FFT(t1, l, -1);
        for (i = 0; i <= n; i++)
            tl += (long long)(t1[i].a + 0.5), y[i] = tl % 10, tl /= 10;
        for (i = n; i >= 0; i--)
        {
            if (y[i] > a[i])
                return;
            else if (y[i] < a[i])
                break;
        }
        for (x[0]++, i = 0; x[i] > 9; i++)
            x[i + 1] += x[i] / 10, x[i] %= 10;
        if (x[d])
            d++;
    }
};

class Polynomial
{
public:
    long long Mod = 998244353;
    int G = 3, iG = 332748118;
    int MS;
    long long* Inv; // 逆元
    long long Sz, * R;
    long long InvSz; //	NTT
    int N, M;
    long long* A1, * B1;

    Polynomial(int sum)
    {
        MS = std::max(64, sum * 2);
        Inv = new long long[MS];
        R = new long long[MS];
        A1 = new long long[MS];
        B1 = new long long[MS];

        for (int i = 0; i < MS; i++)
        {
            Inv[i] = 0;
            R[i] = 0;
            A1[i] = 0;
            B1[i] = 0;
        }
        InvSz = 0;
        Sz = 0;
        N = 0;
        M = 0;
        Init(MS);
    }

    void get(std::vector<int>& A, std::vector<int>& B);
    std::vector<int>& add();
    std::vector<int>& sub();
    std::vector<int>& mul();
    static int cmp(std::vector<int>& A, std::vector<int>& B);

    void Clear();
    long long qPow(long long b, int e);
    void Init(int N);
    void InitFNTT(int N);
    void FNTT(long long* A, int Ty);
    void PolyMul(long long* A, long long* B, int deg);
    void PolyInv(long long* A, int N, long long* B);
    void PolyLn(long long* A, int N, long long* B);
    void PolyExp(long long* A, int N, long long* B);
};

class Integer
{

protected:
    std::vector<int> data;
    int flag = 1;

public:
    Integer();
    Integer(long long a);
    Integer(int flag, std::vector<int>& a);

    void check();
    friend std::istream& operator>>(std::istream& in, Integer& a);
    friend std::ostream& operator<<(std::ostream& out, Integer& a);
    long long abs(long long a);
    Integer& operator=(Integer& a);
    Integer& operator=(long long a);
    Integer& operator=(std::string& a);
    int operator==(long long a);
    int operator==(Integer& a);
    int operator!=(long long a);
    int operator!=(Integer& a);
    int operator<(Integer& a);
    int operator<=(Integer& a);
    int operator>(Integer& a);
    int operator>=(Integer& a);
    Integer& operator+(Integer& a);
    Integer& operator-(Integer& a);
    Integer& operator-();
    Integer& operator/(Integer& a);
    Integer& operator%(Integer& a);
    Integer& operator*(Integer& a);
    std::vector<int>& Polymul_tool(std::vector<int> a, std::vector<int> b);
    int size() { return data.size(); }
    void clear()
    {
        flag = 1;
        data.clear();
        data.push_back(0);
    }
};
long long Integer::abs(long long a)
{
    if (a < 0)
        return -a;
    else
        return a;
}

Integer::Integer()
{
    flag = 1;
    data.push_back(0);
}

Integer::Integer(long long a)
{
    if (a < 0)
        flag = -1;
    else
        flag = 1;
    a = abs(a);
    while (a)
    {
        data.push_back(a % 10);
        a /= 10;
    }
}

Integer::Integer(int flag, std::vector<int>& a)
{
    data = a;
    this->flag = 1;
}

int Integer::operator<(Integer& a)
{
    if (flag == 1) {
        if (a.flag == 1) {
            if (Polynomial::cmp(data, a.data) < 0)
                return 1;
            else
                return 0;
        }
        else
        {
            return 0;
        }
    }
    else {
        if (a.flag == 1) {
            return 1;
        }
        else {
            if (Polynomial::cmp(data, a.data) > 0)
                return 1;
            else
                return 0;
        }
    }
    return 0;
}

int Integer::operator<=(Integer& a)
{
    if (flag == 1) {
        if (a.flag == 1) {
            if (Polynomial::cmp(data, a.data) > 0)
                return 0;
            else
                return 1;
        }
        else
        {
            return 0;
        }
    }
    else {
        if (a.flag == 1) {
            return 1;
        }
        else {
            if (Polynomial::cmp(data, a.data) < 0)
                return 0;
            else
                return 1;
        }
    }
    return 0;
}

int Integer::operator>(Integer& a)
{
    if (flag == 1) {
        if (a.flag == 1) {
            if (Polynomial::cmp(data, a.data) > 0)
                return 1;
            else
                return 0;
        }
        else
        {
            return 1;
        }
    }
    else {
        if (a.flag == 1) {
            return 0;
        }
        else {
            if (Polynomial::cmp(data, a.data) < 0)
                return 1;
            else
                return 0;
        }
    }
    return 0;
}

int Integer::operator>=(Integer& a)
{
    if (flag == 1) {
        if (a.flag == 1) {
            if (Polynomial::cmp(data, a.data) < 0)
                return 0;
            else
                return 1;
        }
        else
        {
            return 1;
        }
    }
    else {
        if (a.flag == 1) {
            return 0;
        }
        else {
            if (Polynomial::cmp(data, a.data) > 0)
                return 0;
            else
                return 1;
        }
    }
    return 0;
}

void Integer::check()
{
    for (int i = data.size() - 1; i >= 1; i--)
        if (data[i] == 0)
            data.pop_back();
        else
            break;

    if (data.size() == 1 && data[0] == 0)
        flag = 1;
    if (data.size() == 0)
    {
        data.push_back(0);
        flag = 1;
    }
}

int Integer::operator==(long long a)
{
    Integer tool(a);
    return (*this) == tool;
}

int Integer::operator==(Integer& a)
{
    if (Polynomial::cmp(data, a.data) == 0)
        return 1;
    else
        return 0;
}

int Integer::operator!=(long long a)
{
    if ((*this) == a)
        return 0;
    else
        return 1;
}

int Integer::operator!=(Integer& a)
{
    if ((*this) == a)
        return 0;
    else
        return 1;
}

Integer& Integer::operator=(long long a)
{
    if (a < 0)
        flag = -1;
    else
        flag = 1;
    a = abs(a);
    data.clear();
    while (a)
    {
        data.push_back(a % 10);
        a /= 10;
    }
    return (*this);
}

Integer& Integer::operator=(std::string& date)
{

    if (date[0] == '-')
    {
        flag = -1;
        for (int i = date.size() - 1; i >= 1; i--)
        {
            data.push_back(date[i] - 48);
        }
    }
    else
    {
        flag = 1;
        data.clear();
        for (int i = date.size() - 1; i >= 0; i--)
        {
            data.push_back(date[i] - 48);
        }
    }
    check();
    return *this;
}

Integer& Integer::operator=(Integer& a)
{
    data = a.data;
    flag = a.flag;
    return *this;
}

Integer& Integer::operator+(Integer& a)
{
    if (flag == 1 && a.flag == 1)
    {
        Integer* p;
        p = new Integer;
        p->flag = 1;
        Polynomial poly(data.size() + a.data.size() - 1);
        poly.get(data, a.data);
        p->data = poly.add();
        p->check();
        return (*p);
    }
    else if (flag == -1 && a.flag == -1)
    {
        Integer* p;
        p = new Integer;
        p->flag = -1;
        Polynomial poly(data.size() + a.data.size() - 1);
        poly.get(data, a.data);
        p->data = poly.add();
        p->check();
        return (*p);
    }
    else if (flag == 1 && a.flag == -1)
    {
        Integer b = a;
        b.flag = 1;
        return (*this) - b;
    }
    else if (flag == -1 && a.flag == 1)
    {
        Integer b = (*this);
        b.flag = 1;
        return a - b;
    }
    return a;
}

Integer& Integer::operator-(Integer& a)
{
    if (flag == 1 && a.flag == 1)
    {
        Integer* p;
        p = new Integer;
        Polynomial poly(data.size() + a.data.size() - 1);
        p->flag = poly.cmp(this->data, a.data);
        if (p->flag == 1)
        {
            poly.get(data, a.data);
            p->data = poly.sub();
            return (*p);
        }
        else
        {
            poly.get(a.data, data);
            p->data = poly.sub();
            return (*p);
        }
    }
    else if (flag == -1 && a.flag == -1)
    {
        Integer* p;
        p = new Integer;
        Polynomial poly(data.size() + a.data.size() - 1);
        p->flag = poly.cmp(this->data, a.data);
        if (p->flag == 1)
        {
            poly.get(data, a.data);
            p->data = poly.sub();
            p->flag *= -1;
            return (*p);
        }
        else
        {
            poly.get(a.data, data);
            p->data = poly.sub();
            p->flag *= -1;
            return (*p);
        }
    }
    else
    {
        Integer b = a;
        b.flag *= -1;
        return (*this) + b;
    }
}

Integer& Integer::operator-()
{
    Integer* p;
    p = new Integer;
    p->data.clear();
    p->flag *= -1;
    p->data = data;
    return *p;
}

std::vector<int>& Integer::Polymul_tool(std::vector<int> a, std::vector<int> b)
{

    int deg = a.size() + b.size() - 1;

    for (int i = a.size() - 1; i >= 0; i--)
        if (a[i] == 0)
            a.pop_back();
        else
            break;
    for (int i = b.size() - 1; i >= 0; i--)
        if (b[i] == 0)
            b.pop_back();
        else
            break;

    Polynomial test(a.size() + b.size());
    test.get(a, b);
    std::vector<int>* c;
    c = &test.mul();

    return *c;
}

Integer& Integer::operator*(Integer& a)
{
    Integer* p;
    p = new Integer;
    p->data.clear();
    p->flag = this->flag * a.flag;

    p->data = Polymul_tool(this->data, a.data);

    return *p;
}

Integer& Integer::operator/(Integer& a)
{

    Integer* p;
    p = new Integer;
    p->data.clear();
    p->flag = this->flag * a.flag;

    int cmp_temp = Polynomial::cmp(data,a.data);
    if (cmp_temp == 0) {
        (*p)=1;
        return (*p);
    }
    else if (cmp_temp == -1) {
        return (*p);
    }

    int len = 64;
    int len_2 = 10 + this->data.size() + a.data.size();

    while (len < len_2)
        len *= 2;
    Poly_div test(len * 2);

    test.get(this->data, a.data);
    test.div(p->data);
    p->check();
    return (*p);
}

Integer& Integer::operator%(Integer& a)
{
    return (*this) - (*this) / a * a;
}

std::istream& operator>>(std::istream& in, Integer& a)
{
    a.data.clear();
    std::string date;
    std::cin >> date;
    if (date[0] == '-')
    {
        a.flag = -1;
        for (int i = date.size() - 1; i >= 1; i--)
        {
            a.data.push_back(date[i] - 48);
        }
    }
    else
    {
        a.flag = 1;
        a.data.clear();
        for (int i = date.size() - 1; i >= 0; i--)
        {
            a.data.push_back(date[i] - 48);
        }
    }
    a.check();
    return in;
}

std::ostream& operator<<(std::ostream& out, Integer& a)
{
    a.check();
    if (a.data.size() == 0)
        out << '0';
    else if (a.flag == 1)
    {
        for (int i = a.data.size() - 1; i >= 0; i--)
            out << a.data[i];
    }
    else
    {
        out << "-";
        for (int i = a.data.size() - 1; i >= 0; i--)
            out << a.data[i];
    }
    return out;
}

// Polynomial的函数
void Polynomial::get(std::vector<int>& A, std::vector<int>& B)
{ // A B都是规范化之后的
    N = A.size();
    for (int i = A.size() - 1; i >= 0; i--)
    {
        A1[i] = A[i];
    }
    M = B.size();
    for (int i = B.size() - 1; i >= 0; i--)
    {
        B1[i] = B[i];
    }
}

int Polynomial::cmp(std::vector<int>& A, std::vector<int>& B)
{

    if (A.size() > B.size())
        return 1;
    else if (A.size() < B.size())
        return -1;
    else
    {
        for (int i = A.size() - 1; i >= 0; i--)
        {
            if (A[i] > B[i])
                return 1;
            else if (A[i] < B[i])
                return -1;
        }
        return 0;
    }
    return 0;
}

std::vector<int>& Polynomial::mul()
{
    std::vector<int>* P;
    P = new std::vector<int>;
    if (N == 0 || M == 0)
    {
        P->push_back(0);
        return (*P);
    }

    PolyMul(A1, B1, N + M - 2);
    long long tool = 0;
    for (int i = 0; i < N + M - 1; i++)
    {
        long long a = tool + (A1[i] + Mod) % Mod;
        P->push_back(a % 10);
        tool = a / 10;
    }
    while (tool)
    {
        P->push_back(tool % 10);
        tool /= 10;
    }
    for (int i = P->size() - 1; i >= 1; i--)
    {
        if ((*P)[i] == 0)
            P->pop_back();
        else
            break;
    }
    Clear();
    return (*P);
}

std::vector<int>& Polynomial::sub()
{
    std::vector<int>* P;
    P = new std::vector<int>;
    if (N == 0 && M == 0)
    {
        P->push_back(0);
        return (*P);
    }
    int tool1 = 0;
    int Z = std::max(N, M);
    for (int i = 0; i < Z; i++)
    {
        int tool2 = A1[i] - B1[i] + tool1;
        if (tool2 >= 0)
        {
            tool1 = 0;
            P->push_back(tool2);
        }
        else
        {
            tool1 = -1;
            P->push_back(10 + tool2);
        }
    }
    for (int i = P->size() - 1; i >= 1; i--)
    {
        if ((*P)[i] == 0)
            P->pop_back();
        else
            break;
    }
    Clear();
    return (*P);
}

std::vector<int>& Polynomial::add()
{
    std::vector<int>* P;
    P = new std::vector<int>;
    if (N == 0 && M == 0)
    {
        P->push_back(0);
        return (*P);
    }
    int tool1 = 0;
    int Z = std::max(N, M);
    for (int i = 0; i < Z; i++)
    {
        int tool2 = A1[i] + B1[i] + tool1;
        P->push_back(tool2 % 10);
        tool1 = tool2 / 10;
    }

    if (tool1 > 0)
    {
        P->push_back(tool1);
    }
    for (int i = P->size() - 1; i >= 1; i--)
    {
        if ((*P)[i] == 0)
            P->pop_back();
        else
            break;
    }
    Clear();
    return (*P);
}

void Polynomial::Clear()
{
    for (int i = 0; i < MS; i++)
    {
        R[i] = 0;
        A1[i] = 0;
        B1[i] = 0;
    }
}

long long Polynomial::qPow(long long b, int e)
{ // 快速幂
    long long a = 1;
    for (; e; e >>= 1, b = b * b % Mod)
        if (e & 1)
            a = a * b % Mod;
    return a;
}

void Polynomial::Init(int N)
{ // 求逆元
    Inv[1] = 1;
    for (int i = 2; i < N; ++i)
        Inv[i] = -(Mod / i) * Inv[Mod % i] % Mod;
    return;
}

void Polynomial::InitFNTT(int N)
{
    int Bt = 0;
    for (; 1 << Bt <= N; ++Bt)
        ;
    if (Sz == (1 << Bt))
        return;
    Sz = 1 << Bt;
    InvSz = -(Mod - 1) / Sz;
    for (int i = 1; i < Sz; ++i)
        R[i] = R[i >> 1] >> 1 | (i & 1) << (Bt - 1);
}

void Polynomial::FNTT(long long* A, int Ty)
{
    for (int i = 0; i < Sz; ++i)
        if (R[i] < i)
            std::swap(A[R[i]], A[i]);
    for (int j = 1, j2 = 2; j < Sz; j <<= 1, j2 <<= 1)
    {
        long long gn = qPow(~Ty ? G : iG, (Mod - 1) / j2), g, X, Y;
        for (int i = 0, k; i < Sz; i += j2)
        {
            for (k = 0, g = 1; k < j; ++k, g = g * gn % Mod)
            {
                X = A[i + k], Y = g * A[i + j + k] % Mod;
                A[i + k] = (X + Y) % Mod, A[i + j + k] = (X - Y) % Mod;
            }
        }
    }
    if (!~Ty)
        for (int i = 0; i < Sz; ++i)
            A[i] = A[i] * InvSz % Mod;
}

void Polynomial::PolyMul(long long* A, long long* B, int deg)
{ // 多项式乘法  A=A*B
    InitFNTT(deg);
    FNTT(A, 1);
    FNTT(B, 1);
    for (int i = 0; i < Sz; i++)
    {
        A[i] = (A[i] * B[i]) % Mod;
    }
    FNTT(A, -1);
    FNTT(B, -1);
}

void Polynomial::PolyInv(long long* A, int N, long long* B)
{ // 多项式求inv 对x^N取mod  A*B=1
    // long long tA[MS], tB[MS];
    long long* tA = new long long[MS];
    long long* tB = new long long[MS];
    B[0] = qPow(A[0], Mod - 2);
    for (int L = 1; L < N; L <<= 1)
    {
        int L2 = L << 1, L4 = L << 2;
        InitFNTT(L4);
        memcpy(tA, A, 8 * L2);
        memset(tA + L2, 0, 8 * (Sz - L2));
        memcpy(tB, B, 8 * L);
        memset(tB + L, 0, 8 * (Sz - L));
        FNTT(tA, 1), FNTT(tB, 1);
        for (int i = 0; i < Sz; ++i)
            tB[i] = (2 - tB[i] * tA[i]) % Mod * tB[i] % Mod;
        FNTT(tB, -1);
        for (int i = 0; i < L2; ++i)
            B[i] = tB[i];
    }
    delete[] tA;
    delete[] tB;
}

void Polynomial::PolyLn(long long* A, int N, long long* B)
{ // 多项式求ln 对x^N取mod A[0]=1 B=ln(A)
    long long* tA = new long long[MS];
    long long* tB = new long long[MS];
    PolyInv(A, N - 1, tB);
    InitFNTT(N * 2 - 3);
    for (int i = 1; i < N; ++i)
        tA[i - 1] = i * A[i] % Mod;
    memset(tA + N - 1, 0, 8 * (Sz - N + 1));
    memset(tB + N - 1, 0, 8 * (Sz - N + 1));
    FNTT(tA, 1), FNTT(tB, 1);
    for (int i = 0; i < Sz; ++i)
        tA[i] = (long long)tA[i] * tB[i] % Mod;
    FNTT(tA, -1);
    B[0] = 0;

    for (int i = 1; i < N; ++i)
        B[i] = (long long)tA[i - 1] * Inv[i] % Mod;

    delete[] tA;
    delete[] tB;
}

void Polynomial::PolyExp(long long* A, int N, long long* B)
{ // 多项式求exp 对x^N取mod A[0]=0  B=e^A
    long long* tA = new long long[MS];
    long long* tB = new long long[MS];
    B[0] = 1;
    for (int L = 1; L < N; L <<= 1)
    {
        int L2 = L << 1, L4 = L << 2;
        memset(B + L, 0, 8 * (L2 - L));
        PolyLn(B, L2, tB);
        InitFNTT(L4);
        memcpy(tA, B, 8 * L);
        memset(tA + L, 0, 8 * (Sz - L));
        for (int i = 0; i < L2; ++i)
            tB[i] = ((!i) - tB[i] + A[i]) % Mod;
        memset(tB + L2, 0, 8 * (Sz - L2));
        FNTT(tA, 1), FNTT(tB, 1);
        for (int i = 0; i < Sz; ++i)
            tA[i] = tA[i] * tB[i] % Mod;
        FNTT(tA, -1);
        for (int i = 0; i < L2; ++i)
            B[i] = tA[i];
    }
    delete[] tA;
    delete[] tB;
}

#endif // !Integer


int main()
{
    Integer a, b,c;
    cin >> a >> b;
    c=a / b;
    cout << c << endl;
    cout<<a-b*c<<endl;
    return 0;
}
