#include <string>
#include <stdint.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <chrono>
using namespace std;

string p;
const uint64_t BASE = 0x100000000;
uint32_t P[128];
uint32_t R[128] = {0};
int P_bits = 0;
int P_words = 0;
int R_bits = 0;
int R_words = 0;
uint32_t R2[128] = {0};
uint32_t P_[128] = {0};
uint32_t ZERO[128] = {0};
uint32_t ONE[128] = {1};
uint32_t TWO[128] = {2};

void mod_mul(uint32_t res[128], uint32_t a[128], uint32_t b[128]);
void mod_mul_trivial(uint32_t res[128], uint32_t a[128], uint32_t b[128]);
void mod_mul_mont(uint32_t res[128], uint32_t a[128], uint32_t b[128]);
#ifndef mod_mul
// #define mod_mul (mod_mul_mont)
#define mod_mul (mod_mul_trivial)
#endif

void mod_pow(uint32_t res[128], uint32_t a[128], uint32_t b[128]);
void mod_pow_trivial(uint32_t res[128], uint32_t a[128], uint32_t b[128]);
void mod_pow_mont(uint32_t res[128], uint32_t a[128], uint32_t b[128]);
#ifndef mod_pow
// #define mod_pow (mod_pow_mont)
#define mod_pow (mod_pow_trivial)
#endif

uint32_t POW2[32] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
                     1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288,
                     1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912,
                     1073741824, 2147483648};


int getBits(const uint32_t a[128])
{
    for (int i = 127; i >= 0; i--)
    {
        if (a[i])
            for (int j = 31; j >= 0; j--)
                if (a[i] & POW2[j])
                    return i * 32 + j;
    }

    return 0;
}

void remove_leading_zeros(string &s)
{
    size_t end = s.find_last_not_of('0');
    if (end != string::npos)
        s.erase(end + 1);
    else
        s = "0";
}

bool bigger(uint32_t a[128], uint32_t b[128])
{
    for (int i = 127; i >= 0; i--)
    {
        if (a[i] > b[i])
            return true;
        else if (a[i] < b[i])
            return false;
    }
    return false;
}

bool equal(uint32_t a[128], uint32_t b[128])
{
    for (int i = 0; i < 128; i++)
        if (a[i] != b[i])
            return false;
    return true;
}

void str_div2(string &s)
{
    int carry = 0;
    for (int i = s.length() - 1; i >= 0; i--)
    {
        int x = s[i] - '0';
        s[i] = (x + carry * 10) / 2 + '0';
        carry = x % 2;
    }

    remove_leading_zeros(s);
}

void str_mul2(string &s)
{
    int carry = 0;
    for (std::string::size_type i = 0; i < s.length(); i++)
    {
        int x = s[i] - '0';
        s[i] = (x * 2 + carry) % 10 + '0';
        carry = (x * 2 + carry) / 10;
    }

    if (carry)
        s = s + "1";
}

void str_add1(string &s)
{
    int carry = 1;
    for (std::string::size_type i = 0; i < s.length(); i++)
    {
        int x = s[i] - '0';
        s[i] = (x + carry) % 10 + '0';
        carry = (x + carry) / 10;
    }
    if (carry)
        s = s + "1";
}

void str2bi(uint32_t res[128], string &s)
{
    int count = 0;
    while (s != "0")
    {
        if ((s[0] - '0') & 1)
            res[count / 32] += POW2[count % 32];
        str_div2(s);
        count++;
    }
}

string bi2str(uint32_t res[128])
{
    string s = "0";
    for (int i = 127; i >= 0; i--)
        for (int j = 31; j >= 0; j--)
        {
            str_mul2(s);
            if (res[i] & POW2[j])
                str_add1(s);
        }

    reverse(s.begin(), s.end());
    return s;
}

void add(uint32_t res[128], uint32_t a[128], uint32_t b[128])
{
    uint64_t carry = 0;
    uint32_t temp[128] = {0};
    for (int i = 0; i < 128; i++)
    {
        uint64_t sum = carry + a[i] + b[i];
        temp[i] = sum & 0xffffffff;
        carry = sum >> 32;
    }

    for (int i = 0; i < 128; i++)
        res[i] = temp[i];
}

void sub(uint32_t res[128], uint32_t a[128], uint32_t b[128])
{
    uint64_t borrow = 0;
    uint32_t temp[128] = {0};
    for (int i = 0; i < 128; i++)
    {
        uint64_t temp1 = b[i] + borrow;
        if (a[i] < temp1)
        {
            temp[i] = BASE + a[i] - temp1;
            borrow = 1;
        }
        else
        {
            temp[i] = a[i] - temp1;
            borrow = 0;
        }
    }

    for (int i = 0; i < 128; i++)
        res[i] = temp[i];
}

void mul(uint32_t res[128], uint32_t a[128], uint32_t b[128])
{
    uint32_t temp[256] = {0};

    for (int i = 0; i < 128; i++)
    {
        uint64_t carry = 0;
        for (int j = 0; j < 128; j++)
        {
            uint64_t sum = (uint64_t)a[i] * b[j] + temp[i + j] + carry;
            temp[i + j] = sum & 0xffffffff;
            carry = sum >> 32;
        }
        if (carry)
            temp[i + 128] += carry;
    }

    for (int i = 0; i < 128; i++)
        res[i] = temp[i];
}

void div(uint32_t res[128], uint32_t a[128], uint32_t b[128])
{
    uint32_t temp[128] = {0};
    int n = getBits(b) / 32;
    int m = getBits(a) / 32 - n; // iterate times

    uint32_t d[128] = {0};
    d[0] = BASE / (b[n] + (uint64_t)1);
    uint32_t u_[129] = {0}, v_[128] = {0};
    uint64_t carry = 0;
    // scale a
    for (int i = 0; i < 128; i++)
    {
        uint64_t temp = (uint64_t)a[i] * d[0] + carry;
        u_[i] = temp & 0xffffffff;
        carry = temp >> 32;
    }
    if (carry)
        u_[128] = carry;
    // scale b
    mul(v_, b, d);

    int j = m;
    while (j >= 0)
    {
        uint32_t tem[129] = {0};
        for (int i = 0; i <= n + 1; i++)
            tem[i] = u_[i + j];
        uint64_t tem2 = (tem[n + 1] * BASE + tem[n]) / v_[n];
        tem2 = min(tem2, BASE - 1);

        uint32_t q_hat = static_cast<uint32_t>(tem2 & 0xffffffff);
        uint32_t qv[128] = {0};
        for (int i = 0; i < 128; i++)
        {
            uint64_t temp = (uint64_t)v_[i] * q_hat + carry;
            qv[i] = temp & 0xffffffff;
            carry = temp >> 32;
        }

        while (bigger(qv, tem))
        {
            q_hat--;
            sub(qv, qv, v_);
        }

        sub(tem, tem, qv);
        for (int i = 0; i <= n + 1; i++)
            u_[i + j] = tem[i];

        temp[j] = q_hat;

        j--;
    }

    for (int i = 0; i < 128; i++)
        res[i] = temp[i];
}

void mod(uint32_t res[128], uint32_t a[128])
{
    uint32_t temp[128] = {0};
    div(temp, a, P);
    mul(temp, temp, P);
    sub(res, a, temp);
}

void mod_add(uint32_t res[128], uint32_t a[128], uint32_t b[128])
{
    add(res, a, b);
    mod(res, res);
}

void mod_sub(uint32_t res[128], uint32_t a[128], uint32_t b[128])
{
    if (bigger(a, b))
    {
        sub(res, a, b);
        mod(res, res);
    }
    else
    {
        sub(res, b, a);
        mod(res, res);
        sub(res, P, res);
    }
}

void mod_mul_trivial(uint32_t res[128], uint32_t a[128], uint32_t b[128])
{
    mul(res, a, b);
    mod(res, res);
}

void mod_div(uint32_t res[128], uint32_t a[128], uint32_t b[128])
{
    div(res, a, b);
    mod(res, res);
}

void mod_pow_trivial(uint32_t res[128], uint32_t a[128], uint32_t b[128])
{
    uint32_t temp[128] = {1};
    uint32_t base[128] = {0};
    for (int i = 0; i < 128; i++)
        base[i] = a[i];

    for (int i = 0; i < 128; i++)
    {
        if (bigger(temp, ONE) && b[i] == 0)
            break;
        for (int j = 0; j < 32; j++)
        {
            if (b[i] & POW2[j])
                mod_mul(temp, temp, base);

            mod_mul(base, base, base);
        }
    }
    for (int i = 0; i < 128; i++)
        res[i] = temp[i];
}

void mod_R(uint32_t res[128], uint32_t a[128])
{
    uint32_t temp[128] = {0};
    for (int i = 0; i < R_bits / 32; i++)
        temp[i] = a[i];
    for (int i = 0; i < R_bits % 32; i++)
        if (a[R_bits / 32] & POW2[i])
            temp[R_bits / 32] += POW2[i];

    for (int i = 0; i < 128; i++)
        res[i] = temp[i];
}

void div_R(uint32_t res[128], uint32_t a[129])
{
    uint32_t temp[128] = {0};

    int word_shift = R_bits / 32;
    int bits_shift = R_bits % 32;
    for (int i = 0; i + word_shift < 129; i++)
        temp[i] = a[i + word_shift];

    if (bits_shift != 0)
    {
        for (int i = 0; i < 128 - 1; i++)
            temp[i] = temp[i] >> bits_shift | (temp[i + 1] << (32 - bits_shift));
        temp[127] = temp[127] >> bits_shift;
    }

    for (int i = 0; i < 128; i++)
        res[i] = temp[i];
}

void inv_fermat(uint32_t res[128], uint32_t a[128])
{
    uint32_t temp[128];
    sub(temp, P, TWO);
    mod_pow(res, a, temp);
}

void inv_exculid(uint32_t a[128], uint32_t b[128], uint32_t x[128], uint32_t y[128])
{
    if (equal(b, ZERO))
    {
        x[0] = 1;
        y[0] = 0;
        return;
    }

    uint32_t temp[128] = {0};
    div(temp, a, b);
    mul(temp, b, temp);
    sub(temp, a, temp);

    inv_exculid(b, temp, y, x);

    uint32_t temp2[128] = {0};
    div(temp2, a, b);
    mul(temp2, temp2, x);

    if (bigger(y, temp2))
    {
        sub(temp2, y, temp2);
        mod(temp2, temp2);
    }
    else
    {
        sub(temp2, temp2, y);
        mod(temp2, temp2);
        sub(temp2, P, temp2);
    }

    for (int i = 0; i < 128; i++)
        y[i] = temp2[i];
}

void inv_exculid2(uint32_t res[128], uint32_t a[128])
{
    uint32_t x1[128] = {0}, y1[128] = {1}, x2[128] = {0}, y2[128] = {0};
    for (int i = 0; i < 128; i++)
        x1[i] = a[i];
    for (int i = 0; i < 128; i++)
        x2[i] = P[i];

    while (!(equal(x1, ONE) || equal(x2, ONE)))
    {
        if ((x1[0] & 1) == 0 && (y1[0] & 1) == 0)
        {
            div(x1, x1, TWO);
            div(y1, y1, TWO);
        }
        else if ((x1[0] & 1) == 0 && (y1[0] & 1) == 1)
        {
            div(x1, x1, TWO);
            add(y1, y1, P);
            div(y1, y1, TWO);
        }

        if ((x2[0] & 1) == 0 && (y2[0] & 1) == 0)
        {
            div(x2, x2, TWO);
            div(y2, y2, TWO);
        }
        else if ((x2[0] & 1) == 0 && (y2[0] & 1) == 1)
        {
            div(x2, x2, TWO);
            add(y2, y2, P);
            div(y2, y2, TWO);
        }

        if (bigger(x1, x2))
        {
            mod_sub(x1, x1, x2);
            mod_sub(y1, y1, y2);
        }
        else
        {
            mod_sub(x2, x2, x1);
            mod_sub(y2, y2, y1);
        }
    }

    if (equal(x1, ONE))
        for (int i = 0; i < 128; i++)
            res[i] = y1[i];
    else
        for (int i = 0; i < 128; i++)
            res[i] = y2[i];
}

void exculid(uint32_t a[128], uint32_t b[128], uint32_t x[128], uint32_t y[128])
{
    if (equal(b, ZERO))
    {
        x[0] = 1;
        y[0] = 0;
        return;
    }

    uint32_t temp[128] = {0};
    div(temp, a, b);
    mul(temp, b, temp);
    sub(temp, a, temp);

    exculid(b, temp, y, x);

    uint32_t temp2[128] = {0};
    div(temp2, a, b);
    mul(temp2, temp2, x);

    if (bigger(y, temp2))
    {
        sub(temp2, y, temp2);
        mod_R(temp2, temp2);
    }
    else
    {
        sub(temp2, temp2, y);
        mod_R(temp2, temp2);
        sub(temp2, R, temp2);
    }

    for (int i = 0; i < 128; i++)
        y[i] = temp2[i];
}

void pre_cal()
{
    P_words = P_bits / 32;
    R_words = P_words + 1;
    R[R_words] = 1;
    R_bits = getBits(R);

    mod(R2, R);
    mul(R2, R2, R2);
    mod(R2, R2);

    uint32_t negP[128] = {0};
    sub(negP, R, P);

    uint32_t x[128] = {0},
             y[128] = {0};
    exculid(negP, R, x, y);
    for (int i = 0; i < 128; i++)
        P_[i] = x[i];
}

void REDC(uint32_t res[128], uint32_t a[128], uint32_t b[128])
{
    uint32_t T[128] = {0};

    mul(T, a, b);

    uint32_t m[128] = {0};
    uint32_t t[128] = {0};
    mod_R(m, T); // m=T mod R

    mul(m, m, P_); // m=(T mod R) * P_
    mod_R(m, m);   // m=((T mod R) * P_) mod R
    mul(t, m, P);  // t=m*P

    uint64_t carry = 0;
    uint32_t te[129] = {0};
    for (int i = 0; i < 128; i++)
    {
        uint64_t sum = carry + t[i] + T[i];
        te[i] = sum & 0xffffffff;
        carry = sum >> 32;
    }
    if (carry)
        te[128] = carry;

    div_R(t, te);

    if (bigger(t, P) || equal(t, P))
        sub(t, t, P);

    for (int i = 0; i < 128; i++)
        res[i] = t[i];
}

void mod_mul_mont(uint32_t res[128], uint32_t a[128], uint32_t b[128])
{
    uint32_t x_[128] = {0}, y_[128] = {0};

    REDC(x_, a, R2);
    REDC(y_, b, R2);
    REDC(x_, x_, y_);
    REDC(x_, x_, ONE);

    for (int i = 0; i < 128; i++)
        res[i] = x_[i];
}

void mod_pow_mont(uint32_t res[128], uint32_t a[128], uint32_t b[128])
{
    uint32_t temp[128] = {0};
    uint32_t base[128] = {0};
    REDC(temp, ONE, R2);
    REDC(base, a, R2);

    for (int i = 0; i < 128; i++)
    {
        if (bigger(temp, ONE) && b[i] == 0)
            break;
        for (int j = 0; j < 32; j++)
        {
            if (b[i] & POW2[j])
                REDC(temp, temp, base);

            REDC(base, base, base);
        }
    }

    REDC(res, temp, ONE);
}

int main(void)
{
    int n;
    cin >> n >> p;
    reverse(p.begin(), p.end());
    str2bi(P, p);
    P_bits = getBits(P);

    pre_cal();

    while (n--)
    {
        string a, b;
        cin >> a >> b;
        reverse(a.begin(), a.end());
        reverse(b.begin(), b.end());

        uint32_t A[128] = {0}, B[128] = {0};
        str2bi(A, a);
        str2bi(B, b);

        uint32_t res1[128] = {0};
        mod_add(res1, A, B);
        cout << bi2str(res1) << '\n';

        uint32_t res2[128] = {0};
        mod_sub(res2, A, B);
        cout << bi2str(res2) << '\n';

        uint32_t res3[128] = {0};
        mod_mul(res3, A, B);
        cout << bi2str(res3) << '\n';

        uint32_t x[128] = {0}, y[128] = {0};
        inv_exculid(A, P, x, y);
        cout << bi2str(x) << '\n';

        uint32_t res5[128] = {0};
        mod_pow_mont(res5, A, B);
        cout << bi2str(res5) << '\n';

        if (n)
            cout << '\n';
    }
}