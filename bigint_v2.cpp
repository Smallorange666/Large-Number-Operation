#include <string>
#include <stdint.h>
#include <algorithm>
#include <fstream>
#include <iostream>
using namespace std;

const int W = 64;
const int MAXLEN = 64;
const __uint128_t BASE = (static_cast<__uint128_t>(1) << 64);
const uint64_t MASK = 0xffffffffffffffff;

string p;
uint64_t P[MAXLEN];

uint64_t R[MAXLEN] = {0};
int P_bits = 0;
int P_words = 0;
int R_bits = 0;
int R_words = 0;
uint64_t R2[MAXLEN] = {0};
uint64_t P_[MAXLEN] = {0};
uint64_t ZERO[MAXLEN] = {0};
uint64_t ONE[MAXLEN] = {1};
uint64_t TWO[MAXLEN] = {2};
uint64_t POW2[64] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
                     1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288,
                     1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912,
                     1073741824, 2147483648, 4294967296, 8589934592, 17179869184, 34359738368, 68719476736, 137438953472, 274877906944, 549755813888,
                     1099511627776, 2199023255552, 4398046511104, 8796093022208, 17592186044416, 35184372088832, 70368744177664, 140737488355328, 281474976710656, 562949953421312,
                     1125899906842624, 2251799813685248, 4503599627370496, 9007199254740992, 18014398509481984, 36028797018963968, 72057594037927936, 144115188075855872, 288230376151711744, 576460752303423488,
                     1152921504606846976, 2305843009213693952, 4611686018427387904, 9223372036854775808ULL};

#ifndef DONLINE_JUDGE
ifstream in("input.txt");
ofstream out("myanswer.txt");
#define cin in
#define cout out
#endif

int getBits(const uint64_t a[MAXLEN])
{
    for (int i = MAXLEN - 1; i >= 0; i--)
    {
        if (a[i])
            for (int j = W - 1; j >= 0; j--)
                if (a[i] & POW2[j])
                    return i * W + j;
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

bool bigger(uint64_t a[MAXLEN], uint64_t b[MAXLEN])
{
    for (int i = MAXLEN - 1; i >= 0; i--)
    {
        if (a[i] > b[i])
            return true;
        else if (a[i] < b[i])
            return false;
    }
    return false;
}

bool equal(uint64_t a[MAXLEN], uint64_t b[MAXLEN])
{
    for (int i = 0; i < MAXLEN; i++)
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

void str2bi(uint64_t res[MAXLEN], string &s)
{
    int count = 0;
    while (s != "0")
    {
        if ((s[0] - '0') & 1)
            res[count / W] += POW2[count % W];
        str_div2(s);
        count++;
    }
}

string bi2str(uint64_t res[MAXLEN])
{
    string s = "0";
    for (int i = MAXLEN - 1; i >= 0; i--)
        for (int j = W - 1; j >= 0; j--)
        {
            str_mul2(s);
            if (res[i] & POW2[j])
                str_add1(s);
        }

    reverse(s.begin(), s.end());
    return s;
}

void add(uint64_t res[MAXLEN], uint64_t a[MAXLEN], uint64_t b[MAXLEN])
{
    uint64_t carry = 0;
    uint64_t temp[MAXLEN] = {0};
    for (int i = 0; i < MAXLEN; i++)
    {
        uint64_t sum;
        if (__builtin_add_overflow(a[i], carry, &sum) || __builtin_add_overflow(sum, b[i], &sum))
        {
            temp[i] = ((__uint128_t)a[i] + b[i] + carry) & MASK;
            carry = 1;
        }
        else
        {
            temp[i] = sum;
            carry = 0;
        }
    }

    for (int i = 0; i < MAXLEN; i++)
        res[i] = temp[i];
}

void sub(uint64_t res[MAXLEN], uint64_t a[MAXLEN], uint64_t b[MAXLEN])
{
    uint64_t borrow = 0;
    uint64_t temp[MAXLEN] = {0};
    for (int i = 0; i < MAXLEN; i++)
    {
        uint64_t gap;
        if (__builtin_sub_overflow(a[i], borrow, &gap) || __builtin_sub_overflow(gap, b[i], &gap))
        {
            temp[i] = (BASE + a[i] - borrow - b[i]) & MASK;
            borrow = 1;
        }
        else
        {
            temp[i] = gap;
            borrow = 0;
        }
    }

    for (int i = 0; i < MAXLEN; i++)
        res[i] = temp[i];
}

void mul(uint64_t res[MAXLEN], uint64_t a[MAXLEN], uint64_t b[MAXLEN])
{
    uint64_t temp[2 * MAXLEN] = {0};

    for (int i = 0; i < MAXLEN; i++)
    {
        uint64_t carry = 0;
        for (int j = 0; j < MAXLEN; j++)
        {
            __uint128_t sum = (__uint128_t)a[i] * b[j] + temp[i + j] + carry;
            temp[i + j] = sum & MASK;
            carry = sum >> W;
        }
        if (carry)
            temp[i + MAXLEN] += carry;
    }

    for (int i = 0; i < MAXLEN; i++)
        res[i] = temp[i];
}

void div(uint64_t res[MAXLEN], uint64_t a[MAXLEN], uint64_t b[MAXLEN])
{
    int n = getBits(b) / W;
    int m = getBits(a) / W - n; // iterate times

    uint64_t d = BASE / ((__uint128_t)b[n] + 1);
    uint64_t u_[MAXLEN + 1] = {0}, v_[MAXLEN] = {0};
    uint64_t carry = 0;
    // scale a
    for (int i = 0; i < MAXLEN; i++)
    {
        __uint128_t temp = (__uint128_t)a[i] * d + carry;
        u_[i] = temp & MASK;
        carry = temp >> W;
    }
    if (carry)
        u_[MAXLEN] = carry;

    // scale b
    for (int i = 0; i < MAXLEN; i++)
    {
        __uint128_t temp = (__uint128_t)b[i] * d + carry;
        v_[i] = temp & MASK;
        carry = temp >> W;
    }

    int j = m;
    while (j >= 0)
    {
        uint64_t tem[MAXLEN + 1] = {0};
        for (int i = 0; i <= n + 1; i++)
            tem[i] = u_[i + j];

        __uint128_t q_hat = (tem[n + 1] * BASE + tem[n]) / v_[n];
        if (q_hat > BASE - 1)
            q_hat = BASE - 1;

        uint64_t qv[MAXLEN] = {0};
        carry = 0;
        for (int i = 0; i < MAXLEN; i++)
        {
            __uint128_t temp = q_hat * v_[i] + carry;
            qv[i] = temp & MASK;
            carry = temp >> W;
        }
        while (bigger(qv, tem))
        {
            q_hat--;
            sub(qv, qv, v_);
        }

        sub(tem, tem, qv);
        for (int i = 0; i <= n + 1; i++)
            u_[i + j] = tem[i];

        res[j] = q_hat;

        j--;
    }
}

void mod(uint64_t res[MAXLEN], uint64_t a[MAXLEN])
{
    uint64_t temp[MAXLEN] = {0};
    div(temp, a, P);
    mul(temp, temp, P);
    sub(res, a, temp);
}

void mod_add(uint64_t res[MAXLEN], uint64_t a[MAXLEN], uint64_t b[MAXLEN])
{
    add(res, a, b);
    mod(res, res);
}

void mod_sub(uint64_t res[MAXLEN], uint64_t a[MAXLEN], uint64_t b[MAXLEN])
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

void mod_div(uint64_t res[MAXLEN], uint64_t a[MAXLEN], uint64_t b[MAXLEN])
{
    div(res, a, b);
    mod(res, res);
}

void mod_R(uint64_t res[MAXLEN], uint64_t a[MAXLEN])
{
    uint64_t temp[MAXLEN] = {0};
    for (int i = 0; i < R_words; i++)
        temp[i] = a[i];
    for (int i = 0; i < MAXLEN; i++)
        res[i] = temp[i];
}

void div_R(uint64_t res[MAXLEN], uint64_t a[MAXLEN + 1])
{
    uint64_t temp[MAXLEN] = {0};
    for (int i = 0; i + R_words < MAXLEN + 1; i++)
        temp[i] = a[i + R_words];
    for (int i = 0; i < MAXLEN; i++)
        res[i] = temp[i];
}

void inv_exculid(uint64_t a[MAXLEN], uint64_t b[MAXLEN], uint64_t x[MAXLEN], uint64_t y[MAXLEN])
{
    if (equal(b, ZERO))
    {
        x[0] = 1;
        y[0] = 0;
        return;
    }

    uint64_t temp[MAXLEN] = {0};
    div(temp, a, b);
    mul(temp, b, temp);
    sub(temp, a, temp);

    inv_exculid(b, temp, y, x);

    uint64_t temp2[MAXLEN] = {0};
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

    for (int i = 0; i < MAXLEN; i++)
        y[i] = temp2[i];
}

void exculid(uint64_t a[MAXLEN], uint64_t b[MAXLEN], uint64_t x[MAXLEN], uint64_t y[MAXLEN])
{
    if (equal(b, ZERO))
    {
        x[0] = 1;
        y[0] = 0;
        return;
    }

    uint64_t temp[MAXLEN] = {0};
    div(temp, a, b);
    mul(temp, b, temp);
    sub(temp, a, temp);

    exculid(b, temp, y, x);

    uint64_t temp2[MAXLEN] = {0};
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

    for (int i = 0; i < MAXLEN; i++)
        y[i] = temp2[i];
}

void pre_cal()
{
    P_words = P_bits / W;
    R_words = P_words + 1;
    R[R_words] = 1;
    R_bits = getBits(R);

    mod(R2, R);
    mul(R2, R2, R2);
    mod(R2, R2);

    uint64_t negP[MAXLEN] = {0};
    sub(negP, R, P);

    uint64_t x[MAXLEN] = {0},
             y[MAXLEN] = {0};
    exculid(negP, R, x, y);
    for (int i = 0; i < MAXLEN; i++)
        P_[i] = x[i];
}

void REDC(uint64_t res[MAXLEN], uint64_t a[MAXLEN], uint64_t b[MAXLEN])
{
    uint64_t T[MAXLEN] = {0};
    mul(T, a, b);

    int s = P_words + 1;
    uint64_t T_[MAXLEN + 1] = {0};
    for (int i = 0; i < MAXLEN; i++)
        T_[i] = T[i];

    for (int i = 0; i < s; i++)
    {
        uint64_t carry = 0;
        __uint128_t mi = ((__uint128_t)T_[i] * P_[0]) & MASK;
        for (int j = 0; j < s; j++)
        {
            __uint128_t temp = mi * P[j] + T_[i + j] + carry;
            T_[i + j] = temp & MASK;
            carry = temp >> W;
        }

        int count = 0;
        while (carry)
        {
            __uint128_t temp = (__uint128_t)T_[i + s + count] + carry;
            T_[i + s + count] = temp & MASK;
            carry = temp >> W;
            count++;
        }
    }

    div_R(res, T_);

    if (bigger(res, P) || equal(res, P))
        sub(res, res, P);
}

void mod_mul_trivial(uint64_t res[128], uint64_t a[128], uint64_t b[128])
{
    mul(res, a, b);
    mod(res, res);
}

void mod_mul_mont(uint64_t res[MAXLEN], uint64_t a[MAXLEN], uint64_t b[MAXLEN])
{
    uint64_t x_[MAXLEN] = {0}, y_[MAXLEN] = {0};

    REDC(x_, a, R2);
    REDC(y_, b, R2);
    REDC(x_, x_, y_);
    REDC(x_, x_, ONE);

    for (int i = 0; i < MAXLEN; i++)
        res[i] = x_[i];
}

void mod_pow_mont(uint64_t res[MAXLEN], uint64_t a[MAXLEN], uint64_t b[MAXLEN])
{
    uint64_t temp[MAXLEN] = {0};
    uint64_t base[MAXLEN] = {0};
    REDC(temp, ONE, R2);
    REDC(base, a, R2);

    for (int i = 0; i < MAXLEN; i++)
    {
        if (bigger(temp, ONE) && b[i] == 0)
            break;
        for (int j = 0; j < W; j++)
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

        uint64_t A[MAXLEN] = {0}, B[MAXLEN] = {0};
        str2bi(A, a);
        str2bi(B, b);

        uint64_t res1[MAXLEN] = {0};
        mod_add(res1, A, B);
        cout << bi2str(res1) << '\n';

        uint64_t res2[MAXLEN] = {0};
        mod_sub(res2, A, B);
        cout << bi2str(res2) << '\n';

        uint64_t res3[MAXLEN] = {0};
        mod_mul_trivial(res3, A, B);
        cout << bi2str(res3) << '\n';

        uint64_t x[MAXLEN] = {0}, y[MAXLEN] = {0};
        inv_exculid(A, P, x, y);
        cout << bi2str(x) << endl;

        uint64_t res5[MAXLEN] = {0};
        mod_pow_mont(res5, A, B);
        cout << bi2str(res5) << '\n';

        if (n)
            cout << '\n';
    }
}