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
//#include <stdint>
using namespace std;

#define BN_64 0

// Файл bn.h
template <int SIZE>
class BN {
#if BN_64
    using rt = unsigned __int128;       // Вычисления
    using irt = __int128;
    using mt = unsigned long long ; // Хранение
    static const int BITS = 64;
#else
    using rt = unsigned long long;       // Вычисления
    using irt = long long;
    using mt = unsigned; // Хранение
    static const int BITS = 32;
#endif
public:
    // The BN representation always keep an extra zero in front of to keep allocs rate low
    static const rt RADIX = (rt)1 << BITS;
    static const rt  MASK = RADIX - 1;
    //static const int SIZE = 96; // В словах
    static const int BSIZE = SIZE * sizeof(mt);

    BN(BN const &orig) { // Создать копию существующего BN
        memset(_b, 0, BSIZE);
        memcpy(_b, orig._b, orig._n * sizeof(mt));
        _n = orig._n;
        _sign = orig._sign;
    }

    explicit BN(double init_double) {
        char buf[40];
        sprintf(buf, "%.0f", init_double);
        init(buf, 10);
    }

    BN(long long init_int = 0) {   // Создать новое BN
        memset(_b, 0, BSIZE);
        _n = 2 + ((init_int > 0 ? init_int : -init_int) >= RADIX);
        assert(_n < SIZE);
        _b[_n - 1] = 0;
        _sign = 0;
        if (init_int < 0) {
            _sign = -1;
            init_int = -init_int;
        }
        _b[0] = (mt)(init_int % RADIX);
        _b[1] = (mt)(init_int / RADIX);
    }


    BN(string const &in, mt radix = 10) {
        if (in.size() > 0 && in[0] == '@') {
            init_hex(in.c_str() + 1);
        } else {
            init(in.c_str(), radix);
        }
    }

    BN(const char *init_string, mt radix = 10) {
        if (*init_string == '@') {
            init_hex(init_string + 1);
        } else {
            init(init_string, radix);
        }
    }


    ~BN() {
    }

    void print(const char *s) const {
        printf("%s=(s=%d n=%d b=[", s, _sign, _n);
        for (int i = 0; i < _n; i++) {
            printf("%08x ", _b[i]);
        }
        printf("])\n");
    }
    BN & mul_to_fast(mt n) {
        rt carry = 0;
        int i;
        for (i = 0; i < _n; i++) {
            rt p = _b[i] * (rt)n + carry;
            carry = (rt)(p / RADIX);
            _b[i] = p % RADIX;
        }
        canonify(_b, _n);
        return *this;
    }

    BN & div_to_fast(mt val) {
        mt rem;
        divmod_to_fast(val, rem);
        return *this;
    }

    BN & add_to_fast(mt n) {
        rt carry = n;
        int i;
        for (i = 0; i < _n && carry != 0; i++) {
            rt p = _b[i] + carry;
            carry = (rt)(p / RADIX);
            _b[i] = p % RADIX;
        }
        canonify(_b, _n);
        return *this;
    }

    bool is_zero() const {
        return _n <= 2 && _b[0] == 0;
    }

    BN & neg() {
        _sign = !_sign;
        return *this;
    }

    BN & abs() { // Взять модуль
        _sign = 0;
        return *this;
    }

    // Вернуть <0, если отрицательное, >0, если положительное, 0.
    int sign() const {
        if (_sign) return -1;
        if (is_zero()) return 0;
        return 1;
    }

    // Выдать представление BN в системе счисления radix в виде строки
    string to_string(int radix = 10) const {
        static const char r[] = "0123456789abcdefghijklmnopqrstuvwxyz";
        BN tmp(*this);
        int sign = tmp._sign;
        char *s = (char *)alloca(_n * 34 * sizeof(char));
        char *ps = s;
        tmp._sign = 0;
        if (s == NULL) return NULL;
        while (tmp._n > 2) {
            mt rem = 0;
            tmp.divmod_to_fast(radix, rem);
            *ps++ = r[rem];
        }
        mt q = tmp._b[0];
        if (q == 0) {
            *ps++ = '0';
        } else {
            while (q > 0) {
                *ps++ = r[q % radix];
                q /= radix;
            }
        }
        if (sign) *ps++ = '-';
        size_t len = ps - s;
        for (size_t i = 0; i < len / 2; i++) {
            char c = s[i];
            s[i] = s[len - 1 - i];
            s[len - 1 - i] = c;
        }
        *ps++ = 0;
        string ret(s);
        return ret;

    }

    // Операции, аналогичные +=, -=, *=, /=, %=
    BN & add_to(BN const &r) {
        if (_sign == r._sign) {
            if (_n >= r._n) {
                bn_add_internal_onplace(_b, _n, r._b, r._n);
                canonify(_b, _n);
                return *this;
            }
            mt rb[SIZE]={0}; 
            memcpy(rb, r._b, r._n * sizeof(mt));
            bn_add_internal_onplace(rb, r._n, _b, _n);
            _n = r._n;
            canonify(rb, _n);
            memcpy(_b, rb, _n * sizeof(mt));
            return *this;
        }
        // 7 + 10 = 17
        // -7 + 10 = 3
        // 7 + -10 = -3
        // -7 + -10 = -17
        // here sign t != sign r
        int cmp = bn_mt_cmp(_b, _n, r._b, r._n);
        if (cmp == 0) {
            zeroify();
        } else if (cmp < 0) {
            mt nb[SIZE]={0};
            memcpy(nb, r._b, r._n * sizeof(mt));
            bn_sub_internal_onplace(nb, r._n, _b, _n);
            memcpy(_b, nb, sizeof _b);
            _n = r._n;
            _sign = !_sign; // -7 + 10 -> s=0, 7 + -10 ->s=1;
            canonify(_b, _n);
        } else {   // cmp > 0; 10 + -7 ->s=1; 
            bn_sub_internal_onplace(_b, _n, r._b, r._n);
            canonify(_b, _n);
        }
        return *this;
    }

    BN & twice() {
        bn_add_internal_onplace(_b, _n, _b, _n);
        canonify(_b, _n);
        return *this;
    }

    static BN gcd_mod(BN const &l, BN const &r) {
        //printf("gcd(%s,%s): ", l.to_string(10).c_str(), r.to_string(10).c_str());
        do {
            BN t = l.mod(r);
            l = r;
            r = t;
        } while (!r.is_zero());
        //printf(" = %s\n", l.to_string(10).c_str());
        return l;
    }

    size_t least_zeroes() const {
        size_t ret = 0;
        assert(_n >= 2);
        for (int i = _n-2; i >= 0; i--) {
            int word_zeroes = least_zero_count(_b[i]);
            ret += word_zeroes;
            if (word_zeroes != BITS) {
                break;
            }
        }
        return ret;
    }
    
    BN & sub_to(BN const &r) {
        BN rc;
        memcpy(rc._b, r._b, BSIZE);
        rc._n = r._n;
        rc._sign = !r._sign;
        add_to(rc);
        return *this;
    }

    BN & mul_to(BN const &r) {
        _sign ^= r._sign;
        auto tmp_n = r._n + _n;
        assert(r._n + _n <= SIZE - 1);
        rt tmp_b[SIZE]={0}; 
        for (int i = 0; i < r._n - 1; i++) {
            rt *pb = tmp_b + i;
            for (int j = 0; j < _n - 1; j++, pb++) {
                pb[0] += (rt)r._b[i] * (rt)_b[j];
                pb[1] += (rt)tmp_b[i + j] / RADIX;
                pb[0] %= RADIX;
            }
        }
        for (int i = 0; i < tmp_n - 2; i++) {
            tmp_b[i + 1] += tmp_b[i] / RADIX;
            tmp_b[i] %= RADIX;
        }
        _n = tmp_n;
        for (int i = 0; i < _n; i++) {
            _b[i] = (mt)tmp_b[i];
        }
        canonify(_b,_n);
        return *this;
    }
    
    BN & div_to(BN const &r) {
        if (r.is_zero()) throw "Divide by zero";
        int sign = _sign ^ r._sign;
        BN rem;
        divmod_internal(this, &r, &rem);
        _sign = sign;
        if (_sign && !rem.is_zero()) {
            BN one(1);
            rem.add_to(r);
            sub_to(one);
        }
        return *this;
    }

    BN & mod_to(BN const &r) {
        if (r.is_zero()) throw "Divide by zero";
        int sign = _sign ^ r._sign;
        BN rem;
        divmod_internal(this, &r, &rem);
        _sign = sign;
        if (_sign && !rem.is_zero()) {
            BN one(1);
            rem.add_to(r);
            sub_to(one);
        }
        if (rem._n == 1) {
            assert(rem._b[0] == 0);
            rem._n = 2;
            rem._b[1] = 0;
        }
        memcpy(_b, rem._b, BSIZE);
        _n = rem._n;
        _sign = rem._sign;
        return *this;
    }


    BN &set_two_pow(int n) {  // Присвоить числу степень двойки
        _n = n / BITS + 2;
        ::memset(_b, 0, BSIZE);
        assert(_n < BSIZE);
        int bit = n % BITS;
        _b[_n - 2] |= ((rt)1 << bit);
        _sign = 0;
        return *this;
    }

    // Возвести число в степень degree
    BN & pow_to(int degree) {
        if (degree <= 0) throw "Power out of range";
        BN q(1);
        while (degree > 0) {
            if (degree & 1) {
                q.mul_to(*this);
            }
            if (degree > 1)
                mul_to(*this);
            degree /= 2;
        }
        myclone(q);
        return *this;
    }

    // Извлечь корень степени reciprocal из BN (бонусная функция)
    BN & root_to(int root) {
        if (root < 1) throw "Root < 1";
        BN one(1);
        if (cmp(one) == 0) return *this;
        if (root == 1) return *this;
        if (sign() < 0) throw "root from negative";
        int bt = (bits() - 1) / root + 1;
        BN l; l.set_two_pow(bt - 1);
        BN r; r.set_two_pow(bt);
        BN retval;
        BN c(l);
        BN f = l.pow(root);
        if (BN::cmp(f) == 0) {
            retval = l;
            goto cleanup;
        }
        f = r;
        f.pow_to(root);
        if (cmp(f) == 0) {
            retval = r;
            goto cleanup;
        }
        for (;;) {
            c = r.sub(l);
            int ccmp = cmp(c, one);
            if (ccmp <= 0) {
                retval = l;
                goto cleanup;
            }
            c = l.add(r).div_to_fast(2);
            f = c.pow(root);
            ccmp = cmp(f);
            if (ccmp == 0) {
                retval = c;
                goto cleanup;
            }
            if (ccmp > 0) {
                l = c;
            } else {
                r = c;
            }
        }
    cleanup:
        myclone(retval);
        return *this;
    }

    BN &operator=(BN const &orig) { // Создать копию существующего BN
        assert(orig._n < SIZE);
        memcpy(_b, orig._b, orig._n * sizeof(mt));
        _n = orig._n;
        _sign = orig._sign;
        return *this;
    }

    // Аналоги операций x = l+r (l-r, l*r, l/r, l%r) 
    BN  add(BN const &right) const {
        BN res(*this);
        return res.add_to(right);
    }

    BN  sub(BN const &right) const {
        BN res(*this);
        return res.sub_to(right);
    }

    BN  mul(BN const &right) const {
        BN res(*this);
        return res.mul_to(right);
    }

    BN  div(BN const &right) const {
        BN res(*this);
        return res.div_to(right);
    }

    BN  mod(BN const &right) const {
        BN res(*this);
        return res.mod_to(right);
    }


    BN pow(int p) const {
        BN res(*this);
        return res.pow_to(p);
    }

    // Если первое меньше, вернуть <0, 0, если равны, >0, если больше
    static int cmp(BN const &left, BN const &right) {
        if (!left._sign) {
            return !right._sign ? bn_mt_cmp(left._b, left._n, right._b, right._n) : 1;
        } else {
            return  right._sign ? -bn_mt_cmp(left._b, left._n, right._b, right._n) : -1;
        }
        return 0;
    }

    int cmp(BN const &right) const {
        if (!_sign) {
            return !right._sign ? bn_mt_cmp(_b, _n, right._b, right._n) : 1;
        } else {
            return  right._sign ? -bn_mt_cmp(_b, _n, right._b, right._n) : -1;
        }
        return 0;
    }

    bool operator<(BN const &right) const {
        return cmp(right) < 0;
    }

    bool operator>(BN const &right) const {
        return cmp(right) > 0;
    }

    bool operator<=(BN const &right) const {
        return cmp(right) <= 0;
    }

    bool operator>=(BN const &right) const {
        return cmp(right) >= 0;
    }


    bool operator==(BN const &right) const {
        return cmp(right) == 0;
    }

    bool operator!=(BN const &right) const {
        return cmp(right) != 0;
    }


    BN & operator +=(const BN &oth) {
        add_to(oth);
        return *this;
    }

    BN & operator +=(mt oth) {
        add_to_fast(oth);
        return *this;
    }

    BN &operator *=(const BN &oth) {
        mul_to(oth);
        return *this;
    }

    BN& operator -=(const BN &oth) {
        sub_to(oth);
        return *this;
    }

    BN& operator /=(const BN &oth) {
        div_to(oth);
        return *this;
    }

    BN& operator %=(const BN &oth) {
        mod_to(oth);
        return *this;
    }

    BN operator +(const BN&r) const {
        BN ret(*this);
        return ret += r;
    }

    BN operator +(mt r) const {
        BN ret(*this);
        return ret += r;
    }

    BN operator -(const BN&r) const {
        BN ret(*this);
        return ret -= r;
    }

    BN operator *(const BN&r) const {
        BN ret(*this);
        return ret *= r;
    }

    BN operator /(const BN&r) const {
        BN ret(*this);
        return ret /= r;
    }

    BN operator %(const BN&r) const {
        BN ret(*this);
        return ret %= r;
    }

    int operator %(mt mod) const {
        BN ret(*this);
        ret %= mod;
        return ret._sign ? -(int)ret._b[0] : (int)ret._b[0];
    }

    double get_double() const { // Veeeeeeeryyyy ugly
        string im = to_string(10);
        return atof(im.c_str());
    }

    int chunks() const {
        return _n;
    }

    const void *addr() const {
        return (void *)_b;
    }
    
    size_t size() const {
        return _n * sizeof(mt);
    }
    
    int bits() const {
        return (int)((_n - 2)*BITS) + most_zero_count(_b[_n - 2]);
    }



private:
    int _n;
    int _sign;
    mt _b[SIZE]={0};

   
    void zeroify() {
        if (_n == 2) {
            _b[0] = 0;
        } else {
            _b[0] = _b[1] = 0;
            _n = 2;
        }
        _sign = 0;
    }

    static int most_zero_count(mt x) {
        int n = 0;
        while (x != 0) {
            n++;
            x >>= 1;
        }
        return n;
    }

    static int least_zero_count(mt x) {
        printf("lzc(%u)=",x);
        if (x == 0) {
            printf("32\n");
            return SIZE;
        }
        int n = 0;
        while ((x & 1)==0) {
            n++;
            x >>= 1;
        }
        printf("%d\n",n);
        return n;
    }

    void fast_clone(BN const &orig) { // Создать копию существующего BN
        memcpy(_b, orig._b, orig._n * sizeof(mt));
        _n = orig._n;
        _sign = orig._sign;
    }

    static void canonify(mt *b, int &n) {
        if (b[n - 1] != 0) {
            assert(n + 1 < SIZE);
            b[n] = 0;
            n++;
        } else {
            // leading zero present
            int nn = n;
            while (nn > 0 && b[nn - 1] == 0) {
                nn--;
            }
            if (nn < n - 1 && n > 2) {
                n = nn + 1;
            }
        }
    }

    int divmod_to_fast(mt radix, mt &rem) {
        rt carry = 0;
        for (int i = _n; --i >= 0; ) {
            rt d = _b[i] + carry*RADIX;
            rt r = d / radix;
            rt q = d % radix;
            _b[i] = (mt)r;
            carry = (mt)q;
        }
        rem = (mt)carry;
        canonify(_b, _n);
        return 0;
    }

#if BN_64
    static int nlz(unsigned long long x) {
        int n = 1;
        if (x == 0) {
            return 64;
        }
        if ((x >> 32) == 0) {
            n += 32;
            x <<= 32;
        }
        if ((x >> 48) == 0) {
            n += 16;
            x <<= 16;
        }
        if ((x >> 56) == 0) {
            n += 8;
            x <<= 8;
        }
        if ((x >> 60) == 0) {
            n += 4;
            x <<= 4;
        }
        if ((x >> 62) == 0) {
            n += 2;
            x <<= 2;
        }
        n -= (x >> 63);
        return n;
    }
#else
    static int nlz(unsigned x) {
        int n = 1;
        if (x == 0) {
            return 32;
        }
        if ((x >> 16) == 0) {
            n += 16;
            x <<= 16;
        }
        if ((x >> 24) == 0) {
            n += 8;
            x <<= 8;
        }
        if ((x >> 28) == 0) {
            n += 4;
            x <<= 4;
        }
        if ((x >> 30) == 0) {
            n += 2;
            x <<= 2;
        }
        n -= (x >> 31);
        return n;
    }
#endif

    static bool divmnu(mt *q, mt *r, const mt *u, const mt *v, int m, int n) {
        const rt b = RADIX;
        rt qhat;
        rt rhat;
        // rt p;
        irt s, t, k;
        //printf("argument u=%d; [", m);  for (int i = 0; i < m; i++) printf("%08x ", u[i]); printf("]\n");
        //printf("argument v=%d; [", n);  for (int i = 0; i < n; i++) printf("%08x ", v[i]); printf("]\n");
        if (m < n || n <= 0 || v[n - 1] == 0) {
            return false;
        }
        if (n == 1) {
            k = 0;
            for (int j = m - 1; j >= 0; j--) {
                rt kb = k*b;
                rt kbuj = kb + u[j];
                rt qj = kbuj / v[0];
                q[j] = (mt)qj;
                k = kbuj - qj*v[0];
                //q[j] = (mt)((k*b + u[j]) / v[0]);
                //k = (mt)((k*b + u[j]) - q[j] * v[0]);
            }
            r[0] = (mt)k;
            return true;
        }
        s = nlz(v[n - 1]);
        mt vn[SIZE]; //= { 0 };
        mt un[SIZE]; //= { 0 };
        if (s != 0) {
            //printf("s=%lld\n", s);
            for (int i = n - 1; i > 0; i--) {
                vn[i] = (v[i] << s) | (v[i - 1] >> (BITS - s));
            }
            vn[0] = v[0] << s;
            un[m] = u[m - 1] >> (BITS - s);
            for (int i = m - 1; i > 0; i--) {
                un[i] = (u[i] << s) | (u[i - 1] >> (BITS - s));
            }
            un[0] = u[0] << s;
            //printf("un=%d; [", n);  for (int i = 0; i <= m; i++) printf("%08x ", un[i]); printf("]\n");
        } else {
            //printf("s=0\n");
            for (int i = n - 1; i > 0; i--) {
                vn[i] = v[i];
            }
            vn[0] = v[0];
            un[m] = 0;
            for (int i = m - 1; i > 0; i--) {
                un[i] = u[i];
            }
            un[0] = u[0];
        }
        //printf("original un=%d; [", m);  for (int i = 0; i <= m; i++) printf("%08x ", un[i]); printf("]\n");
        //printf("original vn=%d; [", n);  for (int i = 0; i < n; i++) printf("%08x ", vn[i]); printf("]\n");

        for (int j = m - n; j >= 0; j--) {
            rt tmp = un[j + n] * b + un[j + n - 1];
            qhat = (tmp /        vn[n - 1]);
            rhat = (tmp - qhat * vn[n - 1]);
            
        again:
            //printf("qhat=%08llx rhat=%08llx\n", qhat, rhat);
            if (qhat >= b || qhat*vn[n - 2] > b*rhat + un[j + n - 2]) {
                qhat--;
                rhat += vn[n - 1];
                if (rhat < b) goto again;
            }
            k = 0;
           //printf("before un=%d; [", m);  for (int i = 0; i <= m; i++) printf("%08x ", un[i]); printf("]\n");
            for (int i = 0; i < n; i++) {
                rt p = qhat*vn[i];
                rt plow = p & MASK;
                rt phigh = p >> BITS;
                t = un[i + j] - k - plow;
                rt ut = (rt)t;
                un[i + j] = (mt)(ut & MASK);
                k = phigh - (t >> BITS);
                //printf("t=%llx k=%lld\n", t, k);
            }
            t = un[j + n] - k;
            rt ut = (rt)t;
            un[j + n] = (mt)(ut & MASK);
            //printf("result t=%lld un=%d; [", t, m);  for (int i = 0; i <= m; i++) printf("%08x ", un[i]); printf("]\n");

            q[j] = (mt)qhat;
            if (t < 0) {
                q[j]--;
                k = 0;
                for (int i = 0; i < n; i++) {
                    t = un[i + j];
                    t += vn[i];
                    t += k;
                    //assert(t >= 0);
                    rt ut = t;
                    un[i + j] = (mt)(ut & MASK);
                    assert(t >= 0);
                    k = t >> BITS;
                }
                un[j + n] = (mt)(un[j+n]+ k);

            }
            //printf("result after back_off un=%d; [", m);  for (int i = 0; i <= m; i++) printf("%08x ", un[i]); printf("]\n");
        }
        //printf("result before back shift s=%lld un=%d; [", s, m);  for (int i = 0; i <= m; i++) printf("%08x ", un[i]); printf("]\n");
        if (s == 0) {
            for (int i = 0; i < n; i++) {
                r[i] = un[i];
            }
        } else {
            for (int i = 0; i < n; i++) {
                r[i] = (un[i] >> s) | (un[i + 1] << (BITS - s));
            }
        }
        r[n] = 0;
        //printf("totals result r=%d; [", n);  for (int i = 0; i <= n; i++) printf("%08x ", r[i]); printf("]\n");
        return true;
    }

    static int divmod_internal(BN *l, const BN *r, BN *rem) {
        if (bn_mt_cmp(l->_b, l->_n, r->_b, r->_n) < 0) {
            rem->_sign = 0;
            rem->_n = l->_n;
            memcpy(rem->_b, l->_b, l->_n * sizeof(mt));
            l->zeroify();
            return 0;
        }
        mt Q[SIZE] = { 0 };
        mt R[SIZE] = { 0 };
        int ret = divmnu(Q, R, l->_b, r->_b, l->_n - 1, r->_n - 1);
        memcpy(l->_b, Q, l->_n * sizeof(mt));
        //for (int i = 0; i < l->_n; i++) {
        //    l->_b[i] = Q[i];
        //}
        canonify(l->_b, l->_n);
        rem->_n = r->_n;
        //for (int i = 0; i < rem->_n; i++) {
        //    rem->_b[i] = R[i];
        //}
        memcpy(rem->_b, R, rem->_n * sizeof(mt));
        rem->_b[rem->_n-1]=0;
        canonify(rem->_b, rem->_n);

        return 0;
    }

    static int bn_mt_cmp(mt const *l, int ln, mt const *r, int rn) {
        if (ln > rn || rn < 0) {
            return 1;
        }
        if (ln < rn || ln < 0) {
            return -1;
        }
        int i;
        for (i = ln; --i >= 0; ) {
            if (l[i] < r[i]) {
                return -1;
            } else if (l[i] > r[i]) {
                return 1;
            }
        }
        return 0;
    }

    void myclone(BN const &orig) {
        assert(orig._n < SIZE);
        memcpy(_b, orig._b, orig._n * sizeof(mt));
        _n = orig._n;
        _sign = orig._sign;
    }

    static void bn_add_internal_onplace(mt *l, int ls, mt const*r, int rs) {
        assert(ls >= rs);
        int i;
        rt carry = 0;
        for (i = 0; i < rs; i++) {
            rt d = (rt)l[i] + (rt)r[i] + carry;
            carry = d / RADIX;
            l[i] = (mt)(d %= RADIX);
        }
        for (i = rs; i < ls && carry > 0; i++) {
            rt d = (mt)l[i] + carry;
            carry = d / RADIX;
            l[i] = (mt)(d %= RADIX);
        }
    }

    static void bn_sub_internal_onplace(mt *l, int ls, mt const *r, int rs) {
        assert(ls >= rs);
        int i;
        irt borrow = 0;
        for (i = 0; i < rs; i++) {
            irt d = (irt)l[i] - (irt)r[i] - borrow;
            borrow = 0;
            if ((irt)d < 0) {
                borrow = 1;
                d += RADIX;
            }
            l[i] = (mt)d;
        }
        for (i = rs; i < ls && borrow > 0; i++) {
            irt d = (irt)l[i] - borrow;
            borrow = 0;
            if ((irt)d < 0) {
                borrow = 1;
                d += RADIX;
            }
            l[i] = (mt)d;
        }
    }
    
    // Сдвинуть влево на bits < SIZE бит. Должно быть выделено ровно на 1 чанк больше.
    int lshift_internal(mt *v, int n, size_t bits) {
        assert(n < SIZE);
        assert(n >= 2);
        mt t = v[n-2] >> (BITS - bits);
        if (t == 0) {  // Не надо расширять
            for (int i = n - 2; i >= 0; i--) {
                v[i] = (v[i] << bits) | (v[i - 1] >> (BITS - bits));
            }
        } else {
            mt vn[SIZE] = { 0 };
            for (int i = n - 2; i >= 0; i--) {
                vn[i] = (v[i] << bits) | (v[i - 1] >> (BITS - bits));
            }
            vn[n-1] = t;
            for (int i = 0; i < n; i++) {
                v[i] = vn[i];
            }
            n++;
        }
        return n;
        //         auto leading_zeroes = nlz(v[n - 2]);
//         auto new_n = n;
//         if (leading_zeroes < (int)bits) {
//             new_n++;
//         }
//         mt vn[SIZE] = { 0 };
//         vn[0] = v[0] << bits;
//         for (int i = 0; i < new_n; i++) {
//             v[i] = vn[i];
//         }
//         return new_n;
    }

    // Сдвинуть вправо на bits < SIZE бит. 
    int rshift_internal(mt *v, int n, size_t bits) {
        assert(n < SIZE);
        //auto leading_zeroes = nlz(v[n - 2]);
        mt vn[SIZE] = { 0 };
        for (int i = 0; i < n-1; i++) {
            vn[i] = (v[i] >> bits) | (v[i + 1] << (BITS - bits));
        }
        for (int i = 0; i < n; i++) {
            v[i] = vn[i];
        }
        if (v[n - 2] == 0) {
            return n - 1;
        }
        return n;
    }

    void init(const char *init_string, mt radix = 10) {
        memset(_b, 0, BSIZE);
        while (*init_string == ' ' || *init_string == '\t')
            init_string++;
        _sign = 0;
        _n = 2;
        _b[0] = 0; _b[1] = 0;
        if (*init_string == '-') {
            _sign = !_sign;
            init_string++;
        }
        if (radix <= 10) {
            while (*init_string >= '0' && *init_string < (int)('0' + radix)) {
                mul_to_fast(radix);
                add_to_fast(*init_string - '0');
                init_string++;
            }
        } else {
            for (;;) {
                if (*init_string >= '0' && *init_string <= '9') {
                    mul_to_fast(radix);
                    add_to_fast(*init_string - '0');
                    init_string++;
                } else if (*init_string >= 'a' && *init_string < (int)('a' - 10 + radix)) {
                    mul_to_fast(radix);
                    add_to_fast(*init_string - 'a' + 10);
                    init_string++;
                } else if (*init_string >= 'A' && *init_string < (int)('A' - 10 + radix)) {
                    mul_to_fast(radix);
                    add_to_fast(*init_string - 'A' + 10);
                    init_string++;
                } else {
                    break;
                }
            }

        }
    }
 

    void init_hex(const char *init_string) {
        memset(_b, 0, BSIZE);
        _sign = 0;
        _n = 0;
        if (*init_string == '-') {
            _sign = !_sign;
            init_string++;
        }
        int len = (int)strlen(init_string);
        auto fromhex = [](int c) -> unsigned {if (c >= 'a') return c - 'a' + 10; else return c - '0'; };
        for (int p = len - 8 ; p >= 0; p -= 8, _n++) {
            unsigned z = 0;
            z |= fromhex(init_string[p]) << 28;
            z |= fromhex(init_string[p+1]) << 24;
            z |= fromhex(init_string[p+2]) << 20;
            z |= fromhex(init_string[p+3]) << 16;
            z |= fromhex(init_string[p+4]) << 12;
            z |= fromhex(init_string[p+5]) << 8;
            z |= fromhex(init_string[p+6]) << 4;
            z |= fromhex(init_string[p+7]);
            _b[_n] = (mt)z;
        }
        if (len % 8 == 0) {
            _b[_n] = 0;
            _n++;
            return;
        }
        unsigned z = 0;
        for (int p = 0; p < len % 8; p++) {
            z <<= 4;
            z |= fromhex(init_string[p]);
        }
        _b[_n] = z;
        _n++;
        _b[_n] = 0;
        _n++;
        assert(_n < SIZE);
    }

    void init_hex4(const char *init_string) {
        memset(_b, 0, BSIZE);
        _sign = 0;
        _n = 0;
        if (*init_string == '-') {
            _sign = !_sign;
            init_string++;
        }
        int len = (int)strlen(init_string);
        auto fromhex = [](int c) -> unsigned {if (c >= 'a') return c - 'a' + 10; else return c - '0'; };
        for (int p = len - 4; p >= 0; p -= 4, _n++) {
            unsigned z = 0;
            z |= fromhex(init_string[p]) << 12;
            z |= fromhex(init_string[p + 1]) << 8;
            z |= fromhex(init_string[p + 2]) << 4;
            z |= fromhex(init_string[p + 3]);
            _b[_n] = (mt)z;
        }
        if (len % 4 == 0) {
            _b[_n] = 0;
            _n++;
            return;
        }
        unsigned z = 0;
        for (int p = 0; p < len % 4; p++) {
            z <<= 4;
            z |= fromhex(init_string[p]);
        }
        _b[_n] = z;
        _n++;
        _b[_n] = 0;
        _n++;
        assert(_n < SIZE);
    }
public:
};



