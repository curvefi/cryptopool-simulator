#include <cassert>
#include <cstdio>
#include <map>
#include <unordered_map>
#include <string>
#include <ctime>
#include <optional>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <utility>
#include <sys/fcntl.h>
#include <sys/mman.h>
#include <stdexcept>
#include <cmath>
#include <memory.h>
#include <fstream>
#include <iostream>
#include "json.hpp"
using nlohmann::json;
//#include "bn_fixed.h"
#define DEBUG 0
static int trace = DEBUG;
#if DEBUG > 1
#define T printf("Entering %s\n", __PRETTY_FUNCTION__)
#define DBG printf("Now in %s line %d\n", __PRETTY_FUNCTION__ , __LINE__)
#define E printf("Leaving %s\n", __PRETTY_FUNCTION__)
#else
#define T
#define DBG
#define E
#endif
#if DEBUG > 0
#define P(x)      if (trace) printf("[%d] %s::%s=%s ", __LINE__, __FUNCTION__, #x, x.to_string(10).c_str())
#define DP(x)      if (trace) printf("[%d] %s::%s=%.12Lf ", __LINE__, __FUNCTION__, #x, (long double)x)
#define P256(x)    if (trace) { printf("[%d] %s::%s={", __LINE__, __FUNCTION__, #x);    for (size_t i = 0; i < x.size(); i++) printf("%s ", x[i].to_string(10).c_str()); printf("}\n"); }
#define EOL     if (trace) printf("\n")
#else
#define P(x)
#define DP(x)
#define P256(x)
#define EOL
#endif

using std::vector, std::string, std::pair, std::sort, std::map, std::min, std::max;

using u64 = unsigned long long;
using money = long double;
static const int MAX_ARRAY = 3;


static void print_clock(string const &mesg, clock_t start, clock_t end) {
    printf("%s %.3lf sec\n", mesg.c_str(), double(end - start) / CLOCKS_PER_SEC);
}

struct trade_data {
    u64 t = 0;          // 0
    money open = 0;    // 1
    money high = 0;    // 2
    money low = 0;     // 3
    money close = 0;   // 4
    money volume = 0;  // 5
    pair<int,int> pair1 = {0,0};
    void print() const {
        printf("{ open: %.6Lf, high: %.6Lf low: %.6Lf close: %.6Lf t: %llu, volume: %.6Lf pair: (%d,%d)} ",
               this->open, this->high, this->low, this->close, this->t, this->volume, this->pair1.first, this->pair1.second);
    }
};

struct trade_one {
    u64 t;
    trade_data trade;
};

money mabs(money val) noexcept {
    return val >= 0 ? val : -val;
}

void debug_print(string const &head, vector<trade_data> const &t, int count) {
    printf("%s:\n", head.c_str());
    size_t first = 0, last = count;
    if (count < 0) {
        first = t.size() + count;
        last = t.size();
    }
    for (auto i = first; i < last; i++) {
        t[i].print(); printf("\n");
    }
    printf("\n");
}

void debug_print(string const &head, vector<trade_one> const &t, int count) {
    printf("%s:\n", head.c_str());
    size_t first = 0, last = count;
    if (count < 0) {
        first = t.size() + count;
        last = t.size();
    }
    for (auto i = first; i < last; i++) {
        t[i].trade.print(); printf("\n");
    }
    printf("\n");
}

struct mapped_file {
    int fd;
    unsigned char *base = nullptr;
    size_t size = 0;
    string name;
    explicit mapped_file() {}
    bool map(string const &name) {
        fd = open(name.c_str(), O_RDONLY);
        if (fd < 0) return false;
        lseek(fd, 0, SEEK_END);
        size = lseek(fd, 0, SEEK_CUR);
        base = (unsigned char *)::mmap(nullptr, size, PROT_READ, MAP_NOCACHE|MAP_FILE|MAP_SHARED, fd, 0);
        printf("mapped_file::open: base=%p\n", base);
        if (base == MAP_FAILED) {
            perror(name.c_str());
            base = nullptr;
            return false;
        }
        this->name = name;
        return true;
    }
    ~mapped_file() {
        if (base != nullptr) {
            ::munmap(base, size);
        }
        if (fd >= 0) {
            close(fd);
        }
    }
};


// py: returns list of dicts ['t'->u64, 'open'->float, 'high'->float, 'low'->float, 'close'->float, 'volume'->float]
// c++ returns vector of struct datum
vector<trade_data> get_data(std::string const &fname) {
    auto start_time = clock();
    auto name_to_open = "download/" + fname + "-1m.json";
    printf("parsing %s\n", name_to_open.c_str());
    mapped_file mf;
    if (!mf.map(name_to_open)) {
        printf("failed to open '%s'\n", name_to_open.c_str());
        abort();
    }
    vector<trade_data> ret;
    auto p = mf.base;
    if (*p == '[') p++; // skip initial '[';
    auto scan_double = [] (unsigned char *p, long double *d) {
        if (*p == '"') p++;
        *d = atof((char *)p);
        while (*p != '"') p++;
        p++; // skip '"';
        while (*p == ',' || *p == ' ') p++;
        return p;
    };
    auto scan_u64 = [] (unsigned char *p, u64 *d) {
        u64 ret = 0;
        while (*p >= '0' && *p <= '9') {
            ret = ret * 10 + *p - '0';
            p++;
        }
        while (*p == ',' || *p == ' ') p++;
        *d = ret;
        return p;
    };
    while (*p != ']') {
        // [1503443580000, "3984.00000000", "3984.00000000", "3984.00000000", "3984.00000000", "0.46619400", 1503443639999, "1857.31689600", 2, "0.00000000", "0.00000000", "11761.90492277"],
        if (*p == '[') {
            trade_data d;
            p++;
            p = scan_u64(p, &d.t); d.t /= 1000;
            p = scan_double(p, &d.open);
            p = scan_double(p, &d.high);
            p = scan_double(p, &d.low);
            p = scan_double(p, &d.close);
            p = scan_double(p, &d.volume);
            ret.push_back(d);
            while (*p != ']') p++;
            p++; // skip ']'
        } else p++;
    }
    auto end_time = clock();
    printf("%s: load %zu elements\n", name_to_open.c_str(), ret.size());
    print_clock("parsing took", start_time, end_time);
    return ret;

}

auto get_price_vector(int n, vector<trade_data> const &data) {
    vector<money> p(n);
    p[0] = 1.L;
    for (auto const &d: data) {
        if (d.pair1.first == 0) {
            p[d.pair1.second] =  d.close;
        }
        bool zeros = false;
        for (auto x: p) {
            zeros |= x == 0;
        }
        if (!zeros) {
            for (auto q: p) {
                printf("%.16Lf ", q);
            }
            printf("\n");
            return p;
        }
    }
    return p;
}

bool get_all(json const &jin, int last_elems, vector<money> & price_vector, mapped_file &mf) {
    // 0 - usdt
    // 1 - btc
    // 2 - eth

    std::vector<string> names{"btcusdt", "ethusdt", "ethbtc"};
    vector<pair<int,int>> pairs{{0,1}, {0,2}, {1, 2}};
    auto d0 = get_data(names[0]);
    auto d1 = get_data(names[1]);
    auto d2 = get_data(names[2]);
    map<string, vector<trade_data>> all_trades{{names[0], d0}, {names[1], d1}, {names[2], d2}};
    u64 min_time = 1ull << 63;
    u64 max_time = 0;
    for (auto const &t: all_trades) {
        min_time = min(min_time, t.second.front().t);
        max_time = max(max_time, t.second.back().t);
    }
    vector<trade_one> out;

    for (size_t i = 0; i < names.size(); i++) {
        auto &trades = all_trades[names[i]];
        for (auto &trade: trades) {
            if (trade.t >= min_time && trade.t <= max_time) {
                trade.pair1 = pairs[i];
                out.push_back({trade.t + (trade.pair1.first + trade.pair1.second) * 15, trade});
            }
        }
    }
    clock_t start_time = clock();
    //debug_print("out first 5", out, 5);
    //debug_print("out last 5", out, -5);
    sort(out.begin(), out.end(), [](trade_one const &l, trade_one const &r) {
        return l.t < r.t;
    });
    clock_t end_time = clock();
    //debug_print("sorted out first 5", out, 5);
    //debug_print("sorted out last 5", out, -5);
    vector<trade_data> ret;
    for (auto &q: out) {
        ret.push_back(q.trade);
    }
    //printf("total %zu elements\n", ret.size());
    print_clock("sorting took", start_time, end_time);
    if (last_elems > 0) {
        printf("Trimming: use last %d elements\n", last_elems);
        ret.erase(ret.begin(), ret.begin() + ret.size() - last_elems);
    }
    price_vector = get_price_vector(3, ret);
    string tmp_name = "_tmp." + std::to_string(getpid());
    FILE *tmp = fopen(tmp_name.c_str(), "w+");
    if (tmp == nullptr) {
        printf("Temp file '%s' is not available\n", tmp_name.c_str());
        return false;
    }
    printf("Using temp file '%s' as interprocedural connect\n", tmp_name.c_str());
    fwrite(&ret[0], sizeof (trade_data), ret.size(), tmp);
    //unlink(tmp_name.c_str()); // Does not work in Windows
    mf.map(tmp_name);
    return true;
}

money geometric_mean(money const *x, size_t N) {
    switch (N) {
        case 2:
            return sqrtl(x[0] * x[1]);
        case 3:
            // Newton process should converged without sort
            //sort(x.begin(),x.end(), [](money const &l, money const &r) {
            //    return l > r;
            //});
            // {0,1,2} {0,2,1} {1,0,2} {1,2,0} {2,0,1} {2,1,0}
            auto median = [&x] {
                money D = x[0];
                if (x[0] >= x[1]) {
                    // {1,0,2} {2,0,1} {2,1,0}
                    if (x[0] >= x[2]) {
                        // {2,0,1} {2,1,0}
                        if (x[1] >= x[2]) D = x[2];
                        else D = x[1];
                    } // else {1,0,2}
                } else {
                    // {0,1,2} {0,2,1} {1,2,0}
                    if (x[0] < x[2]) {
                        // {0,1,2} {0,2,1}
                        if (x[1] >= x[2]) D = x[2];
                        else D = x[1];
                    } // else {1,2,0}
                }
            };
            auto min_max_mean = [&x] {
                if (x[0] >= x[1]) {
                    // {1,0,2} {2,0,1} {2,1,0}
                    if (x[0] >= x[2]) {
                        // {2,0,1} {2,1,0}
                        if (x[1] >= x[2]) return sqrtl(x[0]*x[1]);
                        else return sqrtl(x[0]*x[2]);
                    }
                    return sqrtl(x[1]*x[2]);
                } else {
                    // {0,1,2} {0,2,1} {1,2,0}
                    if (x[0] < x[2]) {
                        // {0,1,2} {0,2,1}
                        if (x[1] >= x[2]) return sqrtl(x[0]*x[1]);
                        else return sqrtl(x[0]*x[2]);
                    }
                    return sqrtl(x[1]*x[2]); // else {1,2,0}
                }

            };
            //money MAX = max(x[0], max(x[1], x[2]));
            money prod = x[0] * x[1] * x[2];
            //money MIN = min(x[0], min(x[1], x[2]));
            //money D = sqrtl(MAX * MIN);
            auto D = min_max_mean();
            for (int i = 0; i < 255; i++) {
                money D_prev = D;
                D = (D + D + prod / D / D) *0.333333333333333333333333L;
                auto diff = mabs(D - D_prev);
                if (diff <= 1E-18 or diff * 1E18L < D) {
                    return D;
                }
            }
            break;
    }
    throw std::logic_error("geometric_mean: Did not converge");
}

auto reduction_coefficient(money const *x, size_t N, money gamma) {
    money K = 1.L;
    money S = 0.L;
    for (size_t i = 0; i < N; i++) S += x[i]; // = sum(x)
    for (size_t i = 0; i < N; i++)  {
        K *= N * x[i] / S;
    }
    if (gamma > 0) {
        K = gamma / (gamma + 1.L - K);
    }
    return K;
}

void print(money const *x, size_t N) {
    printf("[");
    for (size_t i = 0; i < N; i++) {
        printf("%.16Lf ", x[i]);
    }
    printf("]\n");
}

auto newton_D(money A, money gamma, money const *xx, size_t N, money D0) {
    money D = D0;
    money S = 0;
    money x[MAX_ARRAY];
    for (size_t i = 0; i < N; i++) {
        S += x[i] = xx[i];
    }
    sort(x+0, x+N, [](money l, money r) { return l > r; });
    auto NN = 1.L;
    for (size_t j = 0; j < N; j++) { // XXX or just set A to be A*N**N?
        NN *= N;
    }
    A *= NN;
    //for (size_t j = 0; j < N; j++) { // XXX or just set A to be A*N**N?
    //    A = A * N;
    //}

    for (int i = 0; i < 255; i++) {
        money D_prev = D;

        money K0 = NN;
        for (auto const &_x: x) {
            K0 = K0 * _x / D;
        }

        money _g1k0 = mabs((gamma + 1.L - K0));

        // # D / (A * N**N) * _g1k0**2 / gamma**2
        money mul1 = D / gamma * _g1k0 / gamma * _g1k0 / A;

        // # 2*N*K0 / _g1k0
        money mul2 = 2.L * N * K0 / _g1k0;

        money neg_fprime = (S + S * mul2) + mul1 * N / K0 - mul2 * D;
        assert (neg_fprime > 0); //   # Python only: -f' > 0

        // # D -= f / fprime
        D = (D * neg_fprime + D * S - D * D) / neg_fprime - D * (mul1 / neg_fprime) * (1.L - K0) / K0;

        if (D < 0) {
            D = mabs(D) / 2.L;
        }
        if (mabs(D - D_prev) <= max(1e-16L, D / 1e14L)) {
            return D;
        }
    }
    throw std::logic_error("Newton_D: did not converge");
}

auto newton_D_3(money A, money gamma, money const *xx, money D0) {
    money D = D0;
    money S = 0;
    money x[3];
    for (size_t i = 0; i < 3; i++) {
        S += x[i] = xx[i];
    }
    if (x[0] < x[1]) std::swap(x[0],x[1]);
    if (x[1] < x[2]) std::swap(x[1],x[2]);
    if (x[0] < x[1]) std::swap(x[0],x[1]);
    // sort(x+0, x+3, [](money l, money r) { return l > r; });
    A *= 27.L;
    //for (size_t j = 0; j < N; j++) { // XXX or just set A to be A*N**N?
    //    A = A * N;
    //}

    for (int i = 0; i < 255; i++) {
        money D_prev = D;

        money K0 = 27.L;
        for (auto const &_x: x) {
            K0 = K0 * _x / D;
        }

        money _g1k0 = mabs((gamma + 1.L - K0));

        // # D / (A * N**N) * _g1k0**2 / gamma**2
        money mul1 = D / gamma * _g1k0 / gamma * _g1k0 / A;

        // # 2*N*K0 / _g1k0
        money mul2 = 2.L * 3.L * K0 / _g1k0;

        money neg_fprime = (S + S * mul2) + mul1 * 3.L / K0 - mul2 * D;
        assert (neg_fprime > 0); //   # Python only: -f' > 0

        // # D -= f / fprime
        D = (D * neg_fprime + D * S - D * D) / neg_fprime - D * (mul1 / neg_fprime) * (1.L - K0) / K0;

        if (D < 0) {
            D = mabs(D) / 2.L;
        }
        if (mabs(D - D_prev) <= max(1e-16L, D / 1e14L)) {
            return D;
        }
    }
    throw std::logic_error("Newton_D: did not converge");
}

auto newton_y(money A, money gamma, money const *x, size_t N, money D, int i) {
    money y = D / N;
    money K0_i = 1.;
    money S_i = 0.;
    money x_sorted[N - 1];
    for (int j = 0, cnt = 0; j < N; j++) {
        if (j != i) x_sorted[cnt++] = x[j];
    }
    sort(x_sorted+0, x_sorted+N-1);
    money max_x_sorted = x_sorted[N-2];
    money convergence_limit = max(max_x_sorted / 1E14L, D / 1E14L);
    convergence_limit = max(convergence_limit, 1E-16L);
    for (auto const &_x: x_sorted) {
        y = y * D / (_x * N); //  # Small _x first
        S_i += _x;
        K0_i = K0_i * _x * N / D; //  # Large _x first
    }
    auto NN = 1.;
    for (size_t j = 0; j < N; j++) { // in range(N):  # XXX or just set A to be A*N**N?
        NN *= N;
    }
    A = A * NN;
    auto g2a = (gamma * gamma * A);
    for (size_t j = 0; j < 255; j++) {
        money y_prev = y;

        money K0 = K0_i * y * N / D;
        money K0_1 = 1.L - K0;
        money S = S_i + y;

        money _g1k0 = mabs((gamma + K0_1));

        // D / (A * N**N) * _g1k0**2 / gamma**2
        //money mul1 = D / gamma * _g1k0 / gamma * _g1k0 / A;
        money mul1 = D * _g1k0 * _g1k0 / g2a;
        // 2*K0 / _g1k0
        money mul2 = 1.L + (K0 + K0) / _g1k0;

        // money yfprime = y + S * mul2 + mul1 - D * mul2;
        money yfprime = y + mul1 + (S - D) * mul2;
        money fprime = yfprime / y;
        assert (fprime  > 0) ;  //# Python only: f' > 0

        // y -= f / f_prime;  y = (y * fprime - f) / fprime
        // y = (yfprime + D - S) / fprime + mul1 / fprime * K0_1 / K0;
        y = ((yfprime + D - S) + mul1 * K0_1 / K0) / fprime;
        if (j > 100) { //  # Just logging when doesn't converge
            printf("%zu %.6Lf %.16Lf ", j, y, D);
            print(x, N);
        }
        if (y < 0 or fprime < 0) {
            y = y_prev / 2.L;
        }

        if (mabs(y - y_prev) <= max(convergence_limit, y / 1e14L)) {
            return y;
        }
    }
    throw std::logic_error("Did not converge");
}

auto newton_y_3(money A, money gamma, money const *x, money D, int i) {
    money y = D / 3.;
    money K0_i = 1.;
    money S_i = 0.;
    money x_sorted[2];
    switch (i) {
        case 0:
            if (x[1] < x[2]) {
                x_sorted[0] = x[1];
                x_sorted[1] = x[2];
            } else {
                x_sorted[1] = x[1];
                x_sorted[0] = x[2];
            }
            break;
        case 1:
            if (x[0] < x[2]) {
                x_sorted[0] = x[0];
                x_sorted[1] = x[2];
            } else {
                x_sorted[1] = x[0];
                x_sorted[0] = x[2];
            }
            break;
        case 2:
            if (x[0] < x[1]) {
                x_sorted[0] = x[0];
                x_sorted[1] = x[1];
            } else {
                x_sorted[1] = x[0];
                x_sorted[0] = x[1];
            }
            break;
        default:
            abort();
    }
    money max_x_sorted = x_sorted[1];
    money convergence_limit = max(max_x_sorted / 1E14L, D / 1E14L);
    convergence_limit = max(convergence_limit, 1E-16L);
    y = y * D / (x_sorted[0] * 3);
    S_i += x_sorted[0];
    K0_i = K0_i * x_sorted[0] * 3 / D;
    y = y * D / (x_sorted[1] * 3);
    S_i += x_sorted[1];
    K0_i = K0_i * x_sorted[1] * 3 / D;
    A = A * 27.L;
    auto g2a = (gamma * gamma * A);
    for (size_t j = 0; j < 255; j++) {
        money y_prev = y;

        money K0 = K0_i * y * 3 / D;
        money K0_1 = 1.L - K0;
        money S = S_i + y;

        money _g1k0 = mabs((gamma + K0_1));

        // D / (A * N**N) * _g1k0**2 / gamma**2
        //money mul1 = D / gamma * _g1k0 / gamma * _g1k0 / A;
        money mul1 = D * _g1k0 * _g1k0 / g2a;
        // 2*K0 / _g1k0
        money mul2 = 1.L + (K0 + K0) / _g1k0;

        // money yfprime = y + S * mul2 + mul1 - D * mul2;
        money yfprime = y + mul1 + (S - D) * mul2;
        money fprime = yfprime / y;
        assert (fprime  > 0) ;  //# Python only: f' > 0

        // y -= f / f_prime;  y = (y * fprime - f) / fprime
        // y = (yfprime + D - S) / fprime + mul1 / fprime * K0_1 / K0;
        y = ((yfprime + D - S) + mul1 * K0_1 / K0) / fprime;
        if (j > 100) { //  # Just logging when doesn't converge
            printf("%zu %.6Lf %.16Lf ", j, y, D);
            print(x, 3);
        }
        if (y < 0 or fprime < 0) {
            y = y_prev * 0.5L;
        }

        if (mabs(y - y_prev) <= max(convergence_limit, y * 1e-14L)) {
            return y;
        }
    }
    throw std::logic_error("Did not converge");
}

money solve_x(money A, money gamma, money const *x, size_t N, money D, int i) {
    if (N == 3) return newton_y_3(A, gamma, x, D, i);
    else return newton_y(A, gamma, x, N, D, i);
}

auto solve_D(money A, money gamma, money const *x, size_t N) {
    auto D0 = N * geometric_mean(x, N); //  # <- fuzz to make sure it's ok XXX
    return N == 3 ? newton_D_3(A, gamma, x, D0) : newton_D(A, gamma, x, N, D0);
}

struct Curve {
    Curve(money A, money gamma, money D, int n, vector<money> const &p) {
        this->A = A;
        this->gamma = gamma;
        this->n = n;
        if (!p.empty()) {
            this->p = p;
        } else {
            this->p.resize(n, 1.L);
        }
        this->x.resize(n);
        for(int i = 0; i < n; i++) {
            x[i] = D / n / p[i];
        }
    }

    auto xp(money *ret, size_t N) const {
        for (int i = 0; i < N; i++) {
            ret[i] = x[i] * p[i];
            assert(x[i] > 0);
        }
    }

    auto D() const {
        money xp[MAX_ARRAY];
        this->xp(xp,n);
#if 0
        for (size_t i = 0; i < n; i++) {
            if (xp[i] <= 0) {
                throw std::logic_error("Curve::D(): x <= 0");
            }
        }
#endif
        auto ret = solve_D(A, gamma, xp, n);
        return ret;
    }

    money y(money x, int i, int j) {
        money xp[MAX_ARRAY];
        this->xp(xp,n);
        xp[i] = x * this->p[i];
        auto yp = solve_x(A, gamma, xp, n, this->D(), j);
        auto ret = yp / this->p[j];
        return ret;
    }

    money A;
    money gamma;
    size_t n;
    vector<money> p;
    vector<money> x;

};


struct Trader {
    Trader(money A, money gamma, money D, int n, vector<money> const &p0,
           money mid_fee,
           money out_fee,
           money price_threshold,
           money const &fee_gamma,
           money adjustment_step,
           int ma_half_time,
           bool log = true) : curve(A, gamma, D, n, p0) {
        this->p0 = p0;
        this->price_oracle = this->p0;
        this->last_price = this->p0;
        // this->curve = Curve(A, gamma, D, n, p0);
        this->dx = D * 1e-8L;
        this->mid_fee = mid_fee;
        this->out_fee = out_fee;
        this->D0 = this->curve.D();
        this->xcp_0 = this->get_xcp();
        this->xcp_profit = 1.L;
        this->xcp_profit_real = 1.L;
        this->xcp = this->xcp_0;
        this->price_threshold = price_threshold;
        this->adjustment_step = adjustment_step;
        this->log = log;
        this->fee_gamma = fee_gamma; // || gamma;
        this->total_vol = 0.0;
        this->ma_half_time = ma_half_time;
        this->ext_fee = 0; //   # 0.03e-2
        this->slippage = 0;
        this->slippage_count = 0;
        this->not_adjusted = false;
        this->heavy_tx = 0;
        this->light_tx = 0;
        this->is_light = false;
        this->t = 0;
    }


    auto fee(size_t N) {
        money xp[MAX_ARRAY];
        curve.xp(xp, N);
        auto f = reduction_coefficient(xp, N, fee_gamma);
        return (mid_fee * f + out_fee * (1.L - f));
    }

    money get_xcp() const {
        // First calculate the ideal balance
        //  Then calculate, what the constant-product would be
        auto D = curve.D();
        size_t N = curve.x.size();
        money X[MAX_ARRAY];
        for (size_t i = 0; i < N; i++) {
            X[i] = D  / (N * curve.p[i]);
        }
        return geometric_mean(X, N);
    }

    auto price(int i, int j) {
        auto dx_raw = dx  / curve.p[i];
        auto curve_res = curve.y(curve.x[i] + dx_raw, i, j);
        auto ret = dx_raw  / (curve.x[j] - curve_res);
        return ret;
    }

    auto step_for_price(money dp, pair<int, int> p, int sign) {
        auto p0 = price(p.first, p.second);
        dp = p0 * dp;
        money x0[curve.x.size()];
        copy_money(x0, &curve.x[0], curve.x.size());
        auto step = dx / curve.p[p.first];
        while (true) {
            curve.x[p.first] = x0[p.first] + sign * step;
            auto dp_ = mabs(p0 - price(p.first, p.second));
            if (dp_ >= dp or step >= curve.x[p.first] / 10.L) {
                copy_money(&curve.x[0], x0, curve.x.size());
                return step;
            }
            step += step;
        }
    }

    void update_xcp(bool only_real=false) {
        auto xcp = get_xcp();
        xcp_profit_real = xcp_profit_real * xcp / this->xcp;
        if (not only_real) {
            xcp_profit = xcp_profit * xcp / this->xcp;
        }
        this->xcp = xcp;
    }

    inline void copy_money(money *to, money const *from, int n) {
        ::memcpy(to, from, sizeof(money) * n);
    }
    money buy(money dx, int i, int j, money max_price=1e100L) {
        //"""
        //Buy y for x
        //"""
        try {
            money x_old[MAX_ARRAY];
            copy_money(x_old, &curve.x[0], curve.x.size());
            auto x = curve.x[i] + dx;
            auto y = curve.y(x, i, j);
            auto dy = curve.x[j] - y;
            curve.x[i] = x;
            curve.x[j] = y;
            auto fee = this->fee(curve.x.size());
            curve.x[j] += dy * fee;
            dy = dy * (1.L - fee);
            if ((dx / dy) > max_price or dy < 0) {
                copy_money(&curve.x[0], x_old, curve.x.size());
                return 0;
            }
            update_xcp();
            return dy;
        } catch (...) {
            return 0;
        }
    }

    money sell(money dy, int i, int j, money min_price=0) {
        // """
        // Sell y for x
        // """
        try {
            money x_old[MAX_ARRAY];
            copy_money(x_old, &curve.x[0], curve.x.size());
            auto y = curve.x[j] + dy;
            auto x = curve.y(y, j, i);
            auto dx = curve.x[i] - x;
            curve.x[i] = x;
            curve.x[j] = y;
            auto fee = this->fee(curve.x.size());
            curve.x[i] += dx * fee;
            dx = dx * (1.L - fee);
            if ((dx / dy) < min_price or dx < 0) {
                copy_money(&curve.x[0], x_old, curve.x.size());
                return 0;
            }
            update_xcp();
            return dx;
        } catch (...) {
            return 0;
        }
    }

    void ma_recorder(u64 t, vector<money> const &price_vector) {
        //  XXX what if every block only has p_b being last
        if (t > this->t) {
            money alpha = powl(0.5, ((money)(t - this->t) / this->ma_half_time));
            for (int k = 1; k <= 2; k++) {
                price_oracle[k] = price_vector[k] * (1 - alpha) + price_oracle[k] * alpha;
            }
            this->t = t;
        }
    }

    auto tweak_price(u64 t, int a, int b, money p) {
        ma_recorder(t, last_price);
        if (b > 0) {
            last_price[b] = p * last_price[a];
        } else {
            last_price[a] = last_price[0] / p;
        }

        // # price_oracle looks like [1, p1, p2, ...] normalized to 1e18
        money S = 0;
        for (size_t i = 0; i < price_oracle.size(); i++) {
            auto t = price_oracle[i] / curve.p[i] - 1.L;
            S += t*t;
        }
        auto norm = S;
        auto mxp = (max(price_threshold, adjustment_step));
        norm = sqrt(norm); // .root_to();
        if (norm <= mxp) {
            // Already close to the target price
            is_light = true;
            light_tx += 1;
            return norm;
        }
        if (not not_adjusted and (xcp_profit_real - 1.L > (xcp_profit - 1.L) / 2.L + 1e-5L)) {
            not_adjusted = true;
        }
        if (not not_adjusted) {
            light_tx += 1;
            is_light = true;
            return norm;
        }
        heavy_tx += 1;
        is_light = false;

        vector<money> p_new(price_oracle.size());
        p_new[0] = 1.L;
        for (size_t i = 1; i < price_oracle.size(); i++) {
            auto p_target = curve.p[i];
            auto p_real = price_oracle[i];
            p_new[i] = p_target + adjustment_step * (p_real - p_target) / norm;
        }
        auto old_p = curve.p;
        auto old_profit = xcp_profit_real;
        auto old_xcp = xcp;

        curve.p = p_new;
        update_xcp(true);

        if (2.L * (xcp_profit_real - 1.L) <= (xcp_profit - 1.L)) {
            //  If real profit is less than half of maximum - revert params back
            curve.p = old_p;
            xcp_profit_real = old_profit;
            xcp = old_xcp;
            not_adjusted = false;
            auto val = ((xcp_profit_real - 1.L - (xcp_profit - 1.L) / 2.L));
            // printf("%.10Lf\n", val);
        }
        return norm;
    }


    void simulate(mapped_file const &in) {
        // vector<trade_data> const &mdata
        const money CANDLE_VARIATIVES = 50;
        map<pair<int,int>,money> lasts;
        u64 start_t = 0;
        long double last_time = 0;
        size_t total_elements = in.size / sizeof(trade_data);
        trade_data const *mapped_data = (trade_data const *)in.base;
        trade_data const *mapped_data_ptr = mapped_data;
        for (size_t i = 0; i < total_elements; i++)  {
            // if (i > 10) abort();
            trade_data d = *mapped_data_ptr++;
            if (i == 0) start_t = d.t;
            if (last_time > 0) {
                last_time = d.t - last_time;
            }
            auto a = d.pair1.first;
            auto b = d.pair1.second;
            money vol{0.L};
            auto ext_vol = money(d.volume * price_oracle[b]); //  <- now all is in USD
            int ctr{0};
            money last;
            auto itl = lasts.find({a,b});
            if (itl == lasts.end()) {
                last = price_oracle[b] / price_oracle[a];
            } else {
                last = itl->second;
            }
            auto _high = last;
            auto _low = last;

            //  Dynamic step
            //  f = reduction_coefficient(self.curve.xp(), self.curve.gamma)
            auto candle = min(mabs((d.high - d.low) / d.high), 0.1L);
            candle = max(0.001L, candle);
            auto step1 = step_for_price(candle / CANDLE_VARIATIVES, d.pair1, 1);
            auto step2 = step_for_price(candle / CANDLE_VARIATIVES, d.pair1, -1);
            auto step = min(step1, step2);
            auto max_price = d.high;
            money _dx = 0;
            auto p_before = price(a, b);
            while (last < max_price and vol < ext_vol / 2.L) {
                auto dy = buy(step, a, b, max_price);
                if (dy == 0) {
                    break;
                }
                vol += dy * price_oracle[b];
                _dx += dy;
                last = step / dy;
                max_price = d.high;
                ctr += 1;
            }
            auto p_after = price(a, b);
            if (p_before != p_after) {
                slippage_count += last_time;
                slippage += last_time * _dx * (p_before + p_after) / (2.L * mabs(p_before - p_after) * curve.x[b]);
            }
            _high = last;
            auto min_price = d.low;
            _dx = 0;
            p_before = p_after;
            money prev_vol = vol;
            while (last > min_price and vol < ext_vol / 2.L) {
                auto dx = step / last;
                auto dy = sell(dx, a, b, min_price);
                if (dy == 0) {
                    break;
                }
                _dx += dx;
                vol += dx * price_oracle[b];
                last = dy / dx;
                min_price = d.low;
                ctr += 1;
            }
            p_after = price(a, b);
            if (p_before != p_after) {
                slippage_count += last_time;
                slippage += last_time * _dx * (p_before + p_after) / (2.L * mabs(p_before - p_after) * curve.x[b]);
            }
            _low = last;
            lasts[d.pair1] = last;

            tweak_price(d.t, a, b, (_high + _low) / 2.L);

            total_vol += vol;
            last_time = d.t;
            if (i % 1024 == 0 && log) {
                try {
                    long double last01, last02;
                    auto it01 = lasts.find({0,1});
                    if (it01 == lasts.end()) {
                        last01 = price_oracle[1] / price_oracle[0];
                    } else {
                        last01 = it01->second;
                    }
                    auto it02 = lasts.find({0,2});
                    if (it02 == lasts.end()) {
                        last02 = price_oracle[2] / price_oracle[0];
                    } else {
                        last02 = it02->second;
                    }
                    long double ARU_x = xcp_profit_real;
                    long double ARU_y = (86400.L * 365.L / (d.t - start_t + 1.L));
                    APY = powl(ARU_x, ARU_y) - 1.L;
                    printf("t=%llu %.1Lf%%\ttrades: %d\t"
                           "AMM: %.0Lf, %0.Lf\tTarget: %.0Lf, %.0Lf\t"
                           "Vol: %.4Lf\tPR:%.2Lf\txCP-growth: {%.5Lf}\t"
                           "APY:%.1Lf%%\tfee:%.3Lf%% %c\n",
                           d.t,
                           100.L * i / total_elements, ctr, last01, last02,
                           curve.p[1],
                           curve.p[2],
                           total_vol,
                           (xcp_profit_real - 1.) / (xcp_profit - 1.L),
                           xcp_profit_real,
                           APY * 100.L,
                           fee(curve.p.size()) * 100.L,
                           is_light? '*' : '.');
                } catch (std::exception const &e) {
                    printf("caught '%s'\n", e.what());
                }
            }
        }
    }



    vector<money> p0;
    vector<money> price_oracle;
    vector<money> last_price;
    u64 t;
    money dx;
    money mid_fee;
    money out_fee;
    money D0;
    money xcp, xcp_0;
    money xcp_profit;
    money xcp_profit_real;
    money price_threshold;
    money adjustment_step;
    bool log;
    money fee_gamma;
    money total_vol;
    int ma_half_time;
    money ext_fee;
    money slippage;
    money slippage_count;
    long double APY;
    bool not_adjusted;
    int  heavy_tx;
    int  light_tx;
    bool is_light;
    Curve curve;
};

static bool json_load(string const &name, json &j) {
    try {
        std::ifstream ifl(name);
        if (!ifl) throw std::logic_error("can't open file " + name);
        ifl >> j;
    } catch (std::exception const &ex) {
        printf("json_load: %s\n", ex.what());
        return false;
    }
    return true;
}

static bool json_save(string const &name, json const &j) {
    try {
        std::ofstream ofl(name);
        if (!ofl) throw std::logic_error("can't create file " + name);
        ofl << std::setw(4) << j << "\n";
    } catch (std::exception const &ex) {
        printf("json_load: %s\n", ex.what());
        return false;
    }
    return true;
}

int main(int argc, char **argv) {
    if (argc == 1) {
        printf("Usage: %s [prepare|simulate] in-json-file out-json-file\n", argv[0]);
        return 0;
    }
    string in_json_name = "sample_in.json";
    string out_json_name = "sample_out.json";
#if 0
    json jout;
    jout["datafile"][0] = "btcusdt";
    jout["datafile"][1] = "ethbtc";
    jout["datafile"][2] = "ethusdt";
    jout["configuration"][0]["A"] = 135;
    jout["configuration"][0]["gamma"] = 7e-5;
    jout["configuration"][0]["D"] = 100'000'000.;
    jout["configuration"][0]["mod_fee"] = 4e-4;
    jout["configuration"][0]["out_frr"] = 7e-5;
    jout["configuration"][0]["price_threshold"] = 0.0028;
    jout["configuration"][0]["fee_gamma"] = 0.01;
    jout["configuration"][0]["adjustment_step"] = 0.0015;
    jout["configuration"][0]["ma_half_time"] = 600;
    jout["configuration"][1]["A"] = 135;
    jout["configuration"][1]["gamma"] = 6e-5;
    jout["configuration"][1]["D"] = 100'000'000.;
    jout["configuration"][1]["mod_fee"] = 3e-4;
    jout["configuration"][1]["out_frr"] = 8e-5;
    jout["configuration"][1]["price_threshold"] = 0.0026;
    jout["configuration"][1]["fee_gamma"] = 0.01;
    jout["configuration"][1]["adjustment_step"] = 0.0013;
    jout["configuration"][1]["ma_half_time"] = 600;
    jout["debug"] = 0;
    json_save("sample_in.json", jout);
    return 0;
#endif
#if 1
    json jin;
    if (!json_load("sample_in.json", jin)) {
        return 0;
    }
    printf("Total configurations: %zu\n", jin["configuration"].size());
    //for (auto const &cfg: jin["configuration"]) {
    //    std::cout << std::setw(4) << cfg << "\n";
    // }
#endif
    int LAST_ELEMS = 0;
    clock_t start = clock();
    if (argc > 1 && string(argv[1]) == "trim") {
        LAST_ELEMS = 100000;
        if (argc > 2) LAST_ELEMS = atoi(argv[2]);
    }
    vector<money> price_vector;
    mapped_file test_data;
    if (!get_all(jin, LAST_ELEMS, price_vector, test_data)) {
        return 0;
    }
    //debug_print("test_data first 5", test_data, 5);
    //debug_print("test_data last 5", test_data, -5);
    Trader trader(135, (7e-5), money(100'000'000), 3, price_vector,
                  4e-4, 4.0e-3,
                  0.0028, 0.01L,
                  0.0015, 600);
    clock_t start_simulation = clock();
    printf("Begin simulation\n");
    unlink(test_data.name.c_str()); // Temp file can be deleted in *nix even being open
    trader.simulate(test_data);
    printf("Liquidity density vs that of xyz=k: %f\n", 2 * (double)(trader.slippage) / (double)(trader.slippage_count));
    printf("APY: %Lf%%\n", (trader.APY * 100.L));
    clock_t end = clock();
    print_clock("Total simulation time", start_simulation, end);
}