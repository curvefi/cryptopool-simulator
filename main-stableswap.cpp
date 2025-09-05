#include <cassert>
#include <cstdio>
#include <map>
#include <string>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <utility>
#include <sys/fcntl.h>
#include <sys/mman.h>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "json.hpp"
#include <queue>
#include <pthread.h>
#include <sys/time.h>
#include <sys/resource.h>
#ifndef MAP_NOCACHE
#define MAP_NOCACHE 0
#endif
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

#ifdef __MACH__
#include <mach/mach_init.h>
#include <mach/thread_act.h>
#include <mach/mach_port.h>
static double get_thread_time() {
    mach_port_t thread;
    kern_return_t kr;
    mach_msg_type_number_t count;
    thread_basic_info_data_t info;

    thread = mach_thread_self();

    count = THREAD_BASIC_INFO_COUNT;
    kr = thread_info(thread, THREAD_BASIC_INFO, (thread_info_t) &info, &count);
    double ret = 0;
    if (kr == KERN_SUCCESS && (info.flags & TH_FLAGS_IDLE) == 0) {
        ret += info.user_time.seconds;
        ret += info.user_time.microseconds / 1000000.;
        ret += info.system_time.seconds;
        ret += info.system_time.microseconds / 1000000.;
    }
    else {
        abort();
    }
    mach_port_deallocate(mach_task_self(), thread);
    return ret;
}

#else
static double get_thread_time() {
    struct rusage usage;
    getrusage (RUSAGE_THREAD, &usage);
    double ret = 0.;
    ret += usage.ru_utime.tv_sec;
    ret += usage.ru_utime.tv_usec / 1000000.;
    ret += usage.ru_stime.tv_sec;
    ret += usage.ru_stime.tv_usec / 1000000.;
    return ret;
}
#endif

static double get_total_time() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    double ret = 0.;
    ret += usage.ru_utime.tv_sec;
    ret += usage.ru_utime.tv_usec / 1000000.;
    ret += usage.ru_stime.tv_sec;
    ret += usage.ru_stime.tv_usec / 1000000.;
    return ret;
}

static double get_wall_time() {
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    double ret = 0.;
    ret += tv.tv_sec;
    ret += tv.tv_usec / 1000000.;
    return ret;
}

static void print_clock(string const &mesg, double start, double end) {
    printf("%s %.3lf sec\n", mesg.c_str(), double(end - start));
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
    auto start_time = get_thread_time();
    auto name_to_open = "download/" + fname + ".json";
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
        while (*p != '"' && *p!= ' ' && *p != ',' && *p != ']') p++;
        while (*p == ',' || *p == ' ' || *p == '"') p++;
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
            p = scan_u64(p, &d.t);
            if (d.t > 10000000000) {
                d.t /= 1000;
            }
            p = scan_double(p, &d.open);
            p = scan_double(p, &d.high);
            p = scan_double(p, &d.low);
            if (d.high < d.low) {
                auto _high = d.low;
                d.low = d.high;
                d.high = _high;
            }
            p = scan_double(p, &d.close);
            p = scan_double(p, &d.volume);
            ret.push_back(d);
            while (*p != ']') p++;
            p++; // skip ']'
        } else p++;
    }
    auto end_time = get_thread_time();
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

    size_t N = jin["datafile"].size();
    if (N == 1) N = 2;
    assert(N == 2 || N == 3);
    vector<string> names;
    map<string, vector<trade_data>> all_trades;
    vector<pair<int, int>> pairs;
    if (N == 2) {
        pairs.push_back({0, 1});
        string name = jin["datafile"][0];
        auto d0 = get_data(name);
        all_trades[name] = d0;
        names.resize(1);
        printf("using file '%s'\n", name.c_str());
        names[0] = name;
    } else {
        pairs.push_back({0, 1});
        pairs.push_back({0, 2});
        pairs.push_back({1, 2});
        names.resize(3);
        names[0] = jin["datafile"][0];
        auto d0 = get_data(names[0]);
        names[1] = jin["datafile"][1];
        auto d1 = get_data(names[1]);
        names[2] = jin["datafile"][2];
        auto d2 = get_data(names[2]);
        all_trades[names[0]] = d0;
        all_trades[names[1]] = d1;
        all_trades[names[2]] = d2;
    }
    u64 min_time = 1ull << 63;
    u64 max_time = 0;
    for (auto name: names) {
        for (auto const &t: all_trades[name]) {
            min_time = min(min_time, t.t);
            max_time = max(max_time, t.t);
        }
    }
    vector<trade_one> out;

    for (size_t i = 0; i < names.size(); i++) {
        auto &trades = all_trades[names[i]];
        for (auto &trade: trades) {
            if (trade.t >= min_time && trade.t <= max_time) {
                trade.pair1 = pairs[i];

                trade_data trade_min;
                trade_data trade_max;
                // t, open, high, low, close, volume, pair1

                // (1, 2) min
                // (0, 2) min
                // (0, 1) min
                // (0, 1) max
                // (0, 2) max
                // (1, 2) max

                trade_min.t = trade.t - (trade.pair1.first + trade.pair1.second) * 10 + 10;
                trade_max.t = trade.t + (trade.pair1.first + trade.pair1.second) * 10 - 0;
                trade_min.open = trade.open;
                trade_max.close = trade.close;
                trade_min.pair1 = trade.pair1;
                trade_max.pair1 = trade.pair1;
                trade_min.volume = trade.volume;
                trade_max.volume = trade.volume;

                if (mabs(trade.open - trade.low) + mabs(trade.close - trade.high) < mabs(trade.open - trade.high) + mabs(trade.close - trade.low)) {
                    // shortest path: open -> low -> high -> close
                    trade_min.high = trade.low;
                    trade_min.low = trade.low;
                    trade_min.close = trade.low;
                    trade_max.open = trade.high;
                    trade_max.high = trade.high;
                    trade_max.low = trade.high;
                } else {
                    // shortest path: open -> high -> low -> close
                    trade_min.high = trade.high;
                    trade_min.low = trade.high;
                    trade_min.close = trade.high;
                    trade_max.open = trade.low;
                    trade_max.high = trade.low;
                    trade_max.low = trade.low;
                }

                    
                if (trade.t < 1672643800*10) { // TS HARDSTOP
                    out.push_back({trade.t - (trade.pair1.first + trade.pair1.second) * 10 + 10, trade_min});
                    out.push_back({trade.t + (trade.pair1.first + trade.pair1.second) * 10 - 0, trade_max});
                }
            }
        }
    }
    auto start_time = get_thread_time();
    sort(out.begin(), out.end(), [](trade_one const &l, trade_one const &r) {
        return l.t < r.t;
    });
    auto end_time = get_thread_time();
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
    price_vector = get_price_vector(N, ret);
    string tmp_name = "_tmp." + std::to_string(getpid());
    FILE *tmp = fopen(tmp_name.c_str(), "w+");
    if (tmp == nullptr) {
        printf("Temp file '%s' is not available\n", tmp_name.c_str());
        return false;
    }
    printf("Using temp file '%s' as interprocedural connect\n", tmp_name.c_str());
    fwrite(&ret[0], sizeof(trade_data), ret.size(), tmp);

    fflush(tmp);
    // int fd = fileno(tmp);
    // if (fd >= 0) fsync(fd);
    fclose(tmp);
    //unlink(tmp_name.c_str()); // Does not work in Windows
    mf.map(tmp_name);
    return true;

}

money geometric_mean_2(money const *x) {
    return sqrtl(x[0] * x[1]);
}

money geometric_mean_3(money const *x) {
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
        if (diff <= 1E-12 or diff * 1E12L < D) {
            return D;
        }
    }
    throw std::logic_error("geometric_mean: Did not converge");
}


auto reduction_coefficient_3(money const *x, money gamma) {
    money K = 27.L;
    money S = x[0] + x[1] + x[2];
    for (size_t i = 0; i < 3; i++)  {
        K *= x[i] / S;
    }
    if (gamma > 0) {
        K = gamma * K / (gamma * K + 1.L - K);
    }
    return K;
}

auto reduction_coefficient_2(money const *x, money gamma) {
    money K = 1.L;
    money S = 0.L;
    for (size_t i = 0; i < 2; i++) S += x[i]; // = sum(x)
    for (size_t i = 0; i < 2; i++)  {
        K *= 2 * x[i] / S;
    }
    if (gamma > 0) {
        K = gamma * K / (gamma * K + 1.L - K);
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

auto newton_D_2(money A, money gamma, money const *xx, money D0) {
    // ***
    // This now uses stableswap invariant (because invariants are pluggable)
    // ***
    money S = 0;
    const size_t N = 2;
    money x[2];
    for (size_t i = 0; i < N; i++) {
        S += x[i] = xx[i];
    }
    money D = D0;
    A *= N;  // A is already A * N**(N-1)

    for (int i=0; i < 255; i++) {
        money D_P = D;
        for (auto const &_x: x) {
            D_P = D_P * D / (N * _x);
        }
        money Dprev = D;
        D = (A * S + D_P * N) * D / ((A - 1) * D + (N + 1) * D_P);
        if (mabs(D - Dprev) <= 1e-6) {
            return D;
        }
    }
    return D; // we ignore convergence error in simulation
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
    money rev_gamma = 1.L / gamma;
    money gamma_1 = 1.L + gamma;

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
        D = (D * (neg_fprime + S - D)) / neg_fprime - D * (mul1 / neg_fprime) * (1.L - K0) / K0;

        if (D < 0) {
            D = mabs(D) / 2.L;
        }
        if (mabs(D - D_prev) <= max(1e-16L, D / 1e14L)) {
            return D;
        }
    }
    return D;  // XXX
    throw std::logic_error("Newton_D: did not converge");
}

auto newton_y(money A, money gamma, money const *x, size_t N, money D, int i) {
    // ***
    // This now uses stableswap invariant (because invariants are pluggable)
    // ***
    A = A * N;
    int other = 1 - i;
    money c = D*D / (x[other] * N);
    c = c * D / (A * N);
    money b = x[other] + D / A - D;
    money y_prev = 0;
    money y = D;
    for (size_t k = 0; k < 255; k++) {
        y_prev = y;
        y = (y*y + c) / (2 * y + b);
        if (mabs(y - y_prev) <= 1) {
            return y;
        }
    }
    return y;  // XXX
    throw std::logic_error("Did not converge");
}
    // // get_p: Calculates dx/dy (needs to be multiplied by price_scale for actual value).
    // static T get_p(
    //     const std::array<T, Traits::N>& _xp,
    //     const T& _D,
    //     const std::array<T, Traits::N>& _A_gamma
    // ) {
    //     T ANN = _A_gamma[0] * Traits::N;
    //     T Dr  = _D / T(Traits::N * Traits::N); // D / N_COINS**N_COINS

    //     for (size_t idx = 0; idx < Traits::N; ++idx) {
    //         Dr = Dr * _D / _xp[idx];
    //     }

    //     T xp0_A = ANN * _xp[0] / Traits::A_MULTIPLIER();

    //     // For uint256, PRECISION() = 1e18; for double, PRECISION() = 1.0
    //     return (
    //         Traits::PRECISION() * (xp0_A + Dr * _xp[0] / _xp[1])
    //     ) / (xp0_A + Dr);
    // }


auto get_p_2(money const *x, money D, money A, money gamma) {
    money ANN = A * 2.;
    money Dr = D / 4.;
    for (size_t idx = 0; idx < 2; ++idx) {
        Dr = Dr * D / x[idx];
    }
    money xp0_A = ANN * x[0];
    return
        (xp0_A + Dr * x[0] / x[1]) / (xp0_A + Dr);
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
    money convergence_limit = max(max_x_sorted * 1E-14L, D * 1E-14L);
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
    if (N == 3) {
        auto D0 = 3 * geometric_mean_3(x); //  # <- fuzz to make sure it's ok XXX
        return newton_D_3(A, gamma, x, D0);
    } else {
        auto D0 = 2 * geometric_mean_2(x); //  # <- fuzz to make sure it's ok XXX
        return newton_D_2(A, gamma, x, D0);
    }
}

struct Curve {
    Curve(json const &jconf, vector<money> const &p) {
        this->A = jconf["A"];
        this->gamma = jconf["gamma"];
        money D = jconf["D"];
        this->n = jconf["n"];
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

    auto xp_3(money *ret) const {
        for (int i = 0; i < 3; i++) {
            ret[i] = x[i] * p[i];
            assert(x[i] > 0);
        }
    }

    auto xp_2(money *ret) const {
        for (int i = 0; i < 2; i++) {
            ret[i] = x[i] * p[i];
            assert(x[i] > 0);
        }
    }

    auto D_3() const {
        money xp[3];
        this->xp_3(xp);
        auto ret = solve_D(A, gamma, xp, n);
        return ret;
    }

    auto D_2() const {
        money xp[2];
        this->xp_2(xp);
        auto ret = solve_D(A, gamma, xp, n);
        return ret;
    }

    money y_3(money x, int i, int j) {
        money xp[3];
        this->xp_3(xp);
        xp[i] = x * this->p[i];
        auto yp = solve_x(A, gamma, xp, n, this->D_3(), j);
        auto ret = yp / this->p[j];
        return ret;
    }

    money y_2(money x, int i, int j) {
        money xp[2];
        this->xp_2(xp);
        xp[i] = x * this->p[i];
        auto yp = solve_x(A, gamma, xp, n, this->D_2(), j);
        auto ret = yp / this->p[j];
        return ret;
    }
    
    money p_2(int i, int j) {
        money xp[2];
        this->xp_2(xp);
        auto p = get_p_2(xp, this->D_2(), this->A, this->gamma);
        return p;
    }

    money A;
    money gamma;
    size_t n;
    vector<money> p;
    vector<money> x;

};

struct extra_data {
    money APY = 0;
    money APY_boost = 0;
    money xcp_profit_real = 0;
    money xcp_profit = 0;
    money vpminusone = 0;                 // APY from virtual price (alias of APY for clarity)
    money APY_coin0_raw = -1;         // Raw APY from TVL growth in coin0 units
    money APY_coin0_boost_raw = -1;   // Raw APY excluding donations (boost) in coin0
    money APY_coin0 = -1;             // Baseline-subtracted APY vs HODL
    money APY_coin0_boost = -1;       // Baseline-subtracted APY excluding donations
    money donation_coin0_total = 0;    // Total donation (boost) value measured in coin0 units
    money liq_density = 0;
    money slippage = 0;
    money volume = 0;
    money trade_volume = 0; // sum of external trade volumes (USD)
    money arb_profit_coin0 = 0; // sum of positive arbitrage profits in coin0 units
    unsigned long long n_rebalances = 0; // count of successful tweak_price adjustments
    unsigned long long n_trades = 0; // total number of executed trades
    // Time-averaged absolute difference between pool price_scale and CEX price (coin0 per coin1)
    money avg_cex_diff = 0;
    money max_cex_diff = 0;
    money cex_diff_l2 = 0;
    // Fraction of time |price_scale - p_cex|/p_cex > threshold
    money cex_follow_time_frac = 0;
};

struct simulation_data {
    int num = 0;
    json const *jconf = nullptr;
    vector<money> const *price_vector = nullptr;
    mapped_file const *test_data = nullptr;
    extra_data result;
    size_t total = 0;
    size_t current = 0;
};


struct Trader {
    Trader(json const &jconf, vector<money> const &p0) : curve(jconf, p0) {
        money A = jconf["A"];
        money gamma = jconf["gamma"];
        money D = jconf["D"];
        int n = jconf["n"];
        mid_fee = jconf["mid_fee"];
        out_fee = jconf["out_fee"];
        fee_gamma = jconf["fee_gamma"];
        adjustment_step = jconf["adjustment_step"];
        allowed_extra_profit = jconf["allowed_extra_profit"];
        ma_half_time = jconf["ma_half_time"];
        this->ext_fee = jconf["ext_fee"];
        this->gas_fee = jconf["gas_fee"];
        this->boost_rate = jconf["boost_rate"];
        this->boost_rate = this->boost_rate / (86400L * 365L);
        this->boost_integral = 1.L;
        log = jconf["log"];
        this->p0 = p0;
        this->price_oracle = this->p0;
        this->last_price = this->p0;
        // this->curve = Curve(A, gamma, D, n, p0);
        this->dx = D * 1e-10L;
        this->D0 = n == 3 ? this->curve.D_3() : this->curve.D_2();
        this->xcp_0 = n == 3 ? this->get_xcp_3() : this->get_xcp_2();
        this->xcp_profit = 1.L;
        this->xcp_profit_real = 1.L;
        this->xcp = this->xcp_0;
        this->total_vol = 0.0;
        this->volume = 0;
        this->not_adjusted = false;
        this->heavy_tx = 0;
        this->light_tx = 0;
        this->is_light = false;
        this->t = 0;
        this->n_rebalances = 0;
        this->n_trades = 0;
        // Deviation threshold for time-fraction metric (relative, e.g. 0.1 = 10%)
        this->cex_diff_rel_thresh = 0.1L;
        auto it_thresh = jconf.find("cex_diff_rel_thresh");
        if (it_thresh != jconf.end()) {
            try { this->cex_diff_rel_thresh = it_thresh->get<long double>(); } catch (...) {}
        }

    }


    auto fee_3() {
        money xp[3];
        curve.xp_3(xp);
        auto f = reduction_coefficient_3(xp, fee_gamma);
        return (mid_fee * f + out_fee * (1.L - f));
    }

    auto fee_2() {
        money xp[2];
        curve.xp_2(xp);
        auto f = reduction_coefficient_2(xp, fee_gamma);
        return (mid_fee * f + out_fee * (1.L - f));
    }

    money get_xcp_3() const {
        // First calculate the ideal balance
        //  Then calculate, what the constant-product would be
        auto D = curve.D_3();
        size_t N = 3;
        money X[3];
        for (size_t i = 0; i < 3; i++) {
            X[i] = D  / (3 * curve.p[i]);
        }
        return geometric_mean_3(X);
    }

    money get_xcp_2() const {
        // First calculate the ideal balance
        //  Then calculate, what the constant-product would be
        auto D = curve.D_2();
        size_t N = 2;
        money X[2];
        for (size_t i = 0; i < 2; i++) {
            X[i] = D  / (2 * curve.p[i]);
        }
        return geometric_mean_2(X);
    }


    auto price_3(int i, int j) {
        auto dx_raw = dx  / curve.p[i];
        auto curve_res = curve.y_3(curve.x[i] + dx_raw, i, j);
        auto ret = dx_raw  / (curve.x[j] - curve_res);
        return ret;
    }

    auto price_2(int i, int j) {
        // auto dx_raw = dx  / curve.p[i];
        // auto curve_res = curve.y_2(curve.x[i] + dx_raw, i, j);
        // auto ret = dx_raw  / (curve.x[j] - curve_res);
        // return ret;
        return curve.p_2(i, j) * curve.p[j];
    }

    auto step_for_price_3(money p_min, money p_max, pair<int, int> p, money vol, money ext_vol) {
        // code size cleanup for n=2
        money ret = 0;
        return ret;
    }

    auto step_for_price_2(money p_min, money p_max, pair<int, int> p, money vol, money ext_vol) {
        money x0[2];
        copy_money_2(x0, &curve.x[0]);
        money _dx = 0;
        money _dy = 0;
        money x = 0;
        money y = 0;
        money price = 0;
        auto _from = p.first;
        auto _to = p.second;
        if (p_min > 0) {
            _from = p.second;
            _to = p.first;
        }
        auto step0 = dx / curve.p[_from];  // step in units of currency being sold
        auto step = step0;
        money gas = gas_fee / curve.p[_from];

        money previous_profit = 0;
        // print curve.p[_to] and curve.p[_from]
        // printf("to: %d, from: %d, curve.p[_to]: %Lf, curve.p[_from]: %Lf\n", _to, _from, curve.p[_to], curve.p[_from]);
        // + (step increases)
        // printf("p_min: %Lf, p_max: %Lf, vol: %Lf, ext_vol: %Lf\n", p_min, p_max, vol, ext_vol);
        while (true) {
            auto _dx_prev = _dx;
            auto _dy_prev = _dy;

            _dx += step;
            // printf("step: %Lf, _dx: %Lf, _dx_prev: %Lf\n", step, _dx, _dx_prev);
            // buy  -> x: first, y: second
            // sell -> x: second, y: first

            x = x0[_from] + _dx;
            y = curve.y_2(x, _from, _to);

            curve.x[_from] = x;
            curve.x[_to] = y;
            auto fee_mul = 1.L - this->fee_2();

            _dy = (x0[_to] - y) * fee_mul;
            curve.x[_to] = x0[_to] - _dy;

            if (_from == p.first) {
                price = _dx / _dy;
            }
            else {
                price = _dy / _dx;
            }
            auto v = vol + _dy * curve.p[_to];

            // Needed to prevent resonant trading which doesn't happen in reality
            copy_money_2(&curve.x[0], x0);  // restore the state
            // printf("::: %Lf %Lf %Lf %Lf\n", price, inst_price, p_min, p_max);
            // _from == p.first - buy
            // _from != p.first - sell
            money new_profit;
            if (_from == p.first)
                new_profit = (_dx / price - _dx / p_max) * p_max;
            else
                new_profit = (price - p_min) * _dx;
            // printf("new_profit: %Lf, previous_profit: %Lf\n", new_profit, previous_profit);
            if (new_profit > previous_profit and v <= ext_vol / 2.L) {
                previous_profit = new_profit;
            } else {
                _dx = _dx_prev;
                _dy = _dy_prev;
                break;
            }

            // step += step;
            step *= 2;
        }

        // - (step decreases)
        while (true) {
            auto _dx_prev = _dx;
            auto _dy_prev = _dy;
            if (step < 0) step = -step;
            // step /= 2;
            step *= 0.5;

            if (step < step0) {
                break;
            }

            for (int ctr=0;ctr<2;ctr++) {
                step = -step;
                _dx = _dx_prev + step;

                x = x0[_from] + _dx;
                y = curve.y_2(x, _from, _to);

                curve.x[_from] = x;
                curve.x[_to] = y;
                auto fee_mul = 1.L - this->fee_2();

                _dy = (x0[_to] - y) * fee_mul;
                curve.x[_to] = x0[_to] - _dy;

                if (_from == p.first) {
                    price = _dx / _dy;
                }
                else {
                    price = _dy / _dx;
                }
                auto v = vol + _dy * curve.p[_to];

                copy_money_2(&curve.x[0], x0);  // restore the state

                
                // _from == p.first - buy
                // _from != p.first - sell
                money new_profit;
                if (_from == p.first)
                    new_profit = (_dx / price - _dx / p_max) * p_max;
                else
                    new_profit = (price - p_min) * _dx;

                if (new_profit > previous_profit and v <= ext_vol / 2.L) {
                    previous_profit = new_profit;
                    break;
                } else {
                    _dx = _dx_prev;
                    _dy = _dy_prev;
                }
            }
        }
        // printf("*** p_min=%Lf, p_max=%Lf, _dy=%Lf, y=%Lf\n", p_min, p_max, _dy, curve.x[_to]);

        if (_from == p.first) {
            price = (_dx + gas) / _dy;  // need to buy higher than without gas
            previous_profit = (_dx / price - _dx / p_max) * p_max;
        }
        else {
            price = _dy / (_dx + gas); // need to sell lower than without gas
            previous_profit = (price - p_min) * _dx;
        }

        if (previous_profit <= 0) _dx = 0;
        //clamp dx for max_vol and min_vol

        if (_dx > 0) {
            _dx = min(max(step0, _dx), ext_vol / 2.L / curve.p[_from]);
        }
        return _dx;
    }
    
    void update_xcp_3(bool only_real=false) {
    }

    void update_xcp_2(bool only_real=false) {
        auto _xcp = get_xcp_2();
        auto old_xcp_profit_real = xcp_profit_real;
        // printf("t: %llu, xcp_profit_real: %10.12Lf, _xcp: %10.12Lf, xcp: %10.12Lf\n", t, xcp_profit_real, _xcp, xcp);

        xcp_profit_real = xcp_profit_real * _xcp / xcp;
        // printf("t: %llu, xcp_profit_real: %10.12Lf, \n", t, xcp_profit_real);
        if (not only_real) {
            xcp_profit += xcp_profit_real - old_xcp_profit_real;
        }
        xcp = _xcp;
    }

    inline void static copy_money_3(money *to, money const *from) {
        to[0] = from[0];
        to[1] = from[1];
        to[2] = from[2];
    }

    inline void static copy_money_2(money *to, money const *from) {
        to[0] = from[0];
        to[1] = from[1];
    }

    money exchange_3(money dx, int i, int j, money max_price=1e100L) {
        //"""
        //Buy y for x
        //"""
        try {
            money x_old[3];
            copy_money_3(x_old, &curve.x[0]);
            auto x = curve.x[i] + dx;
            auto y = curve.y_3(x, i, j);

            curve.x[i] = x;
            curve.x[j] = y;
            auto fee_mul = 1.L - this->fee_3();
            auto dy = x_old[j] - y;

            curve.x[j] = x_old[j] - dy * fee_mul;
            if ((dx / dy) > max_price or dy < 0) {
                copy_money_3(&curve.x[0], x_old);
                return 0;
            }
            update_xcp_3();
            return dy;
        } catch (...) {
            return 0;
        }
    }

    money exchange_2(money dx, int i, int j, money max_price=1e100L) {
        //"""
        //Buy y for x
        //"""
        try {
            money x_old[2];
            copy_money_2(x_old, &curve.x[0]);
            auto x = curve.x[i] + dx;
            auto y = curve.y_2(x, i, j);
            curve.x[i] = x;
            curve.x[j] = y;
            auto fee_mul = 1.L - this->fee_2();
            auto dy = x_old[j] - y;
            curve.x[j] = x_old[j] - dy * fee_mul;
            if ((dx / dy) > max_price or dy < 0) {
                copy_money_2(&curve.x[0], x_old);
                return 0;
            }
            update_xcp_2();
            return dy; //return is pre-fee
        } catch (...) {
            return 0;
        }
    }

    void ma_recorder(u64 t, vector<money> const &price_vector) {
        //  XXX what if every block only has p_b being last
        if (t > this->t) {
            money alpha = powl(0.5, ((money)(t - this->t) / this->ma_half_time));
            alpha = min(alpha, 1.L);
            // printf("t: %llu, dt: %llu, alpha: %Lf\n", t, t - this->t, alpha);
            // printf("price_vector: %Lf, price_oracle: %Lf\n", price_vector[1], price_oracle[1]);
            for (int k = 1; k < price_vector.size(); k++) {
                // printf("k: %d, price_vector[k]: %Lf, price_oracle[k]: %Lf\n", k, price_vector[k], price_oracle[k]);
                price_oracle[k] = price_vector[k] * (1 - alpha) + price_oracle[k] * alpha;
                // printf("k: %d, price_oracle[k]: %Lf\n", k, price_oracle[k]);
            }
            this->t = t;
        }
    }

    // auto tweak_price_2(u64 t, int a, int b, money p) {
    //     ma_recorder(t, last_price);
    //     const int N = 2;
    //     if (b > 0) {
    //         last_price[b] = p * last_price[a];
    //     } else {
    //         last_price[a] = last_price[0] / p;
    //     }

    auto tweak_price_2(u64 t, int /*a*/, int /*b*/, money spot_prev) {
        const int N = 2;
    
        // --- Feed the EMA with the pool's own spot (pre-fee marginal price),
        //     coin0 per coin1, computed at the current state.
        money amm_p01 = spot_prev;             // dx/dy (coin0 per coin1)
        // money amm_p01 = price_2(0, 1);
        // Optional: cap like the real pool (avoid extreme oracle jumps)
        money capped_p01 = std::min(amm_p01, 2.L * curve.p[1]);
    
        std::vector<money> spot = {1.L, capped_p01};
        ma_recorder(t, spot);
    

        // # price_oracle looks like [1, p1, p2, ...] normalized to 1e18
        money S = 0;
        for (size_t i = 0; i < N; i++) {
            auto t = price_oracle[i] / curve.p[i] - 1.L;
            S += t*t;
        }
        // printf("\nTime: %llu\n", t);
        // printf("xcp_profit: %10.12Lf, xcp_profit_real: %10.12Lf, allowed_extra_profit: %10.12Lf\n", xcp_profit, xcp_profit_real, allowed_extra_profit);
        // printf("2 * xcp_profit_real - 1.L: %10.12Lf, xcp_profit + 2 * allowed_extra_profit: %10.12Lf\n", 2 * xcp_profit_real - 1.L, xcp_profit + 2 * allowed_extra_profit);

        auto norm = S;
        norm = sqrt(norm); // .root_to();
        auto _adjustment_step = max(adjustment_step, norm / 5);
        if (norm <= _adjustment_step) {
            // Already close to the target price
            is_light = true;
            light_tx += 1;
            return norm;
        }
        // if (not not_adjusted and (xcp_profit_real > sqrt(xcp_profit) * (1.L + allowed_extra_profit))) {
        // THIS: 2 * vp - 1 > xcpp
        // POOL: vp - 1 > (xcpp - 1)/2 ::: 2vp - 2 > xcpp - 1 ::: 2vp-1 > xcpp
        
        if (not not_adjusted and (2 * xcp_profit_real - 1.L > xcp_profit + 2 * allowed_extra_profit)) {
            not_adjusted = true;
        }
        if (not not_adjusted) {
            light_tx += 1;
            is_light = true;
            return norm;
        }
        heavy_tx += 1;
        is_light = false;

        money p_new[MAX_ARRAY];
        p_new[0] = 1.L;
        for (size_t i = 1; i < price_oracle.size(); i++) {
            auto p_target = curve.p[i];
            auto p_real = price_oracle[i];
            p_new[i] = p_target + _adjustment_step * (p_real - p_target) / norm;
        }
        money old_p[MAX_ARRAY];
        copy_money_2(old_p, &curve.p[0]);

        auto old_profit = xcp_profit_real;
        auto old_xcp = xcp;

        copy_money_2(&curve.p[0],p_new);
        // printf("t: %llu, rebalance\n", t);
        // printf("old_profit: %10.12Lf, old_xcp: %10.12Lf\n", old_profit, old_xcp);
        update_xcp_2(true);
        // printf("t1: %llu, xcp_profit_real: %10.12Lf, xcp_profit: %10.12Lf\n", t, xcp_profit_real, xcp_profit);
        if (2 * xcp_profit_real - 1.0 <= xcp_profit ) {
            //  If real profit is less than half of maximum - revert params back
            copy_money_2(&curve.p[0], old_p);
            xcp_profit_real = old_profit;
            xcp = old_xcp;
            not_adjusted = false;
            // auto val = ((xcp_profit_real - 1.L - (xcp_profit - 1.L) / 2.L));
            // printf("%.10Lf\n", val);
        } else {
            // successful heavy rebalance
            this->n_rebalances += 1;
        }
        return norm;
    }


    auto tweak_price_3(u64 t, int a, int b, money p) {
        // code size cleanup for n=3
        money ret = 0;
        return ret;
    }


    void simulate(mapped_file const *in, simulation_data *simdata, extra_data *extdata) {
        // vector<trade_data> const &mdata
        const money CANDLE_VARIATIVES = 20;
        map<pair<int, int>, money> lasts;
        size_t N = price_oracle.size();
        u64 start_t = 0;
        u64 end_t = 0;
        long double last_time = 0;
        long double last_time_tweak_price = 0;
        size_t total_elements = in->size / sizeof(trade_data);
        printf("in->size: %zu\n", in->size);
        printf("sizeof(trade_data): %zu\n", sizeof(trade_data));
        printf("total_elements: %zu\n", total_elements);
        simdata->total = total_elements;
        auto mapped_data = (trade_data const *) in->base;
        auto mapped_data_ptr = mapped_data;
        money slippage = 0;
        money antislippage = 0;
        money slippage_count = 0;
        money _slippage = 0; // initialize to avoid using garbage when price doesn't move
        money spot_prev = price_2(0, 1);
        money last_prices = price_2(0, 1);
        money previous_price_scale = curve.p[1];
        // Track TVL growth in coin0 units and HODL baseline
        // TVL in coin0 units: sum_i x[i] * p[i] (p[0] == 1)
        auto tvl_in_coin0 = [&](vector<money> const &x, vector<money> const &p) -> long double {
            long double v = 0;
            for (size_t i = 0; i < x.size(); i++) v += (long double)(x[i] * p[i]);
            return v;
        };
        vector<money> x_start = curve.x; // initial LP balances by coin
        long double tvl_start = tvl_in_coin0(curve.x, curve.p);
        long double donation_coin0_total = 0.0L;

        FILE *out_file;
        FILE *trades_file = nullptr;
        if (log) {
            out_file = fopen("states.json", "w");
            fprintf(out_file, "[");
            // Open per-configuration trades log (JSON Lines)
            char trades_name[64];
            snprintf(trades_name, sizeof(trades_name), "trades-%d.jsonl", simdata->num);
            trades_file = fopen(trades_name, "w");
        }

        long double arb_profit_sum_coin0 = 0.0L;
        // Integral accumulator for |p_cex - price_scale| over time (seconds)

        long double total_cex_diff = 0.0L;
        long double cex_diff_squared = 0.0L;
        long double max_cex_diff = 0.0L;
        // Accumulator: sum of dt where relative deviation exceeds threshold
        unsigned long long prev_t_for_dev = 0ULL;
        long double sum_dt_dev_exceeds = 0.0L;

        for (size_t i = 0; i < total_elements; i++) {
            simdata->current = i;
            // if (i > 10) abort();
            trade_data d = *mapped_data_ptr++;
            if (i == 0) start_t = d.t;
            if (i == 0) last_time_tweak_price = d.t;
            if (i == 0) this->t = d.t;
            end_t = d.t;
            if (last_time > 0) {
                last_time = d.t - last_time;
            }
            // Track time in seconds between events for deviation fraction metric
            unsigned long long dt_for_dev = 0ULL;
            if (prev_t_for_dev != 0ULL) dt_for_dev = d.t - prev_t_for_dev;
            prev_t_for_dev = d.t;
            auto p_diff = fabsl(d.close - curve.p[1]);
            total_cex_diff += p_diff;
            cex_diff_squared += p_diff * p_diff;
            if (p_diff > max_cex_diff) max_cex_diff = p_diff;
            // If relative deviation exceeds threshold, add this dt
            if (d.close > 0) {
                long double rel_dev = fabsl(((long double)curve.p[1] - (long double)d.close) / (long double)d.close);
                if (rel_dev > (long double)this->cex_diff_rel_thresh) {
                    sum_dt_dev_exceeds += (long double)dt_for_dev;
                }
            }
            auto a = d.pair1.first;
            auto b = d.pair1.second;
            money vol{0.L};
            auto ext_vol = money(d.volume / price_oracle[a]); //  <- now all is in USD
            // printf("ext_vol: %Lf, d.volume: %Lf, price_oracle[b]: %Lf\n", static_cast<long double>(ext_vol), static_cast<long double>(d.volume), static_cast<long double>(price_oracle[b]));

            int ctr{0};
            money last;
            auto itl = lasts.find({a, b});
            if (itl == lasts.end()) {
                last = price_oracle[b] / price_oracle[a];
            } else {
                last = itl->second;
            }
            auto _high = last;
            auto _low = last;
            _slippage = 0; // reset per-iteration before any accumulation

            auto max_price = d.high * (1 - ext_fee);
            auto min_price = d.low * (1 + ext_fee);
            money _dx = 0;
            auto p_before = N == 3 ? price_3(a, b) : price_2(a, b);
            bool trade_happened = false;

            if ((max_price != 0) & (max_price > p_before)) {
                auto step = N == 3 ? step_for_price_3(0, max_price, d.pair1, vol, ext_vol) : step_for_price_2(0, max_price, d.pair1, vol, ext_vol);
                if (step > 0) {
                    // printf("+++ %Lf %Lf %d %d\n", curve.x[a], curve.x[b], a, b);
                    // Capture pre-exchange balances to compute taker amounts
                    auto bal_a_before = curve.x[a];
                    auto bal_b_before = curve.x[b];
                    long double vp_pre_log = xcp_profit_real;
                    long double xcp_profit_pre_log = xcp_profit;

                    auto dy = N == 3 ? exchange_3(step, a, b) : exchange_2(step, a, b);
                    trade_happened = true;
                    // printf("+++ %Lf %Lf\n", curve.x[a], curve.x[b]);
                    vol += step * price_oracle[a];
                    _dx += dy;
                    last = N == 3 ? price_3(a, b) : price_2(a, b);
                    ctr += 1;

                    // Approximate arbitrage profit in coin0 units when buying coin b on pool and selling b on CEX
                    if (price_oracle.size() == 2) {
                        long double dx_in = (long double)(curve.x[a] - bal_a_before); // in coin a
                        long double dy_after_fee = (long double)(bal_b_before - curve.x[b]); // in coin b
                        if (dx_in > 0 && dy_after_fee > 0) {
                            long double cex_price = (long double)d.high; // coin0 per coin1
                            long double f_buy = 1.0L - (long double)this->ext_fee;
                            long double gas_coin0 = (long double)this->gas_fee;
                            long double profit_coin0;
                            if (a == 0) {
                                // Buy coin1 on pool with coin0; sell coin1 on CEX
                                profit_coin0 = dy_after_fee * cex_price * f_buy - dx_in - gas_coin0;
                            } else {
                                // Buy coin0 on pool with coin1; sell coin0 on CEX
                                // Revenue is in coin0 already; apply CEX sell fee only
                                profit_coin0 = dy_after_fee * f_buy - dx_in * cex_price - gas_coin0;
                            }
                            if (profit_coin0 > 0) arb_profit_sum_coin0 += profit_coin0;
                        }
                    }

                    // Debug trade log: upper branch (sell coin a, buy coin b)
                    if (trades_file && price_oracle.size() == 2) {
                        long double dx_in_log = (long double)(curve.x[a] - bal_a_before);
                        long double dy_out_log = (long double)(bal_b_before - curve.x[b]);
                        long double cex_price_log = (long double)d.high;
                        long double pool_price_dbg = (N == 3 ? price_3(a, b) : price_2(a, b));
                        long double vp_post_log = xcp_profit_real;
                        long double xcp_profit_post_log = xcp_profit;
                        long double gas_coin0_log = (long double)this->gas_fee;
                        long double f_buy_log = 1.0L - (long double)this->ext_fee;
                        long double profit_coin0_log = 0.0L;
                        if (a == 0) {
                            profit_coin0_log = dy_out_log * cex_price_log * f_buy_log - dx_in_log - gas_coin0_log;
                        } else {
                            profit_coin0_log = dy_out_log * f_buy_log - dx_in_log * cex_price_log - gas_coin0_log;
                        }
                        fprintf(trades_file,
                                "{\"t\": %llu, \"branch\": \"upper\", \"from\": %d, \"to\": %d, \"dx\": %.12Lf, \"dy\": %.12Lf, \"cex_price\": %.12Lf, \"pool_price_before\": %.12Lf, \"pool_price_after\": %.12Lf, \"profit_coin0\": %.12Lf, \"x0_before\": %.12Lf, \"x1_before\": %.12Lf, \"x0_after\": %.12Lf, \"x1_after\": %.12Lf, \"vp_pre\": %.12Lf, \"vp_post\": %.12Lf, \"xcp_profit_pre\": %.12Lf, \"xcp_profit_post\": %.12Lf}\n",
                                d.t, a, b, dx_in_log, dy_out_log, cex_price_log, last_prices, pool_price_dbg, profit_coin0_log,
                                bal_a_before, bal_b_before, curve.x[0], curve.x[1], vp_pre_log, vp_post_log, xcp_profit_pre_log, xcp_profit_post_log);
                    }
                }
            }

            auto p_after = N == 3 ? price_3(a, b) : price_2(a, b);
            // auto _fee = N == 3 ? fee_3() : fee_2();
            // printf("!!! %d %d\n", a, b);
            // printf("!!!1 %Lf %Lf\n", p_after, last);

            if (p_before != p_after) {
                auto v = _dx / curve.x[b];
                _slippage = (_dx * (p_before + p_after)) / (2.L * (mabs(p_before - p_after)) * curve.x[b]);
                volume += v;
            }
            if (_slippage > 0) {
                slippage_count += last_time;
                antislippage += last_time * _slippage;
                slippage += last_time / _slippage;
            }

            _high = last;
            _dx = 0;
            p_before = p_after;
            
            if ((min_price != 0) & (min_price < p_before)) {
                auto step = N == 3 ? step_for_price_3(min_price, 0, d.pair1, vol, ext_vol) : step_for_price_2(min_price, 0, d.pair1, vol, ext_vol);
                if (step > 0) {
                    // printf("=== %Lf %Lf %d %d\n", curve.x[a], curve.x[b], a, b);
                    // Capture pre-exchange balances to compute taker amounts
                    auto bal_a_before = curve.x[a];
                    auto bal_b_before = curve.x[b];
                    auto log_vp_pre = xcp_profit_real;
                    auto log_xcp_profit_pre = xcp_profit;
                    auto dy = N == 3 ? exchange_3(step, b, a) : exchange_2(step, b, a);
                    trade_happened = true;
                    auto log_vp_post = xcp_profit_real;
                    auto log_xcp_profit_post = xcp_profit;

                    // printf("=== %Lf %Lf %d %d\n", curve.x[a], curve.x[b], a, b);
                    // printf("!===! %Lf %Lf\n", step, dy);
                    vol += dy * price_oracle[a];
                    _dx += step;
                    last = N == 3 ? price_3(a, b) : price_2(a, b);
                    ctr += 1;
                    

                    // Approximate arbitrage profit in coin0 units when buying coin a (possibly coin0) on pool and selling a on CEX
                    if (price_oracle.size() == 2) {
                        // In this path we call exchange(..., from=b, to=a)
                        long double dx_in = (long double)(curve.x[b] - bal_b_before); // spent coin b
                        long double dy_after_fee = (long double)(bal_a_before - curve.x[a]); // received coin a
                        if (dx_in > 0 && dy_after_fee > 0) {
                            long double cex_price = (long double)d.low; // coin0 per coin1
                            long double f_sell = 1.0L - (long double)this->ext_fee;
                            long double gas_coin0 = (long double)this->gas_fee;
                            long double profit_coin0;
                            if (a == 0) {
                                // Buy coin0 on pool with coin1; sell coin0 on CEX
                                profit_coin0 = dy_after_fee * f_sell - dx_in * cex_price - gas_coin0;
                            } else {
                                // Buy coin1 on pool with coin0; sell coin1 on CEX
                                profit_coin0 = dy_after_fee * cex_price * f_sell - dx_in - gas_coin0;
                            }
                            if (profit_coin0 > 0) arb_profit_sum_coin0 += profit_coin0;
                        }
                    }

                    // Debug trade log: lower branch (sell coin b, buy coin a)
                    if (trades_file && price_oracle.size() == 2) {
                        long double dx_in_log = (long double)(curve.x[b] - bal_b_before);
                        long double dy_out_log = (long double)(bal_a_before - curve.x[a]);
                        long double cex_price_log = (long double)d.low;
                        long double pool_price_dbg = (N == 3 ? price_3(a, b) : price_2(a, b));
                        long double gas_coin0_log = (long double)this->gas_fee;
                        long double f_sell_log = 1.0L - (long double)this->ext_fee;
                        long double profit_coin0_log = 0.0L;
                        if (a == 0) {
                            profit_coin0_log = dy_out_log * f_sell_log - dx_in_log * cex_price_log - gas_coin0_log;
                        } else {
                            profit_coin0_log = dy_out_log * cex_price_log * f_sell_log - dx_in_log - gas_coin0_log;
                        }
                        fprintf(trades_file,
                                "{\"t\": %llu, \"branch\": \"lower\", \"from\": %d, \"to\": %d, \"dx\": %.12Lf, \"dy\": %.12Lf, \"cex_price\": %.12Lf, \"pool_spot_before\": %.12Lf, \"pool_spot_after\": %.12Lf, \"profit_coin0\": %.12Lf, \"x0_before\": %.12Lf, \"x1_before\": %.12Lf, \"x0_after\": %.12Lf, \"x1_after\": %.12Lf, \"vp_pre\": %.12Lf, \"vp_post\": %.12Lf, \"xcp_profit_pre\": %.12Lf, \"xcp_profit_post\": %.12Lf}\n",
                                d.t, b, a, dx_in_log, dy_out_log, cex_price_log, last_prices, pool_price_dbg, profit_coin0_log,
                                bal_a_before, bal_b_before, curve.x[0], curve.x[1], log_vp_pre, log_vp_post, log_xcp_profit_pre, log_xcp_profit_post);
                    }
                }
            }

            p_after = N == 3 ? price_3(a, b) : price_2(a, b);
            // _fee = N == 3 ? fee_3() : fee_2();
            // printf("!!!2 %Lf %Lf\n", p_after, last);

            if (p_before != p_after) {
                auto v = _dx / curve.x[b];
                _slippage = (_dx * (p_before + p_after)) / (2.L * (mabs(p_before - p_after)) * curve.x[b]);
                volume += v;
            }
            if (_slippage > 0) {
                slippage_count += last_time;
                antislippage += last_time * _slippage;
                slippage += last_time / _slippage;
            }

            _low = last;
            lasts[d.pair1] = last;
            

            // Boost with special donations to the pool
            if (this->boost_rate > 0) {
                auto _boost = (1.L + last_time * this->boost_rate);
                // donation amount in coin0 units before applying boost
                long double tvl_before = tvl_in_coin0(curve.x, curve.p);
                if (_boost > 1.L) {
                    // printf("boost: %10.12Lf, boost_rate: %10.12Lf\n", _boost, this->boost_rate);
                    // Donation equals the increase of each balance converted to coin0 units
                    // i.e., sum_i x[i] * (boost-1) * p[i] == tvl_before * (boost-1)
                    long double donation0 = (long double)curve.x[0] * ((long double)_boost - 1.0L);
                    long double donation1 = (long double)curve.x[1] * ((long double)_boost - 1.0L);
                    long double donation_coin0 = donation0 + donation1 * (curve.p.size() > 1 ? (long double)curve.p[1] : 1.0L);
                    donation_coin0_total += donation_coin0;
                    // if (trades_file) {
                    //     fprintf(trades_file,
                    //             "{\"t\": %llu, \"type\": \"donation\", \"boost\": %.12Lf, \"amt0\": %.12Lf, \"amt1\": %.12Lf, \"coin0_value\": %.12Lf}\n",
                    //             d.t, _boost, donation0, donation1, donation_coin0);
                    // }
                }
                curve.x[0] = curve.x[0] * _boost;
                curve.x[1] = curve.x[1] * _boost;
                if (N == 3) {
                    curve.x[2] = curve.x[2] * _boost;
                }
                xcp_profit_real *= _boost;
                xcp *= _boost;
                this->boost_integral *= _boost;
            }

            long double norm = 0;
            // only tweak_price every N seconds or on trade
            if (d.t - last_time_tweak_price >= 600 || trade_happened) {
                previous_price_scale = curve.p[1];
                money log_oracle_pre = price_oracle[1];
                money cur_get_p = curve.p_2(0, 1);
                money log_vp_pre = xcp_profit_real;
                money log_xcp_profit_pre = xcp_profit;
                if (N == 2) norm = tweak_price_2(d.t, a, b, last_prices);
                else        norm = tweak_price_3(d.t, a, b, (_high + _low) / 2.L);
                // spot_prev = price_2(0, 1) * ps_pre / curve.p[1]; 
                last_prices = cur_get_p * previous_price_scale;
                money log_vp_post = xcp_profit_real;
                money log_xcp_profit_post = xcp_profit;
                money log_ps_post = curve.p[1];
                money log_oracle_post = price_oracle[1];
                last_time_tweak_price = d.t;
                if (trades_file) {
                    fprintf(trades_file,
                            "{\"t\": %llu, \"type\": \"tweak_price\", \"trade_happened\": %d, \"p_cex\": %.12Lf, \"ps_pre\": %.12Lf, \"ps_post\": %.12Lf, \"oracle_pre\": %.12Lf, \"oracle_post\": %.12Lf, \"vp_pre\": %.12Lf, \"vp_post\": %.12Lf, \"xcp_profit_pre\": %.12Lf, \"xcp_profit_post\": %.12Lf}\n",
                            d.t, trade_happened, d.high, previous_price_scale, log_ps_post, log_oracle_pre, log_oracle_post, log_vp_pre, log_vp_post, log_xcp_profit_pre, log_xcp_profit_post);
                }

            }

            total_vol += vol;
            this->n_trades += (unsigned long long)ctr;
            last_time = d.t;
            long double ARU_x = xcp_profit_real;
            long double ARU_y = (86400.L * 365.L / (d.t - start_t + 1.L));
            APY = powl(ARU_x, ARU_y) - 1.L;
            APY_boost = powl(xcp_profit_real / this->boost_integral, ARU_y) - 1.L;
            if (i % 1024 == 0 && log) {
                try {
                    long double last01, last02 = 0.0;
                    auto it01 = lasts.find({0, 1});
                    if (it01 == lasts.end()) {
                        last01 = price_oracle[1] / price_oracle[0];
                    } else {
                        last01 = it01->second;
                    }
                    if (N == 3) {
                        auto it02 = lasts.find({0, 2});
                        if (it02 == lasts.end()) {
                            last02 = price_oracle[2] / price_oracle[0];
                        } else {
                            last02 = it02->second;
                        }
                    }
                    if (N == 3) {
                        printf("t=%llu %.1Lf%%\ttrades: %d\t"
                               "AMM: %.3Lf, %0.3Lf\tTarget: %.3Lf, %.3Lf\t"
                               "Vol: %.4Lf\tPR:%.2Lf\txCP-growth: {%.10Lf}\t"
                               "APY:%.1Lf%%\tfee:%.3Lf%% %c\n",
                               d.t,
                               100.L * i / total_elements, ctr, last01, last02,
                               curve.p[1],
                               curve.p[2],
                               total_vol,
                               (xcp_profit_real - 1.) / (xcp_profit - 1.L),
                               xcp_profit_real,
                               APY * 100.L,
                               (curve.p.size() == 3 ? fee_3() : fee_2()) * 100.L,
                               is_light ? '*' : '.');
                    } else if (N == 2) {
                        printf("t=%llu %.1Lf%%\ttrades: %d\tAMM: %.5Lf\tTarget: %.5Lf\tVol: %.4Lf\tPR:%.2Lf\txCP-growth: {%.10Lf}\tAPY:%.1Lf%%\tfee:%.3Lf%% %c\n",
                                d.t,
                                100.L * i / total_elements,
                                ctr,
                                last01,
                                curve.p[1],
                                total_vol,
                                (xcp_profit_real - 1.) / (xcp_profit - 1.L),
                                xcp_profit_real,
                                APY * 100.L,
                                fee_2() * 100.L,
                                is_light ? '*' : '.');

                    }
                } catch (std::exception const &e) {
                    printf("caught '%s'\n", e.what());
                }
            }

            if (log) {
                fprintf(out_file, "{\"t\": %llu, \"token0\": %.6Lf, \"token1\": %.6Lf, \"price_oracle\": %.6Lf, \"price_scale\": %.6Lf, \"profit\": %.6Lf, \"open\": %.6Lf, \"high\": %.6Lf, \"low\": %.6Lf, \"close\": %.6Lf, \"snapshots\": ",
                        d.t,
                        curve.x[0],
                        curve.x[1],
                        price_oracle[b] / price_oracle[a],
                        curve.p[1],
                        xcp_profit_real - 1.0,
                        d.open, d.high, d.low, d.close);
                fputc('[', out_file);
                fputc(']', out_file);
                fputc('}', out_file);
                if (i < total_elements - 1) {
                    fprintf(out_file, ",\n");
                }
            }

            if (slippage > 1e20) {
                printf("*** Slippage is too high %.5Lf\n", slippage);
            }
        }
        // Compute APY metrics based on coin0 growth vs HODL baseline
        long double duration_s = (end_t > start_t) ? ((long double)end_t - (long double)start_t) : 0.0L;
        long double exponent = (duration_s > 0.0L) ? ((86400.L * 365.L) / duration_s) : 0.0L;
        long double tvl_end = tvl_in_coin0(curve.x, curve.p);
        long double v_hold_end = 0.0L; // value of initial LP holdings if simply held, priced at end
        for (size_t i = 0; i < x_start.size(); i++) v_hold_end += (long double)(x_start[i] * curve.p[i]);
        long double avg_cex_diff = total_cex_diff / duration_s;
        long double cex_follow_time_frac = (duration_s > 0.0L) ? (1 - sum_dt_dev_exceeds / duration_s) : 0.0L;
        // Virtual price APY is what we already compute as APY via xcp_profit_real
        long double apy_vp = xcp_profit_real - 1.0L;

        long double apy_coin0_raw = -1.0L;
        long double apy_coin0_boost_raw = -1.0L;
        long double apy_coin0 = -1.0L;
        long double apy_coin0_boost = -1.0L;

        if (duration_s > 0.0L && tvl_start > 0.0L) {
            if (tvl_end > 0.0L) {
                apy_coin0_raw = powl((long double)(tvl_end / tvl_start), exponent) - 1.0L;
            }
            long double tvl_end_adj_raw = tvl_end - donation_coin0_total;
            if (tvl_end_adj_raw > 0.0L) {
                apy_coin0_boost_raw = powl((long double)(tvl_end_adj_raw / tvl_start), exponent) - 1.0L;
            }
            if (v_hold_end > 0.0L) {
                apy_coin0 = powl((long double)(tvl_end / v_hold_end), exponent) - 1.0L;
                long double tvl_end_adj = tvl_end - donation_coin0_total;
                if (tvl_end_adj > 0.0L) {
                    apy_coin0_boost = powl((long double)(tvl_end_adj / v_hold_end), exponent) - 1.0L;
                }
            }
        }

        extdata->slippage = slippage / slippage_count / 2.L;
        extdata->liq_density = 2.L * antislippage / slippage_count;
        extdata->APY = APY;
        extdata->volume = volume;
        extdata->APY_boost = APY_boost;
        extdata->vpminusone = apy_vp;
        extdata->APY_coin0_raw = apy_coin0_raw;
        extdata->APY_coin0_boost_raw = apy_coin0_boost_raw;
        extdata->APY_coin0 = apy_coin0;
        extdata->APY_coin0_boost = apy_coin0_boost;
        extdata->donation_coin0_total = (money)donation_coin0_total;
        extdata->trade_volume = total_vol;
        extdata->arb_profit_coin0 = (money)arb_profit_sum_coin0;
        extdata->n_rebalances = this->n_rebalances;
        extdata->n_trades = this->n_trades;
        extdata->xcp_profit_real = xcp_profit_real;
        extdata->xcp_profit = xcp_profit;
        extdata->avg_cex_diff = (money)avg_cex_diff;
        extdata->max_cex_diff = (money)max_cex_diff;
        extdata->cex_diff_l2 = sqrt((money)cex_diff_squared);
        extdata->cex_follow_time_frac = (money)cex_follow_time_frac;
        if (log) {
            fprintf(out_file, "]");
            if (trades_file) fclose(trades_file);
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
    money adjustment_step;
    money allowed_extra_profit;
    int log;
    money fee_gamma;
    money total_vol;
    int ma_half_time;
    money ext_fee;
    money gas_fee;
    money boost_rate;
    money boost_integral;
    money volume;
    money slippage;
    money antislippage;
    money slippage_count;
    long double APY;
    long double APY_boost;
    bool not_adjusted;
    int  heavy_tx;
    int  light_tx;
    bool is_light;
    unsigned long long n_rebalances;
    unsigned long long n_trades;
    money cex_diff_rel_thresh;
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


bool simulation(simulation_data *data) {
//        int num, json const *jconf, vector<money> const *price_vector, mapped_file const *test_data) {
    Trader trader(*(data->jconf), *(data->price_vector));
    auto start_simulation = get_thread_time();
    printf("Configuration %d: begin simulation\n", data->num);
    unlink(data->test_data->name.c_str()); // Temp file can be deleted in *nix even being open
    extra_data extdata;
    trader.simulate(data->test_data, data, &extdata);
    data->result = extdata;
    //money liq_density = jout["liq_density"];
    //money APY = jout["APY"];
    printf("Liquidity density vs that of xyz=k: %Lf\n", extdata.liq_density);
    printf("APY-boost: %Lf%%\n", extdata.APY_boost * 100.L);
    printf("APY: %Lf%%\n", extdata.APY * 100.L);
//    json_save(out_json_name, jout);
    auto end = get_thread_time();
    print_clock("Total simulation time", start_simulation, end);
    return true;
}

struct work_queue {
    std::queue<simulation_data> *wq;
    pthread_mutex_t *lock;
    pthread_mutex_t *result_lock;
    json *result;
    int num;
};

void *simulation_thread(void *args) {
    auto data = (work_queue *)args;
    printf("[%d]: simulation thread started\n", data->num);
    for (;;) {
        pthread_mutex_lock(data->lock);
        if (data->wq->empty()) {
            pthread_mutex_unlock(data->lock);
            break;
        }
        simulation_data simdata = data->wq->front();
        data->wq->pop();
        pthread_mutex_unlock(data->lock);
        printf("[%d]: pick up configuration %d\n", data->num, simdata.num);
        simulation(&simdata);
        pthread_mutex_lock(data->result_lock);
        (*(data->result))["configuration"][simdata.num]["Result"]["APY"] = simdata.result.APY;
        (*(data->result))["configuration"][simdata.num]["Result"]["vpminusone"] = simdata.result.vpminusone;
        (*(data->result))["configuration"][simdata.num]["Result"]["APY_coin0_raw"] = simdata.result.APY_coin0_raw;
        (*(data->result))["configuration"][simdata.num]["Result"]["APY_coin0_boost_raw"] = simdata.result.APY_coin0_boost_raw;
        (*(data->result))["configuration"][simdata.num]["Result"]["APY_coin0"] = simdata.result.APY_coin0;
        (*(data->result))["configuration"][simdata.num]["Result"]["APY_coin0_boost"] = simdata.result.APY_coin0_boost;
        (*(data->result))["configuration"][simdata.num]["Result"]["liq_density"] = simdata.result.liq_density;
        (*(data->result))["configuration"][simdata.num]["Result"]["slippage"] = simdata.result.slippage;
        (*(data->result))["configuration"][simdata.num]["Result"]["volume"] = simdata.result.volume;
        (*(data->result))["configuration"][simdata.num]["Result"]["APY_boost"] = simdata.result.APY_boost;
        (*(data->result))["configuration"][simdata.num]["Result"]["trade_volume"] = simdata.result.trade_volume;
        (*(data->result))["configuration"][simdata.num]["Result"]["n_rebalances"] = simdata.result.n_rebalances;
        (*(data->result))["configuration"][simdata.num]["Result"]["n_trades"] = simdata.result.n_trades;
        (*(data->result))["configuration"][simdata.num]["Result"]["donation_coin0_total"] = simdata.result.donation_coin0_total;
        (*(data->result))["configuration"][simdata.num]["Result"]["arb_profit_coin0"] = simdata.result.arb_profit_coin0;
        (*(data->result))["configuration"][simdata.num]["Result"]["xcp_profit_real"] = simdata.result.xcp_profit_real;
        (*(data->result))["configuration"][simdata.num]["Result"]["avg_cex_diff"] = simdata.result.avg_cex_diff;
        (*(data->result))["configuration"][simdata.num]["Result"]["max_cex_diff"] = simdata.result.max_cex_diff;
        (*(data->result))["configuration"][simdata.num]["Result"]["cex_diff_l2"] = simdata.result.cex_diff_l2;
        (*(data->result))["configuration"][simdata.num]["Result"]["cex_follow_time_frac"] = simdata.result.cex_follow_time_frac;
        (*(data->result))["configuration"][simdata.num]["Result"]["xcp_profit"] = simdata.result.xcp_profit;
        pthread_mutex_unlock(data->result_lock);

    }
    printf("[%d]: simulation thread ended\n", data->num);
    return nullptr;
}

int main(int argc, char **argv) {
    if (argc == 1) {
        printf("Usage: %s [trim] [threads=#] [in-json-file] [out-json-file]\n", argv[0]);
        return 0;
    }
    int LAST_ELEMS = 0;
    if (argc > 1 && std::string(argv[1]).find("trim") != std::string::npos) {
        if (argv[1][4] == 0) LAST_ELEMS = 1000000;
        else                 LAST_ELEMS = atoi(argv[1]+4);
        argc--; argv++;
    }
    int THREADS = 1;
    if (argc > 1 && std::string(argv[1]).find("threads=") != std::string::npos) {
        THREADS = atoi(argv[1]+8);
        argc--; argv++;
    }
    string in_json_name = argc > 1 ? argv[1] : "sample_in.json";
    string out_json_name = argc > 2 ? argv[2] : "sample_out.json";
    double real_time_start = get_total_time();
    json jin;
    if (!json_load(in_json_name, jin)) {
        return 0;
    }
    int configurations = jin["configuration"].size();
    if (configurations <= 0) {
        printf("No configurations found\n");
        return 0;
    }

    printf("Total %d configurations will be processed in %d threads\n", configurations, THREADS);
    //for (auto const &cfg: jin["configuration"]) {
    //    std::cout << std::setw(4) << cfg << "\n";
    // }
    //
    vector<money> price_vector;
    mapped_file test_data;
    if (!get_all(jin, LAST_ELEMS, price_vector, test_data)) {
        return 0;
    }
    //debug_print("test_data first 5", test_data, 5);
    //debug_print("test_data last 5", test_data, -5);
    double time_start = get_total_time();
    double wall_time_start = get_wall_time();
    std::queue<simulation_data> sim_queue;
    vector<pthread_t> threads(THREADS);
    vector<work_queue> thr_data(THREADS);
    pthread_mutex_t queue_mutex;
    pthread_mutex_init(&queue_mutex, nullptr);
    pthread_mutex_t result_mutex;
    pthread_mutex_init(&result_mutex, nullptr);
    json result = jin;
    for (int i = 0; i < configurations; i++) {
        simulation_data cd;
        cd.num = i;
        cd.test_data = &test_data;
        cd.price_vector = &price_vector;
        cd.jconf = &jin["configuration"][i];
        cd.current = 0;
        cd.total = 0;
        sim_queue.push(cd);
    }
    for (int i = 0; i < THREADS; i++) {
        thr_data[i].wq = &sim_queue;
        thr_data[i].lock = &queue_mutex;
        thr_data[i].result_lock = &result_mutex;
        thr_data[i].num = i;
        thr_data[i].result = &result;
        if (pthread_create(&threads[i], nullptr, simulation_thread, &thr_data[i]) != 0) {
            printf("Can't create thread %d!\n", i);
        }
    }

    for (int i = 0; i < THREADS; i++) {
        pthread_join(threads[i], nullptr);
    }
    json_save(out_json_name, result);
    double time_end = get_total_time();
    double wall_time_end = get_wall_time();
    print_clock("Data reading and preprocessing time", real_time_start, time_start);
    print_clock("Total simulation wall time", wall_time_start, wall_time_end);
    print_clock("Total simulation processor time", time_start, time_end);
    printf("Parallelizm ratio %.5lf\n", (time_end - time_start) / (wall_time_end - wall_time_start));
}
