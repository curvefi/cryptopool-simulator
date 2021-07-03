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
    auto start_time = get_thread_time();
    //debug_print("out first 5", out, 5);
    //debug_print("out last 5", out, -5);
    sort(out.begin(), out.end(), [](trade_one const &l, trade_one const &r) {
        return l.t < r.t;
    });
    auto end_time = get_thread_time();
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
    price_vector = get_price_vector(N, ret);
    string tmp_name = "_tmp." + std::to_string(getpid());
    FILE *tmp = fopen(tmp_name.c_str(), "w+");
    if (tmp == nullptr) {
        printf("Temp file '%s' is not available\n", tmp_name.c_str());
        return false;
    }
    printf("Using temp file '%s' as interprocedural connect\n", tmp_name.c_str());
    fwrite(&ret[0], sizeof(trade_data), ret.size(), tmp);
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
        if (diff <= 1E-18 or diff * 1E18L < D) {
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
        K = gamma / (gamma + 1.L - K);
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

auto newton_D_2(money A, money gamma, money const *xx, money D0) {
    money D = D0;
    money S = 0;
    const size_t N = 2;
    money x[2];
    for (size_t i = 0; i < N; i++) {
        S += x[i] = xx[i];
    }
    sort(x+0, x+N, [](money l, money r) { return l > r; });
    auto NN = 1.L;
    for (size_t j = 0; j < N; j++) { // XXX or just set A to be A*N**N?
        NN *= N;
    }
    A *= NN;
    money rev_gamma = 1.L / gamma;
    money gamma_1 = 1.L + gamma;
    //for (size_t j = 0; j < N; j++) { // XXX or just set A to be A*N**N?
    //    A = A * N;
    //}

    for (int i = 0; i < 255; i++) {
        money D_prev = D;

        money K0 = NN;
        for (auto const &_x: x) {
            K0 = K0 * _x / D;
        }

        money _g1k0 = mabs((gamma_1 - K0));

        // # D / (A * N**N) * _g1k0**2 / gamma**2
        money mul1 = D * rev_gamma * _g1k0 * rev_gamma * _g1k0 / A;

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

    money A;
    money gamma;
    size_t n;
    vector<money> p;
    vector<money> x;

};

struct extra_data {
    money APY = 0;
    money liq_density = 0;
    money slippage = 0;
    money volume = 0;
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
        price_threshold = jconf["price_threshold"];
        fee_gamma = jconf["fee_gamma"];
        adjustment_step = jconf["adjustment_step"];
        ma_half_time = jconf["ma_half_time"];
        this->ext_fee = jconf["ext_fee"];
        log = jconf["log"];
        this->p0 = p0;
        this->price_oracle = this->p0;
        this->last_price = this->p0;
        // this->curve = Curve(A, gamma, D, n, p0);
        this->dx = D * 1e-8L;
        this->D0 = n == 3 ? this->curve.D_3() : this->curve.D_2();
        this->xcp_0 = n == 3 ? this->get_xcp_3() : this->get_xcp_2();
        this->xcp_profit = 1.L;
        this->xcp_profit_real = 1.L;
        this->xcp = this->xcp_0;
        this->total_vol = 0.0;
        this->slippage = 0;
        this->antislippage = 0;
        this->slippage_count = 0;
        this->volume = 0;
        this->not_adjusted = false;
        this->heavy_tx = 0;
        this->light_tx = 0;
        this->is_light = false;
        this->t = 0;
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
        auto dx_raw = dx  / curve.p[i];
        auto curve_res = curve.y_2(curve.x[i] + dx_raw, i, j);
        auto ret = dx_raw  / (curve.x[j] - curve_res);
        return ret;
    }

    auto step_for_price_3(money dp, pair<int, int> p, int sign) {
        auto p0 = price_3(p.first, p.second);
        dp = p0 * dp;
        money x0[3];
        copy_money_3(x0, &curve.x[0]);
        auto step = dx / curve.p[p.first];
        while (true) {
            curve.x[p.first] = x0[p.first] + sign * step;
            auto dp_ = mabs(p0 - price_3(p.first, p.second));
            if (dp_ >= dp or step >= curve.x[p.first] / 10.L) {
                copy_money_3(&curve.x[0], x0);
                return step;
            }
            step += step;
        }
    }

    auto step_for_price_2(money dp, pair<int, int> p, int sign) {
        auto p0 = price_2(p.first, p.second);
        dp = p0 * dp;
        money x0[2];
        copy_money_2(x0, &curve.x[0]);
        auto step = dx / curve.p[p.first];
        while (true) {
            curve.x[p.first] = x0[p.first] + sign * step;
            auto dp_ = mabs(p0 - price_2(p.first, p.second));
            if (dp_ >= dp or step >= curve.x[p.first] / 10.L) {
                copy_money_2(&curve.x[0], x0);
                return step;
            }
            step += step;
        }
    }

    void update_xcp_3(bool only_real=false) {
        auto xcp = get_xcp_3();
        xcp_profit_real = xcp_profit_real * xcp / this->xcp;
        if (not only_real) {
            xcp_profit = xcp_profit * xcp / this->xcp;
        }
        this->xcp = xcp;
    }

    void update_xcp_2(bool only_real=false) {
        auto xcp = get_xcp_2();
        xcp_profit_real = xcp_profit_real * xcp / this->xcp;
        if (not only_real) {
            xcp_profit = xcp_profit * xcp / this->xcp;
        }
        this->xcp = xcp;
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

    money buy_3(money dx, int i, int j, money max_price=1e100L) {
        //"""
        //Buy y for x
        //"""
        try {
            money x_old[3];
            copy_money_3(x_old, &curve.x[0]);
            auto x = curve.x[i] + dx;
            auto y = curve.y_3(x, i, j);
            auto dy = curve.x[j] - y;
            curve.x[i] = x;
            curve.x[j] = y;
            auto fee = curve.x.size() == 3 ? this->fee_3() : this->fee_2();;
            curve.x[j] += dy * fee;
            dy = dy * (1.L - fee);
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

    money buy_2(money dx, int i, int j, money max_price=1e100L) {
        //"""
        //Buy y for x
        //"""
        try {
            money x_old[2];
            copy_money_2(x_old, &curve.x[0]);
            auto x = curve.x[i] + dx;
            auto y = curve.y_2(x, i, j);
            auto dy = curve.x[j] - y;
            curve.x[i] = x;
            curve.x[j] = y;
            auto fee = this->fee_2();
            curve.x[j] += dy * fee;
            dy = dy * (1.L - fee);
            if ((dx / dy) > max_price or dy < 0) {
                copy_money_2(&curve.x[0], x_old);
                return 0;
            }
            update_xcp_2();
            return dy;
        } catch (...) {
            return 0;
        }
    }

    money sell_3(money dy, int i, int j, money min_price=0) {
        // """
        // Sell y for x
        // """
        try {
            money x_old[3];
            copy_money_3(x_old, &curve.x[0]);
            auto y = curve.x[j] + dy;
            auto x = curve.y_3(y, j, i);
            auto dx = curve.x[i] - x;
            curve.x[i] = x;
            curve.x[j] = y;
            auto fee = this->fee_3();
            curve.x[i] += dx * fee;
            dx = dx * (1.L - fee);
            if ((dx / dy) < min_price or dx < 0) {
                copy_money_3(&curve.x[0], x_old);
                return 0;
            }
            update_xcp_3();
            return dx;
        } catch (...) {
            return 0;
        }
    }

    money sell_2(money dy, int i, int j, money min_price=0) {
        // """
        // Sell y for x
        // """
        try {
            money x_old[2];
            copy_money_2(x_old, &curve.x[0]);
            auto y = curve.x[j] + dy;
            auto x = curve.y_2(y, j, i);
            auto dx = curve.x[i] - x;
            curve.x[i] = x;
            curve.x[j] = y;
            auto fee = this->fee_2();
            curve.x[i] += dx * fee;
            dx = dx * (1.L - fee);
            if ((dx / dy) < min_price or dx < 0) {
                copy_money_2(&curve.x[0], x_old);
                return 0;
            }
            update_xcp_2();
            return dx;
        } catch (...) {
            return 0;
        }
    }

    void ma_recorder(u64 t, vector<money> const &price_vector) {
        //  XXX what if every block only has p_b being last
        if (t > this->t) {
            money alpha = powl(0.5, ((money)(t - this->t) / this->ma_half_time));
            for (int k = 1; k < price_vector.size(); k++) {
                price_oracle[k] = price_vector[k] * (1 - alpha) + price_oracle[k] * alpha;
            }
            this->t = t;
        }
    }

    auto tweak_price_2(u64 t, int a, int b, money p) {
        ma_recorder(t, last_price);
        const int N = 2;
        if (b > 0) {
            last_price[b] = p * last_price[a];
        } else {
            last_price[a] = last_price[0] / p;
        }

        // # price_oracle looks like [1, p1, p2, ...] normalized to 1e18
        money S = 0;
        for (size_t i = 0; i < N; i++) {
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

        money p_new[MAX_ARRAY];
        p_new[0] = 1.L;
        for (size_t i = 1; i < price_oracle.size(); i++) {
            auto p_target = curve.p[i];
            auto p_real = price_oracle[i];
            p_new[i] = p_target + adjustment_step * (p_real - p_target) / norm;
        }
        money old_p[MAX_ARRAY];
        copy_money_2(old_p, &curve.p[0]);

        auto old_profit = xcp_profit_real;
        auto old_xcp = xcp;

        copy_money_2(&curve.p[0],p_new);
        update_xcp_2(true);

        if (2.L * (xcp_profit_real - 1.L) <= (xcp_profit - 1.L)) {
            //  If real profit is less than half of maximum - revert params back
            copy_money_2(&curve.p[0], old_p);
            xcp_profit_real = old_profit;
            xcp = old_xcp;
            not_adjusted = false;
            auto val = ((xcp_profit_real - 1.L - (xcp_profit - 1.L) / 2.L));
            // printf("%.10Lf\n", val);
        }
        return norm;
    }


    auto tweak_price_3(u64 t, int a, int b, money p) {
        ma_recorder(t, last_price);
        const size_t N = 3;
        if (b > 0) {
            last_price[b] = p * last_price[a];
        } else {
            last_price[a] = last_price[0] / p;
        }

        // # price_oracle looks like [1, p1, p2, ...] normalized to 1e18
        money S = 0;
        for (size_t i = 0; i < N; i++) {
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

        money p_new[MAX_ARRAY];
        p_new[0] = 1.L;
        for (size_t i = 1; i < price_oracle.size(); i++) {
            auto p_target = curve.p[i];
            auto p_real = price_oracle[i];
            p_new[i] = p_target + adjustment_step * (p_real - p_target) / norm;
        }
        money old_p[MAX_ARRAY];
        copy_money_3(old_p, &curve.p[0]);

        auto old_profit = xcp_profit_real;
        auto old_xcp = xcp;

        copy_money_3(&curve.p[0],p_new);
        if (N == 3) update_xcp_3(true);
        else        update_xcp_2(true);

        if (2.L * (xcp_profit_real - 1.L) <= (xcp_profit - 1.L)) {
            //  If real profit is less than half of maximum - revert params back
            copy_money_3(&curve.p[0], old_p);
            xcp_profit_real = old_profit;
            xcp = old_xcp;
            not_adjusted = false;
            auto val = ((xcp_profit_real - 1.L - (xcp_profit - 1.L) / 2.L));
            // printf("%.10Lf\n", val);
        }
        return norm;
    }


    void simulate(mapped_file const *in, simulation_data *simdata, extra_data *extdata) {
        // vector<trade_data> const &mdata
        const money CANDLE_VARIATIVES = 50;
        map<pair<int, int>, money> lasts;
        size_t N = price_oracle.size();
        u64 start_t = 0;
        long double last_time = 0;
        size_t total_elements = in->size / sizeof(trade_data);
        simdata->total = total_elements;
        auto mapped_data = (trade_data const *) in->base;
        auto mapped_data_ptr = mapped_data;
        money _slippage;
        for (size_t i = 0; i < total_elements; i++) {
            simdata->current = i;
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
            auto itl = lasts.find({a, b});
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
            auto step1 = N == 3 ? step_for_price_3(candle / CANDLE_VARIATIVES, d.pair1, 1) : step_for_price_2(candle / CANDLE_VARIATIVES, d.pair1, 1);
            auto step2 = N == 3 ? step_for_price_3(candle / CANDLE_VARIATIVES, d.pair1, -1) : step_for_price_2(candle / CANDLE_VARIATIVES, d.pair1, -1);
            auto step = min(step1, step2);
            auto max_price = d.high;
            money _dx = 0;
            auto p_before = N == 3 ? price_3(a, b) : price_2(a, b);
            while (last < (max_price * (1 - ext_fee)) and vol < ext_vol / 2.L) {
                auto dy = N == 3 ? buy_3(step, a, b, max_price) : buy_2(step, a, b, max_price);
                if (dy == 0) {
                    break;
                }
                vol += dy * price_oracle[b];
                _dx += dy;
                last = step / dy; 
                max_price = d.high;
                ctr += 1;
            }
            auto p_after = N == 3 ? price_3(a, b) : price_2(a, b);
            auto _fee = N == 3 ? fee_3() : fee_2();

            if (p_before != p_after) {
                auto v = _dx / curve.x[b];
                _slippage = (_dx * (p_before + p_after)) / (2.L * (mabs(p_before - p_after) + _fee * p_after) * curve.x[b]);
                volume += v;
            }
            if (_slippage > 0) {
                slippage_count += last_time;
                antislippage += last_time * _slippage;
                slippage += last_time / _slippage;
            }

            _high = last;
            auto min_price = d.low;
            _dx = 0;
            p_before = p_after;
            money prev_vol = vol;
            while (last > (min_price * (1 + ext_fee)) and vol < ext_vol / 2.L) {
                auto dx = step / last;
                auto dy = N == 3 ? sell_3(dx, a, b, min_price) : sell_2(dx, a, b, min_price);
                if (dy == 0) {
                    break;
                }
                _dx += dx;
                vol += dx * price_oracle[b];
                last = dy / dx;
                min_price = d.low;
                ctr += 1;
            }
            p_after = N == 3 ? price_3(a, b) : price_2(a, b);
            _fee = N == 3 ? fee_3() : fee_2();

            if (p_before != p_after) {
                auto v = _dx / curve.x[b];
                _slippage = (_dx * (p_before + p_after)) / (2.L * (mabs(p_before - p_after) + _fee * p_after) * curve.x[b]);
                volume += v;
            }
            if (_slippage > 0) {
                slippage_count += last_time;
                antislippage += last_time * _slippage;
                slippage += last_time / _slippage;
            }

            _low = last;
            lasts[d.pair1] = last;

            if (N == 2) tweak_price_2(d.t, a, b, (_high + _low) / 2.L);
            else        tweak_price_3(d.t, a, b, (_high + _low) / 2.L);
            total_vol += vol;
            last_time = d.t;
            long double ARU_x = xcp_profit_real;
            long double ARU_y = (86400.L * 365.L / (d.t - start_t + 1.L));
            APY = powl(ARU_x, ARU_y) - 1.L;
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
                               "AMM: %.5Lf, %0.5Lf\tTarget: %.0Lf, %.0Lf\t"
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
                               (curve.p.size() == 3 ? fee_3() : fee_2()) * 100.L,
                               is_light ? '*' : '.');
                    } else if (N == 2) {
                        printf("t=%llu %.1Lf%%\ttrades: %d\tAMM: %.5Lf\tTarget: %.5Lf\tVol: %.4Lf\tPR:%.2Lf\txCP-growth: {%.5Lf}\tAPY:%.1Lf%%\tfee:%.3Lf%% %c\n",
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
        }
        extdata->slippage = slippage / slippage_count / 2.L;
        extdata->liq_density = 2.L * antislippage / slippage_count;
        extdata->APY = APY;
        extdata->volume = volume;
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
    int log;
    money fee_gamma;
    money total_vol;
    int ma_half_time;
    money ext_fee;
    money volume;
    money slippage;
    money antislippage;
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
        (*(data->result))["configuration"][simdata.num]["Result"]["liq_density"] = simdata.result.liq_density;
        (*(data->result))["configuration"][simdata.num]["Result"]["slippage"] = simdata.result.slippage;
        (*(data->result))["configuration"][simdata.num]["Result"]["volume"] = simdata.result.volume;

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
        if (argv[1][4] == 0) LAST_ELEMS = 100000;
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
