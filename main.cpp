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
//#include "bn_fixed.h"
#include "bn_ldbl.h"
#define DEBUG 1
static int trace = DEBUG;
#define TROFF { trace = 0; }
#define TRON { trace = DEBUG; }
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
#define DP(x)      if (trace) printf("[%d] %s::%s=%.6f ", __LINE__, __FUNCTION__, #x, x)
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
//using i256 = BN<7>; // for 64 bit internals BN
//using i256 = BN<10>; // for 32 bit internals BN
using i256 = BN;
using i256_init = long long;
static const i256       grain("1000000000000000000");
static const long long  igrain = 1000000000000000000;
static const double     dgrain = 1e18;
static const i256       grain13("10000000000000");
static const i256       grain14("100000000000000");
static const i256       grain15("1000000000000000");
static const i256       grain16("10000000000000000");
static const i256       grain17("100000000000000000");
static const i256       i256_0("0");
static const i256       i256_1("1");
static const i256       i256_2("2");
static const i256       i256_10("10");
static const i256       i256_50("50");
static const i256       i256_100("100");


static void print_clock(string const &mesg, clock_t start, clock_t end) {
    printf("%s %.3lf sec\n", mesg.c_str(), double(end - start) / CLOCKS_PER_SEC);
}

struct trade_data {
    u64 t = 0;          // 0
    double open = 0;    // 1
    double high = 0;    // 2
    double low = 0;     // 3
    double close = 0;   // 4
    double volume = 0;  // 5
    pair<int,int> pair1 = {0,0};
    void print() const {
        printf("{ open: %.6f, high: %.6f low: %.6f close: %.6f t: %llu, volume: %.6f pair: (%d,%d)} ",
               this->open, this->high, this->low, this->close, this->t, this->volume, this->pair1.first, this->pair1.second);
    }
};

struct trade_one {
    u64 t;
    trade_data trade;
};

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
    explicit mapped_file(string const &name) {
        fd = open(name.c_str(), O_RDONLY);
        if (fd < 0) return;
        lseek(fd, 0, SEEK_END);
        size = lseek(fd, 0, SEEK_CUR);
        base = (unsigned char *)::mmap(NULL, size, PROT_READ, MAP_SHARED, fd, 0);
        printf("mapped_file::open: base=%p\n", base);
        if (base == MAP_FAILED) {
            perror(name.c_str());
            base = nullptr;
        }
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
    mapped_file mf(name_to_open);
    if (mf.base == nullptr) {
        printf("failed to open '%s'\n", name_to_open.c_str());
        abort();
    }
    vector<trade_data> ret;
    auto p = mf.base;
    if (*p == '[') p++; // skip initial '[';
    auto scan_double = [] (unsigned char *p, double *d) {
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

auto get_all() {
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
        ret.emplace_back(q.trade);
    }
    //printf("total %zu elements\n", ret.size());
    print_clock("sorting took", start_time, end_time);
    return ret;
}

i256 geometric_mean(vector<i256> x) {
    T;
    auto N = x.size();
    sort(x.begin(),x.end(), [](i256 const &l, i256 const &r) {
        return l > r;
    });
    //P256(x);
    auto D = x[0];
    i256 N1256((i256_init)N-1);
    i256 N256((i256_init)N);
    for (int i = 0; i < 255; i++) {
        auto D_prev = D;
        auto tmp = grain;
        for (auto const &_x: x) {
            tmp = tmp * _x / D;
        }
        D = (D * (N1256 * grain + tmp)) / (N256 * grain);
        auto diff = (D - D_prev).abs();
        if (diff <= i256_1 or diff * grain < D) {
            E;
            //P(D); EOL;
            return D;
        }
    }
    throw logic_error("   Did not converge");
}

auto reduction_coefficient(vector<i256> const &x, i256 const &gamma) {
    T;
    auto N256 = i256((i256_init) x.size());
    auto x_prod = grain;
    auto K = grain;
    i256 S;
    for (auto const &q: x) S += q; // = sum(x)
    for (auto const &x_i: x) {
        x_prod = x_prod * x_i / grain;
        K = K * N256 * x_i / S;
    }
    if (gamma.sign() > 0) {
        K = gamma * grain / (gamma + grain - K);
    }
    E;
    //P(K); EOL;
    return K;
}

void print(vector<i256> const &x) {
    printf("[");
    for (size_t i = 0; i < x.size(); i++) {
        printf("%s ", x[i].to_string(10).c_str());
    }
    printf("]\n");
}

auto newton_D(i256 A, i256 const &gamma, vector<i256> x, i256 const &D0) {
    auto D = D0;
    int i = 0;

    i256 S;
    for (auto const &q: x) S += q;
    sort(x.begin(), x.end(), [](i256 const &l, i256 const &r) { return l > r; });
    auto const N = x.size();
    i256 N256((i256_init) N);
    for (size_t j = 0; j < N; j++) { // XXX or just set A to be A*N**N?
        A = A * N256;
    }

    for (int i = 0; i < 255; i++) {
        auto D_prev = D;

        auto K0 = grain;
        for (auto const &_x: x) {
            K0 = K0 * _x * N256 / D;
        }

        auto _g1k0 = (gamma + grain - K0).abs();

        // # D / (A * N**N) * _g1k0**2 / gamma**2
        auto mul1 = grain * D / gamma * _g1k0 / gamma * _g1k0 / A;

        // # 2*N*K0 / _g1k0
        auto mul2 = (grain + grain) * N256 * K0 / _g1k0;

        auto neg_fprime = (S + S * mul2 / grain) + mul1 * N256 / K0 - mul2 * D / grain;
        assert (neg_fprime.sign() > 0); //   # Python only: -f' > 0

        // # D -= f / fprime
        D = (D * neg_fprime + D * S - D * D) / neg_fprime - D * (mul1 / neg_fprime) / grain * (grain - K0) / K0;

        if (D.sign() < 0) {
            D = D.abs() / i256_2;
        }
        if ((D - D_prev).abs() <= max(i256_100, D / grain14)) {
            return D;
        }
    }
    throw logic_error("Did not converge");
}

auto newton_y(i256 A, i256 const &gamma, vector<i256> const &x, i256 const &D, int i) {
    auto save_trace = trace;
    TROFF;
    T;
    //P256(x);
    auto N = x.size();
    i256 N256((i256_init) N);

    auto y = D / N256;
    i256 K0_i(grain);
    i256 S_i;
    vector<i256> x_sorted(N - 1);
    for (int j = 0, cnt = 0; j < N; j++) {
        if (j != i) x_sorted[cnt++] = x[j];
    }
    sort(x_sorted.begin(), x_sorted.end());
    //P256(x_sorted);
    auto max_x_sorted = x_sorted[0];
    // auto max_x_sorted = x_sorted.back();
    for (auto const &q: x_sorted) max_x_sorted = max(max_x_sorted, q);
    //P(max_x_sorted);
    auto convergence_limit = max(max_x_sorted / grain14, D / grain14);
    convergence_limit = max(convergence_limit, i256_100);
    //P(convergence_limit); EOL;
    for (auto const &_x: x_sorted) {
        y = y * D / (_x * N256); //  # Small _x first
        S_i += _x;
    }
    //P(S_i); EOL;
    for (auto it = x_sorted.rbegin(); it != x_sorted.rend(); it++) {
        auto const _x = *it;
        K0_i = K0_i * _x * N256 / D; //  # Large _x first
    }
    //P(K0_i); EOL;
    for (size_t j = 0; j < N; j++) { // in range(N):  # XXX or just set A to be A*N**N?
        A = A * N256;
    }
    //P(A); EOL;
    for (size_t j = 0; j < 255; j++) {
        auto y_prev = y;
        //P(y_prev);

        auto K0 = K0_i * y * N256 / D;
        //P(K0);
        auto S = S_i + y;
        //P(S); EOL;

        auto _g1k0 = (gamma + grain - K0).abs();

        // D / (A * N**N) * _g1k0**2 / gamma**2
        auto mul1 = grain * D / gamma * _g1k0 / gamma * _g1k0 / A;
        //P(mul1); P(D); P(gamma); P(_g1k0); P(A); EOL;
        // 2*K0 / _g1k0
        auto mul2 = grain + (grain+grain) * K0 / _g1k0;
        //P(mul2);

        auto yfprime = (grain * y + S * mul2 + mul1 - D * mul2);
        //P(yfprime);
        auto fprime = yfprime / y;
        //P(fprime); EOL;
        assert (fprime.sign()  > 0) ;  //# Python only: f' > 0

        // y -= f / f_prime;  y = (y * fprime - f) / fprime
        y = (yfprime + grain * D - grain * S) / fprime + mul1 / fprime * (grain - K0) / K0;
        //P(y); EOL;
        if (j > 100) { //  # Just logging when doesn't converge
            printf("%zu %s %s ", j, y.to_string(10).c_str(), D.to_string(10).c_str());
            print(x);
        }
        if (y.sign() < 0 or fprime.sign() < 0) {
            y = y_prev / i256_2;
            //P(y);
        }

        if ((y - y_prev).abs() <= max(convergence_limit, y / grain14)) {
            E;
            //P(y); EOL;
            trace = save_trace;
            TRON;
            return y;
        }
    }
    E;
    throw logic_error("Did not converge");
}

auto solve_x(i256 const &A, i256 const &gamma, vector<i256> const &x, i256 const &D, int i) {
    T;
    //P(A);
    //P(gamma);
    //P256(x);
    //P(D);
    return newton_y(A, gamma, x, D, i);
}

auto solve_D(i256 const &A, i256 const &gamma, vector<i256> const &x) {
    T;
    //P(A);
    //P(gamma);
    //P256(x);
    auto D0 = i256((i256_init)x.size()) * geometric_mean(x); //  # <- fuzz to make sure it's ok XXX
    //P(D0); EOL;
    return newton_D(A, gamma, x, D0);
}

struct Curve {
    Curve(i256 const &A, i256 const &gamma, i256 const &D, int n, vector<i256> const &p) {
        T;
        this->A = A;
        //P(this->A);
        this->gamma = gamma;
        //P(this->gamma);
        this->n = n;
        if (!p.empty()) {
            this->p = p;
        } else {
            this->p.resize(n, grain);
        }
        this->x.resize(n);
        for(int i = 0; i < n; i++) {
            x[i] = D / i256((i256_init)n) * grain / p[i];
        }
        //P256(x);
    }

    auto xp() const {
        T;
        //P256(x);
        //P256(p);
        vector<i256> ret(x.size());
        for (int i = 0; i < x.size(); i++) {
            ret[i] = x[i] * p[i] / grain;
        }
        //P256(ret);
        return ret;
    }

    auto D() const {
        T;
        auto xp = this->xp();
        //P256(xp);
        for (auto const &x: xp) {
            if (x.sign() <= 0) {
                throw logic_error("Curve::D(): x <= 0");
            }
        }
        return solve_D(A, gamma, xp);
    }

    auto y(i256 x, int i, int j) {
        auto xp = this->xp();
        xp[i] = x * this->p[i] / grain;
        auto yp = solve_x(A, gamma, xp, this->D(), j);
        return i256(yp * grain / this->p[j]);
    }

    i256 A;
    i256 gamma;
    int n;
    vector<i256> p;
    vector<i256> x;

};


struct Trader {
    Trader(i256 const &A, i256 const &gamma, i256 const &D, int n, vector<i256> const &p0,
           double mid_fee,
           double out_fee,
           double price_threshold,
           i256 const &fee_gamma,
           double adjustment_step,
           int ma_half_time,
           bool log = true) : curve(A, gamma, D, n, p0) {
        this->p0 = p0;
        this->price_oracle = this->p0;
        this->last_price = this->p0;
        // this->curve = Curve(A, gamma, D, n, p0);
        this->dx = i256(D / i256("100000000"));
        this->mid_fee = i256((mid_fee * dgrain));
        this->out_fee = i256((out_fee * dgrain));
        this->D0 = this->curve.D();
        this->xcp_0 = this->get_xcp();
        this->xcp_profit = grain;
        this->xcp_profit_real = grain;
        this->xcp = this->xcp_0;
        this->price_threshold = i256(price_threshold * dgrain);
        this->adjustment_step = i256(adjustment_step * dgrain);
        this->log = log;
        this->fee_gamma = fee_gamma; // || gamma;
        this->total_vol = 0.0;
        this->ma_half_time = ma_half_time;
        this->ext_fee = i256_0; //   # 0.03e-2
        this->slippage = i256_0;
        this->slippage_count = 0;
        this->not_adjusted = false;
        this->heavy_tx = 0;
        this->light_tx = 0;
        this->is_light = false;
        this->t = 0;
    }


    auto fee() {
        auto f = reduction_coefficient(curve.xp(), fee_gamma);
        return (mid_fee * f + out_fee * (grain - f)) / grain;
    }

    i256 get_xcp() const {
        // First calculate the ideal balance
        //  Then calculate, what the constant-product would be
        auto D = curve.D();
        i256 N256((i256_init) curve.x.size());
        vector<i256> X(curve.p.size());
        for (size_t i = 0; i < X.size(); i++) {
            X[i] = D * grain / (N256 * curve.p[i]);
        }

        return geometric_mean(X);
    }

    auto price(int i, int j) {
        auto dx_raw = dx * grain / curve.p[i];
        auto ret = dx_raw * grain / (curve.x[j] - curve.y(curve.x[i] + dx_raw, i, j));
        //P(ret);
        return ret;
    }

    auto step_for_price(i256 dp, pair<int, int> p, int sign) {
        T;
        //P(dp); // printf("(%d,%d)\n", p.first, p.second);
        auto p0 = price(p.first, p.second);
        //P(p0);
        dp = p0 * dp / grain;
        //P(dp);
        auto x0 = curve.x;
        //P256(x0);
        auto step = dx * grain / curve.p[p.first];
        //P(step);
        while (true) {
            curve.x[p.first] = x0[p.first] + i256((i256_init) sign) * step;
            auto dp_ = (p0 - price(p.first, p.second)).abs();
            if (dp_ >= dp or step >= curve.x[p.first] / i256_10) {
                curve.x = x0;
                E;
                //P(step); EOL;
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

    i256 buy(i256 const &dx, int i, int j, double max_price=1e100) {
        //"""
        //Buy y for x
        //"""
        try {
            auto x_old = curve.x;
            auto x = curve.x[i] + dx;
            auto y = curve.y(x, i, j);
            auto dy = curve.x[j] - y;
            curve.x[i] = x;
            curve.x[j] = y;
            auto fee = this->fee();
            curve.x[j] += dy * fee / grain;
            dy = dy * (grain - fee) / grain;
            if ((dx * grain / dy).get_double() > max_price or dy.sign() < 0) {
                curve.x = x_old;
                return i256_0;
            }
            update_xcp();
            return dy;
        } catch (...) {
            return i256_0;
        }
    }

    i256 sell(i256 const &dy, int i, int j, double min_price=0) {
        // """
        // Sell y for x
        // """
        try {
            auto x_old = curve.x;
            auto y = curve.x[j] + dy;
            auto x = curve.y(y, j, i);
            auto dx = curve.x[i] - x;
            curve.x[i] = x;
            curve.x[j] = y;
            auto fee = this->fee();
            curve.x[i] += dx * fee / grain;
            dx = dx * (grain - fee) / grain;
            if ((dx * grain / dy).get_double() < min_price or dx.sign() < 0) {
                curve.x = x_old;
                return i256_0;
            }
            update_xcp();
            return dx;
        } catch (...) {
            return i256_0;
        }
    }

    void ma_recorder(u64 t, vector<i256> const &price_vector) {
        //  XXX what if every block only has p_b being last
        if (t > this->t) {
            double alpha = pow(0.5, ((double)(t - this->t) / this->ma_half_time));
            for (int k = 1; k <= 2; k++) {
                price_oracle[k] = i256(price_vector[k].get_double() * (1 - alpha) + price_oracle[k].get_double() * alpha);
            }
            this->t = t;
        }
    }

    auto tweak_price(u64 t, int a, int b, i256 const &p) {
        TROFF;
        ma_recorder(t, last_price);
        if (b > 0) {
            last_price[b] = p * last_price[a] / grain;
        } else {
            last_price[a] = last_price[0] * grain / p;
        }

        // # price_oracle looks like [1, p1, p2, ...] normalized to 1e18
        P256(price_oracle);
        P256(curve.p);
        i256 S;
        for (size_t i = 0; i < price_oracle.size(); i++) {
            auto t = price_oracle[i] * grain / curve.p[i] - grain;
            S += t*t;
        }
        auto norm = S;
        // @@@norm.root_to(2);
        auto mxp = (max(price_threshold, adjustment_step));
        norm.root_to();
        P(norm);
        P(price_threshold);
        P(adjustment_step); EOL;
        if (norm <= mxp) {
            // Already close to the target price
            is_light = true;
            light_tx += 1;
            return norm;
        }
        T;
        if (not not_adjusted and (xcp_profit_real - grain > (xcp_profit - grain) / i256_2 + grain13)) {
            not_adjusted = true;
        }
        if (not not_adjusted) {
            light_tx += 1;
            is_light = true;
            return norm;
        }
        heavy_tx += 1;
        is_light = false;

        vector<i256> p_new(price_oracle.size());
        p_new[0] = grain;
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

        if (i256_2 * (xcp_profit_real - grain) <= xcp_profit - grain) {
            //  If real profit is less than half of maximum - revert params back
            curve.p = old_p;
            xcp_profit_real = old_profit;
            xcp = old_xcp;
            not_adjusted = false;
            printf("%s\n", ((xcp_profit_real - grain - (xcp_profit - grain) / i256_2) / grain).to_string(10).c_str());
        }
        return norm;
    }


    void simulate(vector<trade_data> const &mdata) {
        map<pair<int,int>,i256> lasts;
        TROFF;
        auto t = mdata[0].t;
        for (size_t i = 0; i < mdata.size(); i++)  {
            //if (i > 10) abort();
            auto const &d = mdata[i];
            // auto i = d.t;
            auto a = d.pair1.first;
            auto b = d.pair1.second;
            //if (d.t >= 1614973320) exit(0);
            i256 vol;
            auto ext_vol = i256(d.volume * price_oracle[b].get_double()); //  <- now all is in USD
            //P(ext_vol);
            int ctr{0};
            i256 last;
            auto itl = lasts.find({a,b});
            if (itl == lasts.end()) {
                last = price_oracle[b] * grain / price_oracle[a];
            } else {
                last = itl->second;
            }
            //P(last);
            auto _high = last;
            auto _low = last;

            //  Dynamic step
            //  f = reduction_coefficient(self.curve.xp(), self.curve.gamma)
            //DP(d.high);
            //DP(d.low); EOL;
            //P(grain17);
            auto candle = min(i256(dgrain * fabs((d.high - d.low) / d.high)), grain17);
            //P(candle); EOL;
            candle = max(grain15, candle);
            //P(candle); EOL;
            auto step1 = step_for_price(candle / i256_50, d.pair1, 1);
            auto step2 = step_for_price(candle / i256_50, d.pair1, -1);
            auto step = min(step1, step2);
            //P(step1);
            //P(step2);
            //P(step); EOL;

            auto max_price = dgrain * d.high;
            //DP(max_price); EOL;
            i256 _dx;
            auto p_before = price(a, b);
            //P(p_before); EOL;
            while (last.get_double() < max_price and vol < ext_vol / i256_2) {
                auto dy = buy(step, a, b, max_price);
                //P(dy); EOL;
                if (dy.is_zero()) {
                    break;
                }
                vol += dy * price_oracle[b] / grain;
                //P(vol);
                _dx += dy;
                //P(_dx);
                last = step * grain / dy;
                //P(last);
                max_price = dgrain * d.high;
                ctr += 1;
                //DP(max_price); DP((double)ctr); EOL;
            }
            auto p_after = price(a, b);
            //P(p_before); P(p_after); EOL;
            if (p_before != p_after) {
                slippage_count++;
                slippage += _dx * curve.p[b] / grain * (p_before + p_after) / (i256_2 * (p_before - p_after).abs());
                //DP((double)slippage_count); P(slippage); EOL;
            }
            _high = last;
            auto min_price = dgrain * d.low;
            //DP(min_price); EOL;
            _dx = i256_0;
            p_before = p_after;
            //P(p_before); EOL;
            while (last.get_double() > min_price and vol < ext_vol / i256_2) {
                auto dx = step * grain / last;
                //P(dx); EOL;
                auto dy = sell(dx, a, b, min_price);
                //P(dy); EOL;
                _dx += dx;
                if (dy.is_zero()) {
                    break;
                }
                vol += dx * price_oracle[b] / grain;
                last = dy * grain / dx;
                min_price = dgrain * d.low;
                ctr += 1;
                //P(vol); P(last); DP(min_price); DP((double)ctr); EOL;
            }
            p_after = price(a, b);
            //P(p_before); P(p_after); EOL;
            if (p_before != p_after) {
                slippage_count += 1;
                slippage += _dx * curve.p[b] / grain * (p_before + p_after) / (i256_2 * (p_before - p_after).abs());
            }
            _low = last;
            //P(_low); EOL;
            lasts[d.pair1] = last;

            tweak_price(d.t, a, b, (_high + _low) / i256_2);

            total_vol += vol.get_double();
            //DP(total_vol); EOL;
            if (log) {
                try {
                    double last01, last02;
                    auto it01 = lasts.find({0,1});
                    if (it01 == lasts.end()) {
                        last01 = price_oracle[1].get_double() * dgrain / price_oracle[0].get_double() / dgrain;
                    } else {
                        last01 = it01->second.get_double() / dgrain;
                    }
                    auto it02 = lasts.find({0,2});
                    if (it02 == lasts.end()) {
                        last02 = price_oracle[2].get_double() * dgrain / price_oracle[0].get_double() / dgrain;
                    } else {
                        last02 = it02->second.get_double() / dgrain;
                    }
                    printf("t=%llu %.1f%%\ttrades: %d\t"
                           "AMM: %.0f, %0.f\tTarget: %.0f, %.0f\t"
                           "Vol: %.4f\tPR:%.2f\txCP-growth: {%.5f}\t"
                           "APY:%.1f%%\tfee:%.3f%% %c\n",
                          d.t,
                          100. * i / mdata.size(), ctr, last01, last02,
                          curve.p[1].get_double() / dgrain,
                          curve.p[2].get_double() / dgrain,
                          total_vol / dgrain,
                          (xcp_profit_real - grain).get_double() / (xcp_profit - grain).get_double(),
                          xcp_profit_real.get_double() / dgrain,
                          pow((xcp_profit_real.get_double() / dgrain),(86400. * 365. / (d.t - mdata[0].t + 1)) - 1) * 100.,
                          fee().get_double() / dgrain * 100.,
                          is_light? '*' : '.');
                } catch (exception const &e) {
                    printf("caught '%s'\n", e.what());
                }
            }
        }
    }



    vector<i256> p0;
    vector<i256> price_oracle;
    vector<i256> last_price;
    u64 t;
    i256 dx;
    i256 mid_fee;
    i256 out_fee;
    i256 D0;
    i256 xcp, xcp_0;
    i256 xcp_profit;
    i256 xcp_profit_real;
    i256 price_threshold;
    i256 adjustment_step;
    bool log;
    i256 fee_gamma;
    double total_vol;
    int ma_half_time;
    i256 ext_fee;
    i256 slippage;
    int slippage_count;
    bool not_adjusted;
    int  heavy_tx;
    int  light_tx;
    bool is_light;
    Curve curve;


};

auto get_price_vector(int n, vector<trade_data> const &data) {
    vector<i256> p(n);
    p[0] = grain;
    for (auto const &d: data) {
        if (d.pair1.first == 0) {
            p[d.pair1.second] = i256(dgrain * d.close);
        }
        bool zeros = false;
        for (auto x: p) {
            zeros |= x.is_zero();
        }
        if (!zeros) {
            for (auto q: p) {
                printf("%s ", q.to_string(10).c_str());
            }
            printf("\n");
            return p;
        }
    }
    return p;
}


int main() {
    const int LAST_ELEMS = 100000;
    clock_t start = clock();
    auto test_data = get_all();
    test_data.erase(test_data.begin(), test_data.begin() + test_data.size() - LAST_ELEMS);
    //debug_print("test_data first 5", test_data, 5);
    //debug_print("test_data last 5", test_data, -5);
    Trader trader(i256("135"), i256((7e-5 * dgrain)), i256(grain*i256((i256_init )5'000'000)), 3, get_price_vector(3, test_data),
            4e-4, 4.0e-3,
            0.0028, i256(0.01 * dgrain),
            0.0015, 600);
    trader.simulate(test_data);
    printf("Fraction of light transactions:%.5f\n", (double)(trader.light_tx) / (trader.light_tx + trader.heavy_tx));
    clock_t end = clock();
    print_clock("Total time", start, end);
}
