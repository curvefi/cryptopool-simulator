#!/usr/bin/env python3


def geometric_mean(x):
    N = len(x)
    x = sorted(x, reverse=True)  # Presort - good for convergence
    D = x[0]
    for i in range(255):
        D_prev = D
        tmp = 10 ** 18
        for _x in x:
            tmp = tmp * _x // D
        D = D * ((N - 1) * 10**18 + tmp) // (N * 10**18)
        diff = abs(D - D_prev)
        if diff <= 1 or diff * 10**18 < D:
            return D
    print(x)
    raise ValueError("Did not converge")


def reduction_coefficient(x, gamma):
    N = len(x)
    x_prod = 10**18
    K = 10**18
    S = sum(x)
    for x_i in x:
        x_prod = x_prod * x_i // 10**18
        K = K * N * x_i // S
    if gamma > 0:
        K = gamma * 10**18 // (gamma + 10**18 - K)
    return K


def absnewton(f, fprime, x0, handle_x=False, handle_D=False):
    x = x0
    i = 0
    while True:
        x_prev = x
        _f = f(x)
        _fprime = fprime(x)
        x -= _f * 10**18 // _fprime

        # XXX vulnerable to edge-cases
        # Need to take out of unstable equilibrium if ever gets there
        # Might be an issue in smart contracts

        # XXX TODO fuzz if we can remove
        if handle_x:
            if x < 0 or _fprime < 0:
                x = x_prev // 2
        elif handle_D:
            if x < 0:
                x = -x // 2

        i += 1
        # if i > 1000:  # XXX where do we stop?
        #     print(i, (x - x_prev) / x_prev, x, x_prev)
        if i > 1050:
            raise ValueError("Did not converge")
        if abs(x - x_prev) <= max(100, x // 10**14):
            return x


def newton_D(A, gamma, x, D0):
    D = D0
    i = 0

    S = sum(x)
    x = sorted(x, reverse=True)
    N = len(x)
    for j in range(N):  # XXX or just set A to be A*N**N?
        A = A * N
    # XXX
    A = int(A)

    for i in range(255):
        D_prev = D

        K0 = 10**18
        for _x in x:
            K0 = K0 * _x * N // D

        _g1k0 = abs(gamma + 10**18 - K0)

        # D / (A * N**N) * _g1k0**2 / gamma**2
        mul1 = 10**18 * D // gamma * _g1k0 // gamma * _g1k0 // A

        # 2*N*K0 / _g1k0
        mul2 = (2 * 10**18) * N * K0 // _g1k0

        neg_fprime = (S + S * mul2 // 10**18) + mul1 * N // K0 - mul2 * D // 10**18
        assert neg_fprime > 0  # Python only: -f' > 0

        # D -= f / fprime
        D = (D * neg_fprime + D * S - D**2) // neg_fprime - D * (mul1 // neg_fprime) // 10**18 * (10**18 - K0) // K0

        if D < 0:
            D = -D // 2
        if abs(D - D_prev) <= max(100, D // 10**14):
            return D

    raise ValueError("Did not converge")


def newton_y(A, gamma, x, D, i):
    N = len(x)

    y = D // N
    K0_i = 10**18
    S_i = 0
    x_sorted = sorted(_x for j, _x in enumerate(x) if j != i)
    convergence_limit = max(max(x_sorted) // 10**14, D // 10**14, 100)
    for _x in x_sorted:
        y = y * D // (_x * N)  # Small _x first
        S_i += _x
    for _x in x_sorted[::-1]:
        K0_i = K0_i * _x * N // D  # Large _x first

    for j in range(N):  # XXX or just set A to be A*N**N?
        A = A * N
    # XXX
    A = int(A)

    for j in range(255):
        y_prev = y

        K0 = K0_i * y * N // D
        S = S_i + y

        _g1k0 = abs(gamma + 10**18 - K0)

        # D / (A * N**N) * _g1k0**2 / gamma**2
        mul1 = 10**18 * D // gamma * _g1k0 // gamma * _g1k0 // A

        # 2*K0 / _g1k0
        mul2 = 10**18 + (2 * 10**18) * K0 // _g1k0

        yfprime = (10**18 * y + S * mul2 + mul1 - D * mul2)
        fprime = yfprime // y
        assert fprime > 0  # Python only: f' > 0

        # y -= f / f_prime;  y = (y * fprime - f) / fprime
        y = (yfprime + 10**18 * D - 10**18 * S) // fprime + mul1 // fprime * (10**18 - K0) // K0

        if j > 100:  # Just logging when doesn't converge
            print(j, y, D, x)
        if y < 0 or fprime < 0:
            y = y_prev // 2
        if abs(y - y_prev) <= max(convergence_limit, y // 10**14):
            return y

    raise Exception("Did not converge")


def solve_x(A, gamma, x, D, i):
    return newton_y(A, gamma, x, D, i)


def solve_D(A, gamma, x):
    D0 = len(x) * geometric_mean(x)  # <- fuzz to make sure it's ok XXX
    return newton_D(A, gamma, x, D0)


class Curve:
    def __init__(self, A, gamma, D, n, p=None):
        self.A = A
        self.gamma = gamma
        self.n = n
        if p:
            self.p = p
        else:
            self.p = [10 ** 18] * n
        self.x = [D // n * 10**18 // self.p[i] for i in range(n)]

    def xp(self):
        return [x * p // 10 ** 18 for x, p in zip(self.x, self.p)]

    def D(self):
        xp = self.xp()
        if any(x <= 0 for x in xp):
            raise ValueError
        return solve_D(self.A, self.gamma, xp)

    def y(self, x, i, j):
        xp = self.xp()
        xp[i] = x * self.p[i] // 10 ** 18
        yp = solve_x(self.A, self.gamma, xp, self.D(), j)
        return yp * 10**18 // self.p[j]

    def xcp(self):
        D = self.D()
        N = len(self.x)
        X = [D * 10**18 // (N*p) for p in self.p]
        return geometric_mean(X)


if __name__ == '__main__':
    import numpy as np
    import pylab

    p0 = [10**18, 35000 * 10**18, 2000*10**18]
    curve = Curve(A=0.254, gamma=int(8e-4 * 1e18), D=10**18, n=3, p=p0)

    losses = []
    price_shifts = []
    step = 4.9e-4

    for ps in np.linspace(0.8, 1.2, 1000):
        p0 = curve.p[:]
        curve.p[1] = int(curve.p[1] * ps)
        xcp0 = curve.xcp()
        curve.p[1] = int(curve.p[1] * (1 + step))
        xcp = curve.xcp()
        curve.p = p0
        price_shifts.append(ps)
        losses.append(xcp/xcp0 - 1)

    pylab.plot(price_shifts, losses)
    pylab.show()
