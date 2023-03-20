import numpy as np
from scipy.interpolate import interp1d

from simulation_int_many import Curve, geometric_mean


A = int(5000 * 2**2 * 10000)
gamma = int(1e-4 * 1e18)
D = 2 * 10**18
N = 2
dp = 1e-8
dx = 1e-5


class PCurve(Curve):
    def get_xcp(self):
        # First calculate the ideal balance
        # Then calculate, what the constant-product would be
        D = self.D()
        N = len(self.x)
        X = [D * 10**18 // (N * p) for p in self.p]

        return geometric_mean(X)


def profit_deriv_calc(c, p):
    p0 = int((1 + p) * 1e18)
    p1 = int((1 + p + dp) * 1e18)
    c.p[0] = p0
    P0 = c.get_xcp()
    c.p[0] = p1
    P1 = c.get_xcp()
    return (P1 - P0) / P0 * p0 / (p1 - p0)


def f_profit_deriv(A, gamma):
    c = PCurve(A, gamma, D, N)
    d_price = np.linspace(-8e-1, 8e-1, 6000)
    d_profit = [profit_deriv_calc(c, p) for p in d_price]
    return interp1d(d_price, d_profit, fill_value='extrapolate')


def f_liq_density(A, gamma):
    c = PCurve(A, gamma, D, N)
    d_price = np.linspace(-8e-1, 8e-1, 6000)
    dens = []
    for dp in d_price:
        c.p[0] = int(10**18 * (1 + dp))
        xx = np.array([1 - dx, 1, 1 + dx]) * 1e18
        yy = np.array([c.y(int(x), 0, 1) for x in xx])
        d2y = (yy[2] + yy[0] - 2 * yy[1]) / ((xx[2] - xx[0]) / 2)**2
        dens.append(2 / (xx[1] * d2y))
    return interp1d(d_price, np.array(dens), fill_value='extrapolate')


if __name__ == '__main__':
    import pylab
    d_price = np.linspace(-60e-2, 60e-2, 1000)
    # pylab.plot(d_price, f_liq_density(2 * 3**3 * 10000, int(2e-3 * 1e18))(d_price))
    # pylab.plot(d_price, f_profit_deriv(int(2 * 3**3 * 10000), int(2.8e-5 * 1e18))(d_price))
    # pylab.plot(d_price, f_profit_deriv(int(2 * 3**3 * 10000) * 10, int(2.8e-5 * 1e18 / 10**.5))(d_price), '--')
    # pylab.semilogy(d_price, f_liq_density(int(2 * 3**3 * 10000), int(2.8e-5 * 1e18))(d_price))
    pylab.semilogy(d_price, f_liq_density(int(5000 * 2**2 * 10000), int(1e-4 * 1e18))(d_price), '--')
    pylab.show()
