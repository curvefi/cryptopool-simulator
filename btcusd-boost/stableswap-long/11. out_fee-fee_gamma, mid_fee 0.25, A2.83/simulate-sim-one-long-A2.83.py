#!/usr/bin/env python3

import json
import lzma
import numpy as np
from math import sqrt
from datetime import datetime


PEG_TO = 'price_oracle'  # price_oracle vs price_scale
BORROW_RATE = 0.1


class AMM:
    def __init__(self, collateral, leverage, fee, oracle):
        self.collateral = collateral  # y
        self.leverage = leverage
        self.fee = fee
        self.p_oracle = oracle
        self.debt = collateral * oracle * (leverage - 1) / leverage

    def set_p_oracle(self, p):
        self.p_oracle = p

    def get_x0(self):
        lev_ratio = (self.leverage / (2 * self.leverage - 1))**2
        D = self.p_oracle**2 * self.collateral**2 - 4 * self.p_oracle * self.collateral * self.debt * lev_ratio
        return (self.p_oracle * self.collateral + sqrt(D)) / (2 * lev_ratio)

    def get_invariant(self):
        x0 = self.get_x0()
        return self.collateral * (x0 - self.debt)

    def get_dy(self, i, j, in_amount):
        x0 = self.get_x0()
        y_initial = self.collateral
        x_initial = x0 - self.debt
        if i == 0:
            x = x_initial + in_amount
            y = x_initial * y_initial / x
            return (y_initial - y) * (1 - self.fee)
        else:
            y = y_initial + in_amount
            x = x_initial * y_initial / y
            return (x_initial - x) * (1 - self.fee)

    def exchange(self, i, j, amount):
        x0 = self.get_x0()
        y_initial = self.collateral
        x_initial = x0 - self.debt
        if i == 0:
            x = x_initial + amount
            y = x_initial * y_initial / x
            dy = (y_initial - y) * (1 - self.fee)
            self.debt -= amount
            self.collateral -= dy
            return dy
        else:
            y = y_initial + amount
            x = x_initial * y_initial / y
            dx = (x_initial - x) * (1 - self.fee)
            self.collateral += amount
            self.debt += dx
            return dx

    def get_p(self):
        x0 = self.get_x0()
        return (x0 - self.debt) / self.collateral

    def trade_to_price(self, p):
        x0 = self.get_x0()
        initial_price = (x0 - self.debt) / self.collateral

        if p > initial_price:
            # Decrease debt + decrease collateral
            p *= (1 - self.fee)
            if p <= initial_price:
                return

        elif p < initial_price:
            p /= (1 - self.fee)
            if p >= initial_price:
                return

        else:
            return

        inv = (x0 - self.debt) * self.collateral
        x_after = sqrt(inv * p)
        y_after = x_after / p

        if p > initial_price:
            x_after += (x_after - (x0 - self.debt)) * self.fee  # Buy collateral -> more USD in
        else:
            y_after += (y_after - self.collateral) * self.fee  # Dupm more collateral -> more collateral in

        self.collateral = y_after
        self.debt = x0 - x_after

    def get_value(self):
        x0 = self.get_x0()
        Ip = ((x0 - self.debt) * self.collateral * self.p_oracle) ** 0.5
        return 2 * Ip - x0


class Simulator:
    def __init__(self, filename, ext_fee,
                 log=False, verbose=False):
        """
        filename - OHLC data in the same format as Binance returns
        ext_fee - Fee which arb trader pays to external platforms
        """
        self.filename = filename
        self.ext_fee = ext_fee
        self.log = log
        self.verbose = verbose
        self.ema_time = 0
        self.emas = []
        self.load_prices()

    def load_prices(self):
        with lzma.open(self.filename, 'r') as f:
            self.simulation_data = json.load(f)

    def single_run(self, fee, leverage):
        # Data for prices
        initial_price = 1.0
        amm = AMM(collateral=(1.0 * leverage), leverage=leverage, fee=fee, oracle=initial_price)
        initial_value = amm.get_value() / initial_price**leverage

        losses = []

        t_start = self.simulation_data[0]['t']
        t_end = self.simulation_data[-1]['t']

        ema0 = self.simulation_data[0][PEG_TO]
        V0 = self.simulation_data[0]['token0'] + self.simulation_data[0]['token1'] * self.simulation_data[0]['low']

        t_prev = t_start

        for d in self.simulation_data:
            t = d['t']
            o = d['open']
            high = d['high']
            low = d['low']
            pool_profit = 1 + d['profit']
            ema = (d[PEG_TO] / ema0)**0.5
            amm.set_p_oracle(ema * pool_profit)

            r = (high / low)**0.5
            low = (d['token0'] + d['token1'] * low) / V0
            high = low * r
            high = high * (1 - self.ext_fee)
            low = low * (1 + self.ext_fee)

            if high > amm.get_p():
                amm.trade_to_price(high)

            if low < amm.get_p():
                amm.trade_to_price(low)

            amm.debt *= (1 + BORROW_RATE * (t - t_prev) / (86400 * 365))
            t_prev = t

            if self.log or self.verbose:
                # current_value = amm.get_value() / ((d['close'] / self.simulation_data[0]['open'])**0.5)**leverage
                current_value = amm.get_value() / ema**leverage
                d = datetime.fromtimestamp(t).strftime("%Y/%m/%d %H:%M")
                loss = current_value / initial_value * 100
                if self.log:
                    print(f'{d}\t{o:.2f}\t{ema:.2f}\t{amm.get_p():.2f}\t\t{loss:.2f}%')
                if self.verbose:
                    losses.append([t, loss / 100])

        if losses:
            self.losses = losses

        # current_value = amm.get_value() / ((self.simulation_data[-1]['close'] / self.simulation_data[0]['open'])**0.5)**leverage
        current_value = amm.get_value() / ema**leverage
        loss = (current_value / initial_value) ** ((365 * 86400) / (t_end - t_start)) - 1
        return loss


if __name__ == '__main__':
    simulator = Simulator(
            filename='detailed-output-A2.83-fee0.25.json.xz',
            ext_fee=0.0002,
            log=False, verbose=True)

    simulator.single_run(fee=7e-3, leverage=2)

    import pylab
    losses = np.array(simulator.losses[::100])
    t = losses[:, 0]
    t = (t - t[0]) / 86400
    loss = losses[:, 1] * 100 - 100
    pylab.plot(t, loss)
    pylab.xlabel('t (days)')
    pylab.ylabel('Deposit growth (%)')
    pylab.tight_layout()
    pylab.show()

    # Optimal fee = 7e-3
