#!/usr/bin/env python3

import numpy as np
from math import log10
import json
import itertools
from copy import copy

A_0 = 3.62548
gamma_0 = 0.0278
fee_gamma_0 = 0.01

X = np.logspace(log10(0.5), log10(40), 32)
Xname = "A"
Y = np.logspace(log10(1e-3), log10(5e-1), 32)
Yname = "gamma"

other_params = dict(
    D=40e6,
    adjustment_step=1e-7,
    fee_gamma=0.005,
    ma_half_time=600,
    mid_fee=0.003,
    out_fee=0.03,
    gas_fee=5,
    n=2,
    log=0,
    allowed_extra_profit=1e-10,
    ext_fee=0.0003,
    gamma=1e-5,
    boost_rate=0.0,
    A=10)

config = {
    'configuration': [],
    'datafile': ["btcusdt-2023-2024"],
    'debug': 0}

for x, y in itertools.product(X, Y):
    params = copy(other_params)
    params[Xname] = x
    params[Yname] = y
    params['fee_gamma'] = fee_gamma_0 * (A_0 * gamma_0)**2 / (params['A'] * params['gamma'])**2
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
