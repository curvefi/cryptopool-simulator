#!/usr/bin/env python3

import numpy as np
from math import log10
import json
import itertools
from copy import copy

A_0 = 7
gamma_0 = 2e-5
fee_gamma_0 = 3.7e-5

X = np.logspace(log10(0.3), log10(100), 64)
Xname = "A"
Y = np.logspace(log10(5e-7), log10(5e-3), 64)
Yname = "gamma"

other_params = dict(
    D=40e3,
    adjustment_step=1e-7,
    fee_gamma=fee_gamma_0,
    ma_half_time=600,
    mid_fee=0.0003,
    out_fee=0.002,
    gas_fee=0,
    n=2,
    log=0,
    allowed_extra_profit=1e-10,
    ext_fee=0.0001,
    gamma=2e-3,
    boost_rate=0.0,
    A=100)

config = {
    'configuration': [],
    'datafile': ["cnyusd-1m"],
    'debug': 0}

for x, y in itertools.product(X, Y):
    params = copy(other_params)
    params[Xname] = x
    params[Yname] = y
    # params['fee_gamma'] = fee_gamma_0 * (params['A'] * params['gamma'])**2 / (A_0 * gamma_0)**2
    # params['fee_gamma'] = fee_gamma_0 * (A_0 * gamma_0)**2 / (params['A'] * params['gamma'])**2
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
