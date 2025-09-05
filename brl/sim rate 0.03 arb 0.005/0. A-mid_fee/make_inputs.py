#!/usr/bin/env python3

import numpy as np
from math import log10
import json
import itertools
from copy import copy

N_GRID = 8
X = np.logspace(log10(5), log10(500), N_GRID)
X = np.linspace(5, 500, N_GRID)
Xname = "A"
Y = np.logspace(log10(1e-4), log10(0.05), N_GRID)
Y = np.linspace(1e-4, 0.05, N_GRID)
Yname = "mid_fee"

other_params = dict(
    D=1e6,
    adjustment_step=1e-7,
    fee_gamma=0.003,
    ma_half_time=600,
    mid_fee=0.003,
    out_fee=0.003,
    gas_fee=1,
    n=2,
    log=0,
    allowed_extra_profit=1e-12,
    ext_fee=0.005,
    # 1 -> 10_000
    # 0.5 -> 5_000
    # 0.05 -> 500
    # 0.005 -> 50
    gamma=0,
    boost_rate=0.03,
    A=10)

config = {
    'configuration': [],
    'datafile': ["brlusd-1m"],
    'debug': 0}

for x, y in itertools.product(X, Y):
    params = copy(other_params)
    params[Xname] = x
    params[Yname] = y
    params['out_fee'] = y
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
