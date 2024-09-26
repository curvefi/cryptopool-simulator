#!/usr/bin/env python3

import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(0.5), log10(30), 32)
Xname = "A"
Y = np.logspace(log10(1e-6), log10(1.99e-1), 32)
Yname = "gamma"

other_params = dict(
    D=40e6,
    adjustment_step=1e-7,
    fee_gamma=0.01,
    ma_half_time=600,
    mid_fee=0.007,
    out_fee=0.03,
    gas_fee=5,
    n=2,
    log=0,
    allowed_extra_profit=1e-10,
    ext_fee=0.0003,
    gamma=1e-5,
    boost_rate=0.05,
    A=10)

config = {
    'configuration': [],
    'datafile': ["btcusdt-2023-2024"],
    'debug': 0}

for x, y in itertools.product(X, Y):
    params = copy(other_params)
    params[Xname] = x
    params[Yname] = y
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
