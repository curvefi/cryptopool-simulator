#!/usr/bin/env python3

import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(3e-4), log10(0.03), 32)
Xname = "mid_fee"
Y = np.logspace(log10(1e-6), log10(1e-2), 32)
Yname = "fee_gamma"

other_params = dict(
    D=40e6,
    adjustment_step=1e-7,
    fee_gamma=0.003,
    ma_half_time=600,
    mid_fee=0.003,
    out_fee=0.03,
    gas_fee=5,
    n=2,
    log=0,
    allowed_extra_profit=1e-10,
    ext_fee=0.0003,
    gamma=0.0476,
    boost_rate=0.05,
    A=3.14)

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