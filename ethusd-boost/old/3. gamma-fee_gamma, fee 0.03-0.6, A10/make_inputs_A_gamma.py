#!/usr/bin/env python3

import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(1e-5), log10(4e-2), 32)
Xname = "fee_gamma"
Y = np.logspace(log10(2e-5), log10(4e-3), 32)
Yname = "gamma"

other_params = dict(
    D=40e6,
    adjustment_step=1e-7,
    fee_gamma=0.006,
    ma_half_time=600,
    mid_fee=0.0003,
    out_fee=0.006,
    gas_fee=10,
    n=2,
    log=0,
    allowed_extra_profit=1e-10,
    ext_fee=0.0003,
    gamma=1e-5,
    boost_rate=0.05,
    A=10)

config = {
    'configuration': [],
    'datafile': ["ethusdt-1m"],
    'debug': 0}

for x, y in itertools.product(X, Y):
    params = copy(other_params)
    params[Xname] = x
    params[Yname] = y
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
