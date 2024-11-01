#!/usr/bin/env python3

import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(5), log10(200), 64)
Xname = "A"
Y = np.logspace(log10(1e-5), log10(0.2), 64)
Yname = "gamma"

other_params = dict(
    D=40e3,
    adjustment_step=1e-7,
    fee_gamma=3e-3,
    ma_half_time=600,
    mid_fee=0.0003,
    out_fee=0.0020,
    gas_fee=0,
    n=2,
    log=0,
    allowed_extra_profit=1e-10,
    ext_fee=0.0001,
    gamma=2e-3,
    boost_rate=0.05,
    A=100)

config = {
    'configuration': [],
    'datafile': ["xauusd-1m"],
    'debug': 0}

for x, y in itertools.product(X, Y):
    params = copy(other_params)
    params[Xname] = x
    params[Yname] = y
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
