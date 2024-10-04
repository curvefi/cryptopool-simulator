#!/usr/bin/env python3

import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(1), log10(60), 128)
Xname = "A"
Y = np.logspace(log10(1e-2), log10(0.199), 128)
Yname = "gamma"

other_params = dict(
    D=40e3,
    adjustment_step=1e-7,
    fee_gamma=1e-2,
    ma_half_time=600,
    mid_fee=0.0002,
    out_fee=0.003,
    gas_fee=0,
    n=2,
    log=0,
    allowed_extra_profit=1e-12,
    ext_fee=0.0002,
    gamma=1e-5,
    boost_rate=0.0,
    A=200)

config = {
    'configuration': [],
    'datafile': ["eurusd-1m"],
    'debug': 0}

for x, y in itertools.product(X, Y):
    params = copy(other_params)
    params[Xname] = x
    params[Yname] = y
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
