#!/usr/bin/env python3

import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(5), log10(2000), 64)
Xname = "A"
Y = np.logspace(log10(1e-4), log10(0.05), 64)
Yname = "mid_fee"

other_params = dict(
    D=1e6,
    adjustment_step=1e-7,
    fee_gamma=0.003,
    ma_half_time=600,
    mid_fee=0.003,
    out_fee=0.003,
    gas_fee=0,
    n=2,
    log=0,
    allowed_extra_profit=1e-12,
    ext_fee=0.0001,
    gamma=0,
    boost_rate=0.1,
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
