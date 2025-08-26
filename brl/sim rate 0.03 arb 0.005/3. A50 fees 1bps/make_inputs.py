#!/usr/bin/env python3

import numpy as np
from math import log10
import json
from copy import copy


X = np.logspace(log10(1e-3), log10(0.5), 256)
Xname = "boost_rate"

other_params = dict(
    D=1e6,
    adjustment_step=1e-7,
    fee_gamma=0.003,
    ma_half_time=600,
    mid_fee=0.0001,
    out_fee=0.0001,
    gas_fee=0,
    n=2,
    log=0,
    allowed_extra_profit=1e-12,
    ext_fee=0.005,
    gamma=0,
    boost_rate=0.03,
    A=50)

config = {
    'configuration': [],
    'datafile': ["brlusd-1m"],
    'debug': 0}

for x in X:
    params = copy(other_params)
    params[Xname] = x
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
