#!/usr/bin/env python3

import numpy as np
from math import log10
import json
from copy import copy


X = np.logspace(log10(0.0003), log10(0.95), 128)
param = "mid_fee"

other_params = dict(
    D=40e6,
    adjustment_step=1e-7,
    fee_gamma=0.0005,
    ma_half_time=600,
    mid_fee=0.003,
    out_fee=0.003,
    gas_fee=10,
    n=2,
    log=0,
    allowed_extra_profit=1e-10,
    ext_fee=0.0003,
    gamma=0.06,
    boost_rate=0.05,
    A=2.6)

config = {
    'configuration': [],
    'datafile': ["btcusdt-2023-2024"],
    'debug': 0}

for x in X:
    params = copy(other_params)
    params[param] = x
    params['out_fee'] = x
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)