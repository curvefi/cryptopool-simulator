#!/usr/bin/env python3

import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(1e-6), log10(0.5), 64)
Xname = "fee_gamma"
Y = np.logspace(log10(1e-6), log10(0.1), 64)
Yname = "gamma"

other_params = dict(
    D=50e6,
    adjustment_step=1e-7,
    fee_gamma=0.006,
    ma_half_time=600,
    mid_fee=0.003,
    out_fee=0.025,
    gas_fee=5,
    n=3,
    log=0,
    allowed_extra_profit=1e-10,
    ext_fee=0.0003,
    gamma=1e-5,
    boost_rate=0.0,
    A=2)

config = {
    'configuration': [],
    'datafile': [
        "btcusdt-1m-processed",
        "ethusdt-1m-processed",
        "ethbtc-1m-processed"],
    'debug': 0}

for x, y in itertools.product(X, Y):
    params = copy(other_params)
    params[Xname] = x
    params[Yname] = y
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
